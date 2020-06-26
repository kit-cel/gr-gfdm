/* -*- c++ -*- */
/*
 * Copyright 2019 Johannes Demel.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "short_burst_shaper_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <chrono>
#include <cstring>
#include <exception>

namespace gr {
namespace gfdm {

short_burst_shaper::sptr short_burst_shaper::make(int pre_padding,
                                                  int post_padding,
                                                  gr_complex scale,
                                                  const std::string& length_tag_name,
                                                  bool use_timed_commands,
                                                  double timing_advance,
                                                  double cycle_interval)
{
    return gnuradio::get_initial_sptr(new short_burst_shaper_impl(pre_padding,
                                                                  post_padding,
                                                                  scale,
                                                                  length_tag_name,
                                                                  use_timed_commands,
                                                                  timing_advance,
                                                                  cycle_interval));
}

/*
 * The private constructor
 */
short_burst_shaper_impl::short_burst_shaper_impl(int pre_padding,
                                                 int post_padding,
                                                 gr_complex scale,
                                                 const std::string& length_tag_name,
                                                 bool use_timed_commands,
                                                 double timing_advance,
                                                 double cycle_interval)
    : gr::tagged_stream_block("short_burst_shaper",
                              gr::io_signature::make(1, 1, sizeof(gr_complex)),
                              gr::io_signature::make(1, 1, sizeof(gr_complex)),
                              length_tag_name),
      d_pre_padding(pre_padding),
      d_post_padding(post_padding),
      d_scale(scale),
      d_use_timed_commands(use_timed_commands),
      d_timing_advance(timing_advance),
      d_timing_advance_ticks(double2ticks(timing_advance)),
      d_cycle_interval(cycle_interval),
      d_cycle_interval_ticks(double2ticks(cycle_interval))
{
    if (d_pre_padding < 0) {
        throw std::invalid_argument("Pre-padding length MUST be >= 0!");
    }
    if (d_post_padding < 0) {
        throw std::invalid_argument("Post-padding length MUST be >= 0!");
    }

    GR_LOG_DEBUG(d_logger,
                 "cycle interval: " + std::to_string(d_cycle_interval) +
                     ", ticks: " + std::to_string(d_cycle_interval_ticks));

    enable_update_rate(true);

    message_port_register_out(MSG_OUT_PORT);

    message_port_register_in(MSG_IN_PORT);
    set_msg_handler(MSG_IN_PORT, [this](pmt::pmt_t msg) { this->handle_msg(msg); });
}

/*
 * Our virtual destructor.
 */
short_burst_shaper_impl::~short_burst_shaper_impl() {}

int short_burst_shaper_impl::calculate_output_stream_length(
    const gr_vector_int& ninput_items)
{
    int noutput_items = ninput_items[0] + d_pre_padding + d_post_padding;
    return noutput_items;
}

void short_burst_shaper_impl::send_rx_gain_command(const uint64_t full_secs,
                                                   const double frac_secs,
                                                   const double gain)
{
    auto msg = pmt::dict_add(
        pmt::make_dict(),
        d_command_time_key,
        pmt::cons(pmt::from_uint64(full_secs), pmt::from_double(frac_secs)));
    msg = pmt::dict_add(msg, d_command_gain_key, pmt::from_double(gain));

    message_port_pub(MSG_OUT_PORT, msg);
}

void short_burst_shaper_impl::send_rx_gain_commands(const int packet_len)
{
    double frac_secs = d_frac_secs - 1.0e-4;
    uint64_t full_secs = d_full_secs;
    if (frac_secs < 0.0) {
        full_secs -= 1;
        frac_secs += 1.0;
    }
    send_rx_gain_command(full_secs, frac_secs, 0.0);

    frac_secs = d_frac_secs + 1.0e-4 + packet_len / d_samp_rate;
    full_secs = d_full_secs;
    if (frac_secs >= 1.0) {
        full_secs += 1;
        frac_secs -= 1.0;
    }
    send_rx_gain_command(full_secs, frac_secs, 65.0);
}

void short_burst_shaper_impl::handle_msg(pmt::pmt_t time_msg)
{
    if (pmt::is_dict(time_msg) &&
        (pmt::equal(pmt::dict_ref(time_msg, pmt::mp("src"), pmt::from_long(-1)),
                    pmt::from_long(0)))) {
        d_rx_time =
            pmt::to_long(pmt::dict_ref(time_msg, pmt::mp("rx_time"), pmt::from_long(0)));

    } else if (d_use_timed_commands && pmt::is_tuple(time_msg)) {
        // if (not pmt::is_tuple(time_msg)) {
        //     return; // time messages are tuples!
        // }
        d_full_secs = pmt::to_uint64(pmt::tuple_ref(time_msg, 0));
        d_frac_secs = pmt::to_double(pmt::tuple_ref(time_msg, 1));
        d_tag_offset = pmt::to_uint64(pmt::tuple_ref(time_msg, 2));
        d_samp_rate = pmt::to_double(pmt::tuple_ref(time_msg, 3));
        d_slot_counter = pmt::to_uint64(pmt::tuple_ref(time_msg, 4));

        d_has_new_time_tag = true;
    }
}

int short_burst_shaper_impl::work(int noutput_items,
                                  gr_vector_int& ninput_items,
                                  gr_vector_const_void_star& input_items,
                                  gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    // if (d_use_timed_commands && !d_has_new_time_tag) {
    //     return 0; // Do not emit packet without timing info!
    // }

    std::memset(out, 0, sizeof(gr_complex) * d_pre_padding);

    volk_32fc_s32fc_multiply_32fc(out + d_pre_padding, in, d_scale, ninput_items[0]);

    std::memset(
        out + d_pre_padding + ninput_items[0], 0, sizeof(gr_complex) * d_post_padding);

    if (d_use_timed_commands) {

        uint64_t fts = pc_clock_ticks();
        uint64_t ticks = fts;

        fts -= fts % d_cycle_interval_ticks;
        fts += d_cycle_interval_ticks;
        while (fts <= d_last_tx_ns) {
            fts += d_cycle_interval_ticks;
        }
        fts += d_rx_time % d_cycle_interval_ticks;
        d_last_tx_ns = fts;

        fts += d_timing_advance_ticks;

        uint64_t full_secs = ticks2fullsecs(fts);
        double frac_secs = ticks2fracsecs(fts);

        // uint64_t heartbeat_counter = timespec2ticks(d_full_secs, d_frac_secs);
        // uint64_t tx_counter = timespec2ticks(full_secs, frac_secs);
        // if (heartbeat_counter > tx_counter) {
        //     GR_LOG_DEBUG(d_logger,
        //                  "timestamp: " + std::to_string(d_full_secs) + " . " +
        //                      std::to_string(1000.0 * d_frac_secs) +
        //                      " / system: " + std::to_string(full_secs) + " . " +
        //                      std::to_string(1000.0 * frac_secs));
        // }


        add_item_tag(
            0,
            nitems_written(0),
            d_tx_time_key,
            pmt::make_tuple(pmt::from_uint64(full_secs), pmt::from_double(frac_secs)));

        // send_rx_gain_command(full_secs, frac_secs, 0.0f);
        // const uint64_t eob_ticks = fts + noutput_items + d_pre_padding +
        // d_post_padding; send_rx_gain_command(ticks2fullsecs(eob_ticks),
        // ticks2fracsecs(eob_ticks), 65.0f);

        d_has_new_time_tag = false;
        d_last_full_secs = full_secs;
        d_last_frac_secs = frac_secs;

        std::vector<tag_t> tags;
        get_tags_in_range(
            tags, 0, nitems_read(0), (nitems_read(0) + noutput_items), pmt::mp("time"));
        uint64_t tx_timestamp = 0;
        for (auto t : tags) {
            tx_timestamp = pmt::to_long(t.value);
        }

        // GR_LOG_DEBUG(d_logger, "Timestamp: " + std::to_string(tx_timestamp) + " PC
        // timestamp: " + std::to_string(ticks) + " TX timestamp: " +
        // std::to_string(fts)); GR_LOG_DEBUG(d_logger, "TX timestamp: " +
        // std::to_string(fts));
    }

    // Tell runtime system how many output items we produced.
    return ninput_items[0] + d_pre_padding + d_post_padding;
}

} /* namespace gfdm */
} /* namespace gr */
