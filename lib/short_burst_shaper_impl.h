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

#ifndef INCLUDED_GFDM_SHORT_BURST_SHAPER_IMPL_H
#define INCLUDED_GFDM_SHORT_BURST_SHAPER_IMPL_H

#include <gfdm/short_burst_shaper.h>

static const pmt::pmt_t MSG_IN_PORT = pmt::string_to_symbol("time_tag");
static const pmt::pmt_t MSG_OUT_PORT = pmt::string_to_symbol("command");

namespace gr {
namespace gfdm {

class short_burst_shaper_impl : public short_burst_shaper
{
private:
    const int d_pre_padding;
    const int d_post_padding;
    gr_complex d_scale;
    const bool d_use_timed_commands;
    bool d_has_new_time_tag = false;

    double d_timing_advance;
    uint64_t d_timing_advance_ticks;

    double d_cycle_interval;
    uint64_t d_cycle_interval_ticks;

    uint64_t d_full_secs = 0;
    double d_frac_secs = 0.0;
    uint64_t d_tag_offset = 0;
    double d_samp_rate = 0.0;
    uint64_t d_slot_counter = 0;

    uint64_t d_last_full_secs = 0;
    double d_last_frac_secs = 0.0;

    uint64_t d_last_tx_ns = 0;

    uint64_t d_rx_time = 0;


    uint64_t double2ticks(const double interval) const
    {
        return uint64_t(1000000000ull * interval);
    }

    uint64_t timespec2ticks(const uint64_t full_secs, const double frac_secs) const
    {
        return 1000000000ull * full_secs + double2ticks(frac_secs);
    }

    uint64_t ticks2fullsecs(const uint64_t ticks) const { return ticks / 1000000000ull; }

    double ticks2fracsecs(const uint64_t ticks) const
    {
        return double(ticks % 1000000000ull) / 1000000000.0d;
    }

    uint64_t pc_clock_ticks() const
    {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(
                   std::chrono::high_resolution_clock::now().time_since_epoch())
            .count();
    }

    const pmt::pmt_t d_tx_time_key = pmt::intern("tx_time");
    const pmt::pmt_t d_command_time_key = pmt::mp("time");
    const pmt::pmt_t d_command_gain_key = pmt::mp("gain");

    void send_rx_gain_commands(const int packet_len);
    void send_rx_gain_command(const uint64_t full_secs,
                              const double frac_secs,
                              const double gain);

protected:
    int calculate_output_stream_length(const gr_vector_int& ninput_items);

public:
    short_burst_shaper_impl(int pre_padding,
                            int post_padding,
                            gr_complex scale,
                            const std::string& length_tag_name,
                            bool use_timed_commands,
                            double timing_advance,
                            double cycle_interval);
    ~short_burst_shaper_impl();

    gr_complex scale() const { return d_scale; }
    void set_scale(gr_complex scale) { d_scale = scale; }

    double timing_advance() const { return d_timing_advance; }
    void set_timing_advance(double timing_advance)
    {
        d_timing_advance = timing_advance;
        d_timing_advance_ticks = double2ticks(timing_advance);
    }

    double cycle_interval() const { return d_cycle_interval; }
    void set_cycle_interval(double cycle_interval)
    {
        d_cycle_interval = cycle_interval;
        d_cycle_interval_ticks = double2ticks(cycle_interval);
    }

    void handle_msg(pmt::pmt_t time_msg);

    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_int& ninput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_SHORT_BURST_SHAPER_IMPL_H */
