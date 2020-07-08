/* -*- c++ -*- */
/*
 * Copyright 2017 Johannes Demel.
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

#include "extract_burst_cc_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <cmath>

namespace gr {
namespace gfdm {

extract_burst_cc::sptr extract_burst_cc::make(int burst_len,
                                              int tag_backoff,
                                              std::string burst_start_tag,
                                              bool activate_cfo_correction)
{
    return gnuradio::get_initial_sptr(new extract_burst_cc_impl(
        burst_len, tag_backoff, burst_start_tag, activate_cfo_correction));
}

/*
 * The private constructor
 */
extract_burst_cc_impl::extract_burst_cc_impl(int burst_len,
                                             int tag_backoff,
                                             std::string burst_start_tag,
                                             bool activate_cfo_correction)
    : gr::block("extract_burst_cc",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_burst_len(burst_len),
      d_tag_backoff(tag_backoff),
      d_burst_start_tag(pmt::mp(burst_start_tag)),
      d_activate_cfo_correction(activate_cfo_correction)
{
    set_output_multiple(burst_len);
    set_tag_propagation_policy(TPP_DONT);
}

/*
 * Our virtual destructor.
 */
extract_burst_cc_impl::~extract_burst_cc_impl() {}

void extract_burst_cc_impl::forecast(int noutput_items,
                                     gr_vector_int& ninput_items_required)
{
    ninput_items_required[0] = noutput_items;
}

float extract_burst_cc_impl::get_scale_factor(const pmt::pmt_t& info) const
{
    float scale_factor = 1.0f;
    if (pmt::is_dict(info)) {
        scale_factor = pmt::to_double(
            pmt::dict_ref(info, d_scale_factor_key, pmt::from_double(1.0)));
    }
    return scale_factor;
}

gr_complex extract_burst_cc_impl::get_phase_rotation(const pmt::pmt_t& info) const
{
    const auto phase_rotation = pmt::to_complex(pmt::dict_ref(
        info, d_phase_rotation_key, pmt::from_complex(gr_complex(1.0f, 0.0f))));
    const auto scale = 1.0 / std::abs(phase_rotation);
    return gr_complex(scale * phase_rotation.real(),
                      -1.0f * scale * phase_rotation.imag());
}

void extract_burst_cc_impl::normalize_power_level(gr_complex* p_out,
                                                  const gr_complex* p_in,
                                                  const float norm_factor,
                                                  const int ninput_size)
{
    volk_32f_s32f_multiply_32f(
        (float*)p_out, (const float*)p_in, norm_factor, 2 * ninput_size);
}

void extract_burst_cc_impl::activate_cfo_compensation(bool activate_cfo_compensation)
{
    std::cout << "activate_cfo_compensation="
              << (activate_cfo_compensation ? "True" : "False") << std::endl;
    d_activate_cfo_correction = activate_cfo_compensation;
}

void extract_burst_cc_impl::compensate_cfo(gr_complex* p_out,
                                           const gr_complex* p_in,
                                           const gr_complex phase_increment,
                                           const int ninput_size)
{
    gr_complex initial_phase = gr_complex(1.0f, 0.0f);
    volk_32fc_s32fc_x2_rotator_32fc(
        p_out, p_in, phase_increment, &initial_phase, ninput_size);
}

int extract_burst_cc_impl::general_work(int noutput_items,
                                        gr_vector_int& ninput_items,
                                        gr_vector_const_void_star& input_items,
                                        gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    const int n_out_bursts = noutput_items / d_burst_len;
    const int avail_items = ninput_items[0];
    int consumed_items = avail_items;
    int produced_items = 0;


    std::vector<tag_t> tags;
    get_tags_in_window(tags, 0, 0, avail_items, d_burst_start_tag);
    std::sort(tags.begin(), tags.end(), tag_t::offset_compare);
    const int n_max_bursts = std::min(int(tags.size()), n_out_bursts);
    // if (tags.size() > 0) {
    //     std::string offset_str("");
    //     for (const auto& t : tags) {
    //         offset_str += std::to_string(t.offset) + "\t";
    //     }
    //     GR_LOG_DEBUG(d_logger,
    //                  "Tags: " + std::to_string(tags.size()) + "/" +
    //                      std::to_string(n_out_bursts) +
    //                      " nitems_read=" + std::to_string(nitems_read(0)) +
    //                      "\tavail=" + std::to_string(avail_items) + "\tnout_items=" +
    //                      std::to_string(noutput_items) + " offsets: " + offset_str);
    // }

    // for (int i = 0; i < n_max_bursts; ++i) {
    for (const auto& tag : tags) {
        // const auto& tag = tags[i];
        const int burst_start = tag.offset - nitems_read(0);
        const int actual_start = burst_start - d_tag_backoff;

        const pmt::pmt_t& info = tag.value;
        const uint64_t xcorr_idx = pmt::to_uint64(
            pmt::dict_ref(info, pmt::mp("xcorr_idx"), pmt::from_uint64(0)));
        const uint64_t xcorr_offset = pmt::to_uint64(
            pmt::dict_ref(info, pmt::mp("xcorr_offset"), pmt::from_uint64(0)));
        // const std::string src_str =
        //     pmt::is_symbol(tag.srcid) ? pmt::symbol_to_string(tag.srcid) : "N/A";

        if (avail_items - burst_start >= d_burst_len &&
            produced_items + d_burst_len <= noutput_items) {
            // GR_LOG_DEBUG(d_logger,
            //              "Burst " + std::to_string(tags.size()) + "/" +
            //                  std::to_string(n_out_bursts) +
            //                  "\tburst_idx=" + std::to_string(d_frame_counter) + " @" +
            //                  std::to_string(tag.offset) +
            //                  " xcorr_offset=" + std::to_string(xcorr_offset) +
            //                  " xcorr_idx: " + std::to_string(xcorr_idx) +
            //                  " src: " + src_str);

            const float scale_factor = get_scale_factor(info);
            if (actual_start < 0) {
                const int num_prepend_zeros = std::abs(actual_start);
                memset(out, 0, sizeof(gr_complex) * num_prepend_zeros);
                normalize_power_level(out + num_prepend_zeros,
                                      in,
                                      scale_factor,
                                      d_burst_len - num_prepend_zeros);
            } else {
                normalize_power_level(out, in + actual_start, scale_factor, d_burst_len);
            }

            if (d_activate_cfo_correction) {
                compensate_cfo(out, out, get_phase_rotation(info), d_burst_len);
            }
            auto value = pmt::dict_add(
                info, pmt::intern("burst_idx"), pmt::from_uint64(d_frame_counter));
            add_item_tag(0,
                         nitems_written(0) + produced_items,
                         d_burst_start_tag,
                         value,
                         pmt::string_to_symbol(name()));

            d_last_xcorr_offset = xcorr_offset;
            d_last_xcorr_idx = xcorr_idx;

            d_frame_counter++;
            produced_items += d_burst_len;
            consumed_items = burst_start + d_burst_len;
            out += d_burst_len;
        } else {
            // GR_LOG_DEBUG(d_logger,
            //              "Again " + std::to_string(tags.size()) + "/" +
            //                  std::to_string(n_out_bursts) +
            //                  "\tburst_idx=" + std::to_string(d_frame_counter) + " @" +
            //                  std::to_string(tag.offset) +
            //                  " xcorr_offset=" + std::to_string(xcorr_offset) +
            //                  " xcorr_idx: " + std::to_string(xcorr_idx) +
            //                  " src: " + src_str);

            d_expected_xcorr_idx = xcorr_idx;

            consumed_items = std::max(0, burst_start);
            break;
        }
    }

    // if (tags.size() > 0) {
    //     GR_LOG_DEBUG(d_logger,
    //                  "Return: nitems_read=" + std::to_string(nitems_read(0)) +
    //                      "\tavail=" + std::to_string(avail_items) +
    //                      "\tnout_items=" + std::to_string(noutput_items) +
    //                      "\tconsumed=" + std::to_string(consumed_items) +
    //                      "\tproduced=" + std::to_string(produced_items));
    // }

    consume_each(consumed_items);
    return produced_items;
}

} /* namespace gfdm */
} /* namespace gr */
