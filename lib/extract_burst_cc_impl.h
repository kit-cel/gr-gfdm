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

#ifndef INCLUDED_GFDM_EXTRACT_BURST_CC_IMPL_H
#define INCLUDED_GFDM_EXTRACT_BURST_CC_IMPL_H

#include <gfdm/extract_burst_cc.h>

namespace gr {
namespace gfdm {

class extract_burst_cc_impl : public extract_burst_cc
{
private:
    const pmt::pmt_t d_scale_factor_key = pmt::mp("scale_factor");
    const pmt::pmt_t d_phase_rotation_key = pmt::mp("sc_rot");

    const int d_burst_len;
    const int d_tag_backoff;
    const pmt::pmt_t d_burst_start_tag;
    bool d_activate_cfo_correction;

    uint64_t d_last_tag_offset = 0;
    uint64_t d_last_xcorr_offset = 0;
    uint64_t d_last_xcorr_idx = 0;
    uint64_t d_expected_xcorr_idx = 0;
    uint64_t d_frame_counter = 0;

    float get_scale_factor(const pmt::pmt_t& info) const;
    gr_complex get_phase_rotation(const pmt::pmt_t& info) const;

    void normalize_power_level(gr_complex* p_out,
                               const gr_complex* p_in,
                               const float norm_factor,
                               const int ninput_size);
    void compensate_cfo(gr_complex* p_out,
                        const gr_complex* p_in,
                        const gr_complex phase_increment,
                        const int ninput_size);

public:
    extract_burst_cc_impl(int burst_len,
                          int tag_backoff,
                          std::string burst_start_tag,
                          bool activate_cfo_correction);
    ~extract_burst_cc_impl();

    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    void activate_cfo_compensation(bool activate_cfo_compensation);
    // Where all the action really happens
    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_EXTRACT_BURST_CC_IMPL_H */
