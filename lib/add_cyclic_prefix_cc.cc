/* -*- c++ -*- */
/*
 * Copyright 2016 Johannes Demel.
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

#include <gfdm/add_cyclic_prefix_cc.h>
#include <string.h>
#include <volk/volk.h>
#include <iostream>

namespace gr {
namespace gfdm {

add_cyclic_prefix_cc::add_cyclic_prefix_cc(int block_len,
                                           int cp_len,
                                           int cs_len,
                                           int ramp_len,
                                           std::vector<gfdm_complex> window_taps,
                                           int cyclic_shift)
    : d_block_len(block_len),
      d_cp_len(cp_len),
      d_cs_len(cs_len),
      d_ramp_len(ramp_len),
      d_cyclic_shift(cyclic_shift)
{
    int window_len = block_len + cp_len + cs_len;
    if (window_taps.size() != (unsigned int)window_len &&
        window_taps.size() != (unsigned int)2 * ramp_len) {
        std::stringstream sstm;
        sstm << "number of window taps(" << window_taps.size()
             << ") MUST be equal to 2*ramp_len(";
        sstm << 2 * ramp_len << ") OR block_len+cp_len (" << window_len << ")!";
        throw std::invalid_argument(sstm.str().c_str());
    }

    d_front_ramp =
        volk::vector<gfdm_complex>(window_taps.begin(), window_taps.begin() + ramp_len);
    const unsigned offset = window_taps.size() - ramp_len;
    d_back_ramp =
        volk::vector<gfdm_complex>(window_taps.begin() + offset, window_taps.end());
}

add_cyclic_prefix_cc::~add_cyclic_prefix_cc() {}

void add_cyclic_prefix_cc::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in)
{
    const unsigned cp_start = block_size() - d_cp_len + d_cyclic_shift;
    const unsigned shifted_cp_len = d_cp_len - d_cyclic_shift;
    memcpy(p_out, p_in + cp_start, sizeof(gfdm_complex) * shifted_cp_len);
    const unsigned out_block_start = d_cp_len - d_cyclic_shift;
    memcpy(p_out + shifted_cp_len, p_in, sizeof(gfdm_complex) * block_size());
    const unsigned shifted_cs_len = d_cs_len + d_cyclic_shift;
    memcpy(p_out + shifted_cp_len + block_size(),
           p_in,
           sizeof(gfdm_complex) * shifted_cs_len);

    if (d_ramp_len > 0) {
        const unsigned tail_start = block_size() + d_cp_len + d_cs_len - d_ramp_len;
        volk_32fc_x2_multiply_32fc(p_out, p_out, d_front_ramp.data(), d_ramp_len);
        volk_32fc_x2_multiply_32fc(
            p_out + tail_start, p_out + tail_start, d_back_ramp.data(), d_ramp_len);
    }
}

} /* namespace gfdm */
} /* namespace gr */
