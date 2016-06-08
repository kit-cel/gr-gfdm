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
#include <volk/volk.h>
#include <iostream>
#include <string.h>

namespace gr {
  namespace gfdm {

    add_cyclic_prefix_cc::add_cyclic_prefix_cc(int ramp_len, int cp_len, int block_len, std::vector<gfdm_complex> window_taps)
            : d_ramp_len(ramp_len), d_cp_len(cp_len), d_block_len(block_len)
    {
      set_block_size(block_len);
      int window_len = block_len + cp_len;
      if(window_taps.size() != window_len && window_taps.size() != 2 * ramp_len){
        std::cout << "window_len: " << window_len << ", ramp_len: " << ramp_len << ", taps.size: " << window_taps.size() << std::endl;
        throw std::invalid_argument("ERROR: number of window_taps elements MUST be equal to 2*ramp_len OR n_timeslots*n_subcarriers+cp_len!");
      }
    }

    add_cyclic_prefix_cc::~add_cyclic_prefix_cc()
    {
    }

    void
    add_cyclic_prefix_cc::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      int cp_start = block_size() - d_cp_len;
      std::cout << "block_len: " << block_size() << ", cp_start: " << cp_start << ", cp_len: " << d_cp_len << std::endl;
      memcpy(p_out, p_in + block_size() - d_cp_len, sizeof(gfdm_complex) * d_cp_len);
      memcpy(p_out + d_cp_len, p_in, sizeof(gfdm_complex) * block_size());
    }

  } /* namespace gfdm */
} /* namespace gr */

