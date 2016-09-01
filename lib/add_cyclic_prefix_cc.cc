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
      int window_len = block_len + cp_len;
      if(window_taps.size() != (unsigned int) window_len && window_taps.size() != (unsigned int) 2 * ramp_len){
        std::stringstream sstm;
        sstm << "number of window taps(" << window_taps.size() << ") MUST be equal to 2*ramp_len(";
        sstm << 2 * ramp_len << ") OR block_len+cp_len (" << window_len << ")!";
        throw std::invalid_argument(sstm.str().c_str());
      }

      d_front_ramp = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * ramp_len, volk_get_alignment());
      d_back_ramp = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * ramp_len, volk_get_alignment());
      memcpy(d_front_ramp, &window_taps[0], sizeof(gfdm_complex) * ramp_len);
      memcpy(d_back_ramp, &window_taps[block_len + cp_len - ramp_len], sizeof(gfdm_complex) * ramp_len);
    }

    add_cyclic_prefix_cc::~add_cyclic_prefix_cc()
    {
      volk_free(d_front_ramp);
      volk_free(d_back_ramp);
    }

    void
    add_cyclic_prefix_cc::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      const int cp_start = block_size() - d_cp_len;
      memcpy(p_out, p_in + cp_start, sizeof(gfdm_complex) * d_cp_len);
      memcpy(p_out + d_cp_len, p_in, sizeof(gfdm_complex) * block_size());

      if(d_ramp_len > 0){
        const int tail_start = block_size() + d_cp_len - d_ramp_len;
        volk_32fc_x2_multiply_32fc(p_out, p_out, d_front_ramp, d_ramp_len);
        volk_32fc_x2_multiply_32fc(p_out + tail_start, p_out + tail_start, d_back_ramp, d_ramp_len);
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

