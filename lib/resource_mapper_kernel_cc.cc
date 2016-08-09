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


#include <gfdm/resource_mapper_kernel_cc.h>
#include <string.h>
#include <stdexcept>

namespace gr {
  namespace gfdm {

    resource_mapper_kernel_cc::resource_mapper_kernel_cc(int active_subcarriers, int fft_len, int timeslots, std::vector<int> subcarrier_map, bool per_timeslot)
        : d_active_subcarriers(active_subcarriers), d_fft_len(fft_len), d_timeslots(timeslots), d_per_timeslot(per_timeslot)
    {
      if (active_subcarriers > fft_len){
        std::stringstream sstm;
        sstm << "active_subcarriers(" << active_subcarriers << ") MUST be smaller or equal to fft_len(" << fft_len << ")!";
        std::string err_str = sstm.str();
        throw std::invalid_argument(err_str.c_str());
      }
      if (int(subcarrier_map.size()) != active_subcarriers){
        std::stringstream sstm;
        sstm << "number of subcarrier_map entries(" << subcarrier_map.size() << ") MUST be equal to active_subcarriers(";
        sstm << active_subcarriers << ")!";
        std::string err_str = sstm.str();
        throw std::invalid_argument(err_str.c_str());
      }
      std::sort(subcarrier_map.begin(), subcarrier_map.end());
      if ( !(adjacent_find(subcarrier_map.begin(), subcarrier_map.end()) == subcarrier_map.end()) ){
        throw std::invalid_argument("All entries in subcarrier_map MUST be unique!");
      }
      if (*std::min_element(subcarrier_map.begin(), subcarrier_map.end()) < 0){
        throw std::invalid_argument("All subcarrier indices MUST be greater or equal to ZERO!");
      }
      if (*std::max_element(subcarrier_map.begin(), subcarrier_map.end()) > fft_len){
        throw std::invalid_argument("All subcarrier indices MUST be smaller or equal to fft_len!");
      }
      d_subcarrier_map = subcarrier_map;
    }

    resource_mapper_kernel_cc::~resource_mapper_kernel_cc()
    {
    }

    void
    resource_mapper_kernel_cc::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size)
    {
      if (ninput_size > input_vector_size()){
        std::stringstream sstm;
        sstm << "input vector size(" << ninput_size << ") MUST not exceed active_subcarriers * timeslots(" << d_active_subcarriers * d_timeslots << ")!";
        throw std::invalid_argument(sstm.str().c_str());
      }
      memset(p_out, 0x0, sizeof(gfdm_complex) * output_vector_size());
      if(d_per_timeslot){
        map_per_timeslot(p_out, p_in, ninput_size);
      }
      else{
        map_per_subcarrier(p_out, p_in, ninput_size);
      }
    }

    void
    resource_mapper_kernel_cc::map_per_timeslot(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size)
    {
      int sym_ctr = 0;
      for(int tidx = 0; tidx < d_timeslots; ++tidx){
        for(int i = 0; i < d_active_subcarriers; ++i){
          int sidx = d_subcarrier_map.at(i);
          int out_idx = d_timeslots * sidx + tidx;
          p_out[out_idx] = sym_ctr < ninput_size ? *p_in++ : gfdm_complex(0.0, 0.0);
          sym_ctr++;
        }
      }
    }

    void
    resource_mapper_kernel_cc::map_per_subcarrier(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size)
    {
      int sym_ctr = 0;
      for(int i = 0; i < d_active_subcarriers; ++i){
        int sidx = d_subcarrier_map.at(i);
        for(int tidx = 0; tidx < d_timeslots; ++tidx){
          int out_idx = d_timeslots * sidx + tidx;
          p_out[out_idx] = sym_ctr < ninput_size ? *p_in++ : gfdm_complex(0.0, 0.0);
          sym_ctr++;
        }
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

