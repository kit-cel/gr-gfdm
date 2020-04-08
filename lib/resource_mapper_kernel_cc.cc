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
#include <cstring>
#include <stdexcept>
#include <algorithm>

namespace gr {
  namespace gfdm {

    resource_mapper_kernel_cc::resource_mapper_kernel_cc(int timeslots,
                                                         int subcarriers,
                                                         int active_subcarriers,
                                                         std::vector<int> subcarrier_map,
                                                         bool per_timeslot,
                                                         bool is_mapper)
        : d_timeslots(timeslots),
          d_subcarriers(subcarriers),
          d_active_subcarriers(active_subcarriers),
          d_block_size(timeslots * active_subcarriers),
          d_frame_size(timeslots * subcarriers),
          d_per_timeslot(per_timeslot),
          d_is_mapper(is_mapper)
    {
      if (active_subcarriers > subcarriers){
        throw std::invalid_argument("active_subcarriers(" +
                                    std::to_string(active_subcarriers) +
                                    ") MUST be smaller or equal to subcarriers(" +
                                    std::to_string(subcarriers) +
                                    ")!");
      }
      if (int(subcarrier_map.size()) != active_subcarriers){
        throw std::invalid_argument("number of subcarrier_map entries(" +
                                    std::to_string(subcarrier_map.size()) +
                                    ") MUST be equal to active_subcarriers(" +
                                    std::to_string(active_subcarriers) +
                                    ")!");
      }
      std::sort(subcarrier_map.begin(), subcarrier_map.end());
      if ( !(adjacent_find(subcarrier_map.begin(), subcarrier_map.end()) == subcarrier_map.end()) ){
        throw std::invalid_argument("All entries in subcarrier_map MUST be unique!");
      }
      if (*std::min_element(subcarrier_map.begin(), subcarrier_map.end()) < 0){
        throw std::invalid_argument("All subcarrier indices MUST be greater or equal to ZERO!");
      }
      if (*std::max_element(subcarrier_map.begin(), subcarrier_map.end()) > subcarriers){
        throw std::invalid_argument("All subcarrier indices MUST be smaller or equal to subcarriers!");
      }
      d_subcarrier_map = subcarrier_map;
    }

    resource_mapper_kernel_cc::~resource_mapper_kernel_cc()
    {
    }

    void
    resource_mapper_kernel_cc::map_to_resources(gfdm_complex* p_out,
                            const gfdm_complex* p_in,
                            const size_t ninput_size)
    {
      if (ninput_size > block_size()){
        throw std::invalid_argument("input vector size(" +
                                    std::to_string (ninput_size) +
                                    ") MUST not exceed active_subcarriers * timeslots(" +
                                    std::to_string(d_block_size) +
                                    ")!");
      }
      memset(p_out, 0x0, sizeof(gfdm_complex) * frame_size());
      if(d_per_timeslot){
        map_per_timeslot(p_out, p_in, ninput_size);
      }
      else{
        map_per_subcarrier(p_out, p_in, ninput_size);
      }
    }

    void
    resource_mapper_kernel_cc::demap_from_resources(gfdm_complex* p_out,
                                const gfdm_complex* p_in,
                                const size_t noutput_size)
    {
      if (noutput_size > block_size()){
        throw std::invalid_argument("output vector size(" +
                                    std::to_string(noutput_size) +
                                    ") MUST not exceed active_subcarriers * timeslots(" +
                                    std::to_string(d_block_size) +
                                    ")!");
      }
      memset(p_out, 0x0, sizeof(gfdm_complex) * noutput_size);
      if(d_per_timeslot){
        demap_per_timeslot(p_out, p_in, noutput_size);
      }
      else{
        demap_per_subcarrier(p_out, p_in, noutput_size);
      }
    }

    void
    resource_mapper_kernel_cc::map_per_timeslot(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size)
    {
      int sym_ctr = 0;
      for(int tidx = 0; tidx < d_timeslots; ++tidx){
        for(const auto sidx : d_subcarrier_map){
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
      for(const auto sidx : d_subcarrier_map){
        for(int tidx = 0; tidx < d_timeslots; ++tidx){
          int out_idx = d_timeslots * sidx + tidx;
          p_out[out_idx] = sym_ctr < ninput_size ? *p_in++ : gfdm_complex(0.0, 0.0);
          sym_ctr++;
        }
      }
    }

        void
    resource_mapper_kernel_cc::demap_per_timeslot(gfdm_complex* p_out, const gfdm_complex* p_in, const int noutput_size)
    {
      for(int i = 0; i < noutput_size; ++i){
        int tidx = i / d_active_subcarriers;
        int sidx = d_subcarrier_map.at(i % d_active_subcarriers);
        *p_out++ = p_in[d_timeslots * sidx + tidx];
      }
    }

    void
    resource_mapper_kernel_cc::demap_per_subcarrier(gfdm_complex* p_out, const gfdm_complex* p_in, const int noutput_size)
    {
      int sym_ctr = 0;
      for(const auto sidx : d_subcarrier_map){
        for(int tidx = 0; tidx < d_timeslots; ++tidx){
          int idx = d_timeslots * sidx + tidx;
          *p_out++ = p_in[idx];
          sym_ctr++;
          if(sym_ctr > noutput_size){
            return;
          }
        }
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

