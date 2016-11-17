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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gfdm/advanced_receiver_kernel_cc.h>
#include <volk/volk.h>

namespace gr {
  namespace gfdm {

    advanced_receiver_kernel_cc::advanced_receiver_kernel_cc(int timeslots, int subcarriers, int overlap, std::vector<gr_complex> frequency_taps,
                                                             std::vector<int> subcarrier_map, int ic_iter,
                                                             gr::digital::constellation_sptr constellation):
            d_ic_iter(ic_iter),
            d_constellation(constellation),
            d_subcarrier_map(subcarrier_map)
    {
      d_kernel = receiver_kernel_cc::sptr(new receiver_kernel_cc(timeslots, subcarriers, overlap, frequency_taps));

      //Initialize buffers for temporary subcarrier data
      d_freq_block = (gr_complex *) volk_malloc(sizeof(gr_complex) * d_kernel->block_size(), volk_get_alignment());
      d_ic_time_buffer = (gr_complex *) volk_malloc(sizeof(gr_complex) * d_kernel->block_size(), volk_get_alignment());
      d_ic_freq_buffer = (gr_complex *) volk_malloc(sizeof(gr_complex) * d_kernel->block_size(), volk_get_alignment());
    }

    advanced_receiver_kernel_cc::~advanced_receiver_kernel_cc()
    {
      volk_free(d_ic_freq_buffer);
      volk_free(d_freq_block);
      volk_free(d_ic_time_buffer);
    }

    void advanced_receiver_kernel_cc::generic_work(gr_complex *p_out, const gr_complex *p_in)
    {
      d_kernel->fft_filter_downsample(d_freq_block, p_in);
      d_kernel->transform_subcarriers_to_td(p_out, d_freq_block);

      for (int j = 0; j < d_ic_iter; ++j) {
        map_symbols_to_constellation_points(d_ic_time_buffer, p_out);
        d_kernel->cancel_sc_interference(d_ic_freq_buffer, d_ic_time_buffer, d_freq_block);
        d_kernel->transform_subcarriers_to_td(p_out, d_ic_freq_buffer);
      }
    }

    void advanced_receiver_kernel_cc::map_symbols_to_constellation_points(gr_complex *p_out, const gr_complex *p_in)
    {
      memset(p_out, 0x0, sizeof(gr_complex) * d_kernel->block_size());
      unsigned int symbol_tmp = 0;
      std::vector<gr_complex> const_points = d_constellation->points();
      // perform symbol decision for active subcarriers. All others should be set to '0'!
      for (int j = 0; j < d_subcarrier_map.size(); ++j) {
        int k = d_subcarrier_map.at(j);
        for (int m = 0; m < d_kernel->timeslots(); ++m) {
          symbol_tmp = d_constellation->decision_maker(&p_in[k * d_kernel->timeslots() + m]);
          p_out[k * d_kernel->timeslots() + m] = const_points[symbol_tmp];
        }
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

