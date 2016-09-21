/* -*- c++ -*- */
/* 
 * Copyright 2016 Andrej Rode, Johannes Demel.
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
#include "advanced_receiver_sb_cc_impl.h"
#include <volk/volk.h>

namespace gr
{
  namespace gfdm
  {

    advanced_receiver_sb_cc::sptr
    advanced_receiver_sb_cc::make(int n_timeslots, int n_subcarriers, int overlap, int ic_iter,
                                  std::vector<gr_complex> frequency_taps, gr::digital::constellation_sptr constellation,
                                  std::vector<int> subcarrier_map)
    {
      return gnuradio::get_initial_sptr
              (new advanced_receiver_sb_cc_impl(n_timeslots, n_subcarriers, overlap, ic_iter, frequency_taps,
                                                constellation, subcarrier_map));
    }

    /*
     * The private constructor
     */
    advanced_receiver_sb_cc_impl::advanced_receiver_sb_cc_impl(int n_timeslots, int n_subcarriers, int overlap,
                                                               int ic_iter, std::vector<gr_complex> frequency_taps,
                                                               gr::digital::constellation_sptr constellation,
                                                               std::vector<int> subcarrier_map)
            : gr::sync_block("advanced_receiver_sb_cc",
                             gr::io_signature::make(1, 1, sizeof(gr_complex)),
                             gr::io_signature::make(1, 1, sizeof(gr_complex))),
              d_n_timeslots(n_timeslots),
              d_n_subcarriers(n_subcarriers),
              d_ic_iter(ic_iter),
              d_constellation(constellation),
              d_subcarrier_map(subcarrier_map)
    {
      d_kernel = receiver_kernel_cc::sptr(new receiver_kernel_cc(n_timeslots, n_subcarriers, overlap, frequency_taps));
      set_output_multiple(d_kernel->block_size());

      //Initialize buffers for temporary subcarrier data
      d_freq_block = (gr_complex *) volk_malloc(sizeof(gr_complex) * d_kernel->block_size(), volk_get_alignment());
      d_ic_time_buffer = (gr_complex *) volk_malloc(sizeof(gr_complex) * d_kernel->block_size(), volk_get_alignment());
      d_ic_freq_buffer = (gr_complex *) volk_malloc(sizeof(gr_complex) * d_kernel->block_size(), volk_get_alignment());
    }

    /*
     * Our virtual destructor.
     */
    advanced_receiver_sb_cc_impl::~advanced_receiver_sb_cc_impl()
    {
      volk_free(d_ic_freq_buffer);
      volk_free(d_freq_block);
      volk_free(d_ic_time_buffer);
    }

    int
    advanced_receiver_sb_cc_impl::work(int noutput_items,
                                       gr_vector_const_void_star &input_items,
                                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      const int n_blocks = noutput_items / d_kernel->block_size();

      for (int i = 0; i < n_blocks; ++i) {
        demodulate_block_ic_array(out, in);

        in += d_kernel->block_size();
        out += d_kernel->block_size();
      }
      return noutput_items;
    }

    void advanced_receiver_sb_cc_impl::demodulate_block_ic_array(gr_complex *p_out, const gr_complex *p_in)
    {
      d_kernel->fft_filter_downsample(d_freq_block, p_in);
      d_kernel->transform_subcarriers_to_td(p_out, d_freq_block);

      for (int j = 0; j < d_ic_iter; ++j) {
        map_symbols_to_constellation_points(d_ic_time_buffer, p_out);
        d_kernel->cancel_sc_interference(d_ic_freq_buffer, d_ic_time_buffer, d_freq_block);
        d_kernel->transform_subcarriers_to_td(p_out, d_ic_freq_buffer);
      }
    }

    void advanced_receiver_sb_cc_impl::map_symbols_to_constellation_points(gr_complex *p_out, const gr_complex *p_in)
    {
      memset(p_out, 0x0, sizeof(gr_complex) * d_kernel->block_size());
      unsigned int symbol_tmp = 0;
      std::vector<gr_complex> const_points = d_constellation->points();
      // perform symbol decision for active subcarriers. All others should be set to '0'!
      for (int j = 0; j < d_subcarrier_map.size(); ++j) {
        int k = d_subcarrier_map.at(j);
        for (int m = 0; m < d_n_timeslots; ++m) {
          symbol_tmp = d_constellation->decision_maker(&p_in[k * d_n_timeslots + m]);
          p_out[k * d_n_timeslots + m] = const_points[symbol_tmp];
        }
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

