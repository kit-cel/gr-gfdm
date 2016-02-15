/* -*- c++ -*- */
/* 
 * Copyright 2016 Andrej Rode.
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

#ifndef INCLUDED_GFDM_MODULATOR_CC_IMPL_H
#define INCLUDED_GFDM_MODULATOR_CC_IMPL_H

#include <gnuradio/fft/fft.h>
#include <gfdm/modulator_cc.h>
#include <gnuradio/filter/firdes.h>
#include <pmt/pmt.h>
#include <volk/volk.h>
#include <gfdm/gfdm_utils.h>

namespace gr {
  namespace gfdm {

    class modulator_cc_impl : public modulator_cc
    {
     private:
       int d_ntimeslots;
       int d_nsubcarrier;
       int d_filter_width;
       int d_N;
       int d_fft_len;
       int d_sync_fft_len;
       std::string d_len_tag_key;
       std::vector<gr_complex> d_filter_taps;
       fft::fft_complex *d_sc_fft;
       gr_complex * d_sc_fft_in;
       gr_complex * d_sc_fft_out;
       fft::fft_complex *d_sync_ifft;
       gr_complex * d_sync_ifft_in;
       gr_complex * d_sync_ifft_out;
       fft::fft_complex *d_out_ifft;
       gr_complex * d_out_ifft_in;
       gr_complex * d_out_ifft_out;
      // Nothing to declare in this block.

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);
      virtual void update_length_tags  (int n_produced, int n_ports);

     public:
      modulator_cc_impl(
          int nsubcarrier,
          int ntimeslots,
          double filter_alpha,
          int fft_len,
          int sync_fft_len,
          const std::string& len_tag_key);
      ~modulator_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_MODULATOR_CC_IMPL_H */

