/* -*- c++ -*- */
/* 
 * Copyright 2015 Andrej Rode.
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

#ifndef INCLUDED_GFDM_TRANSMITTER_CVC_IMPL_H
#define INCLUDED_GFDM_TRANSMITTER_CVC_IMPL_H

#include <gfdm/transmitter_cvc.h>
#include <gnuradio/fft/fft.h>
#include <gnuradio/filter/firdes.h>

namespace gr {
  namespace gfdm {

    class transmitter_cvc_impl : public transmitter_cvc
    {
     private:
       int d_nsubcarrier;
       int d_ntimeslots;
       std::vector<gr_complex> d_filtertaps;
       int d_symbols_per_set;
       int d_filter_width;
       fft::fft_complex *d_sc_fft;
       gr_complex * d_sc_fft_in;
       gr_complex * d_sc_fft_out;
       fft::fft_complex *d_out_ifft;
       gr_complex * d_out_ifft_in;
       gr_complex * d_out_ifft_out;
       int mod(int k, int n);


     public:
      transmitter_cvc_impl(
                    int nsubcarrier,
                    int ntimeslots,
                    int filter_width,
                    double filter_alpha);
      ~transmitter_cvc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_TRANSMITTER_CVC_IMPL_H */

