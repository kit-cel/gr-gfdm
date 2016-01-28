/* -*- c++ -*- */
/*
 * Copyright 2016 Andrej Rode
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_GFDM_RECEIVER_H
#define INCLUDED_GFDM_RECEIVER_H
#endif

#include <gfdm/api.h>
#include <gfdm/gfdm_utils.h>
#include <gnuradio/fft/fft.h>
#include <gnuradio/gr_complex.h>
#include <volk/volk.h>

namespace gr {
  namespace gfdm {
    namespace kernel {
      
      class GFDM_API gfdm_receiver
      {
        protected:
          int d_nsubcarrier;
          int d_ntimeslots;
          int d_filter_width;
          int d_N;
          int d_fft_len;
          std::vector<gr_complex> d_filter_taps;
          std::vector< std::vector<gr_complex> > *sc_fdomain;
          fft::fft_complex *d_in_fft;
          gr_complex *d_in_fft_in;
          gr_complex *d_in_fft_out;
          fft::fft_complex *d_sc_ifft;
          gr_complex *d_sc_ifft_in;
          gr_complex *d_sc_ifft_out;

          void filter_superposition(std::vector< std::vector<gr_complex> > &out, const gr_complex in[]);
          void demodulate_subcarrier(gr_complex out[], std::vector< std::vector<gr_complex> > &sc_fdomain);

        public:
          gfdm_receiver(int nsubcarrier, int ntimeslots, double filter_alpha, int fft_len);
          ~gfdm_receiver();
          


      };
    } /* namespace kernel */
  } /* namespace filter */
} /* namespace gr */


