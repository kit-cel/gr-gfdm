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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gfdm/gfdm_utils.h>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace gfdm {

    rrc_filter_sparse::rrc_filter_sparse(
        int ntaps,
        double alpha,
        int filter_width,
        int nsubcarrier,
        int ntimeslots)
    {
      std::vector<float> filtertaps(ntaps);
      std::vector<float> filtertaps_center = gr::filter::firdes::root_raised_cosine(
          1.0,
          1.0,
          double(1.0/nsubcarrier),
          alpha,
          ntaps);
      for (int i=0;i<ntaps;i++)
      {
        filtertaps[i] = filtertaps_center[(i+ntaps/2) % ntaps];
      }
      //Initialize FFT
      fft::fft_real_fwd *filter_fft = new fft::fft_real_fwd(ntaps,1);
      float *in = filter_fft->get_inbuf();
      gr_complex *out = filter_fft->get_outbuf();
      std::memset((void*) in, 0x00, sizeof(float)*ntaps);
      //Copy Filtertaps in FFT Input
      std::memcpy(&in[0], &filtertaps[0], sizeof(float)*ntaps);
      filter_fft->execute();
      d_filter_taps.resize(ntimeslots*filter_width,gr_complex(0.0f, 0.0f));
      // Only works for d_filter_width = 2 needs some rework for d_filter_width other than 2
      std::memcpy(&d_filter_taps[0], out, sizeof(gr_complex)*ntimeslots);
      for (int i=0; i<ntimeslots-1; i++)
      {
        d_filter_taps[i+ntimeslots+1] = std::conj(d_filter_taps[ntimeslots-1-i]);
      }
      delete filter_fft;

    };
    rrc_filter_sparse::~rrc_filter_sparse()
    {

    };
    void
    rrc_filter_sparse::get_taps(std::vector<gr_complex> &out)
    {
      out.resize(d_filter_taps.size());
      std::memcpy(&out[0],&d_filter_taps[0],sizeof(gr_complex)*d_filter_taps.size());
    };

  } /* namespace gfdm */
} /* namespace gr */
