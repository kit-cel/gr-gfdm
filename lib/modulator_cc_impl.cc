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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "modulator_cc_impl.h"

namespace gr {
  namespace gfdm {

    modulator_cc::sptr
    modulator_cc::make(
        const std::string& len_tag_key,
        int nsubcarrier,
        int ntimeslots,
        double filter_alpha,
        int fft_len)
    {
      return gnuradio::get_initial_sptr
        (new modulator_cc_impl(len_tag_key,
                               nsubcarrier,
                               ntimeslots,
                               filter_alpha,
                               fft_len)
         );

    }

    /*
     * The private constructor
     */
    modulator_cc_impl::modulator_cc_impl(
        const std::string& len_tag_key,
        int nsubcarrier,
        int ntimeslots,
        double filter_alpha,
        int fft_len)
      : gr::tagged_stream_block("modulator_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              len_tag_key),
      d_ntimeslots(ntimeslots),
      d_nsubcarrier(nsubcarrier),
      d_N(ntimeslots*nsubcarrier),
      d_fft_len(fft_len)
    {
      if (d_fft_len < d_N)
      {
        throw std::invalid_argument("fft_len must be greater than or equal to nsubcarrier*ntimeslots");
      }
      int d_filter_width = 2;
      std::vector<float> filtertaps = gr::filter::firdes::root_raised_cosine(
          double(d_nsubcarrier),
          1.0,
          double(1.0/d_nsubcarrier),
          filter_alpha,
          d_N);
      fft::fft_real_fwd *filter_fft = new fft::fft_real_fwd(d_N,1);
      float *in = filter_fft->get_inbuf();
      gr_complex *out = filter_fft->get_outbuf();

      //Copy Filtertaps in FFT Input
      std::memcpy((void*) in, &filtertaps[0], sizeof(gr_complex)*d_N);
      
      filter_fft->execute();
      d_filter_taps.reserve(d_filter_width*d_ntimeslots);
      // Only works for d_filter_width = 2 needs some rework for d_filter_width other than 
      std::memcpy(&d_filter_taps[0], out, sizeof(gr_complex)*d_ntimeslots);
      d_filter_taps[d_ntimeslots] = gr_complex(0j);
      for (int i=0; i<d_ntimeslots-1; i++)
      {
        d_filter_taps[i+d_ntimeslots+1] = std::conj(d_filter_taps[d_ntimeslots-1-i]);
      }
      delete filter_fft;

      //Initialize FFT per subcarrier
      d_sc_fft = new fft::fft_complex(d_ntimeslots,true,1);
      d_sc_fft_in = d_sc_fft->get_inbuf();
      d_sc_fft_out = d_sc_fft->get_outbuf();


      //Initialize resulting FFT
      d_out_ifft = new fft::fft_complex(d_fft_len,false,1);
      d_out_ifft_in = d_out_ifft->get_inbuf();
      d_out_ifft_out = d_out_ifft->get_outbuf();

    }

    /*
     * Our virtual destructor.
     */
    modulator_cc_impl::~modulator_cc_impl()
    {
      delete d_sc_fft;
      delete d_out_ifft;
    }

    int
    modulator_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0];
      return noutput_items ;
    }

    int
    modulator_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];

        // Do signal processing!


        return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

