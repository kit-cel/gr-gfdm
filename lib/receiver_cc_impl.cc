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
#include "receiver_cc_impl.h"

namespace gr {
  namespace gfdm {

    receiver_cc::sptr
    receiver_cc::make(int nsubcarrier,
                      int ntimeslots,
                      double filter_alpha,
                      int fft_len,
                      const std::string& len_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new receiver_cc_impl(nsubcarrier,
                              ntimeslots,
                              filter_alpha,
                              fft_len,
                              len_tag_key));
    }

    /*
     * The private constructor
     */
    receiver_cc_impl::receiver_cc_impl(int nsubcarrier,
                                        int ntimeslots,
                                        double filter_alpha,
                                        int fft_len,
                                        const std::string& len_tag_key)
      : gr::tagged_stream_block("receiver_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), len_tag_key),
      d_nsubcarrier(nsubcarrier),
      d_ntimeslots(ntimeslots),
      d_N(ntimeslots*nsubcarrier),
      d_fft_len(fft_len)
    {
      d_filter_width = 2;
      std::vector<float> filtertaps = gr::filter::firdes::root_raised_cosine(
          1.0,
          1.0,
          double(1.0/d_nsubcarrier),
          filter_alpha,
          d_N);
      fft::fft_real_fwd *filter_fft = new fft::fft_real_fwd(d_N,1);
      float *in = filter_fft->get_inbuf();
      gr_complex *out = filter_fft->get_outbuf();
      std::memset((void*) in, 0x00, sizeof(float)*d_N);
      //Copy Filtertaps in FFT Input
      std::memcpy((void*) in, &filtertaps[0], sizeof(float)*d_N);
      filter_fft->execute();
      d_filter_taps.resize(d_ntimeslots*d_filter_width,0j);
      // Only works for d_filter_width = 2 needs some rework for d_filter_width other than 
      std::memcpy(&d_filter_taps[0], out, sizeof(gr_complex)*d_ntimeslots);
      for (int i=0; i<d_ntimeslots-1; i++)
      {
        d_filter_taps[i+d_ntimeslots+1] = std::conj(d_filter_taps[d_ntimeslots-1-i]);
      }
      delete filter_fft;

      //Initialize input FFT
      d_in_fft = new fft::fft_complex(d_fft_len,true,1);
      d_in_fft_in = d_in_fft->get_inbuf();
      d_in_fft_out = d_in_fft->get_outbuf();
      
      //Initialize IFFT per subcarrier
      d_sc_ifft = new fft::fft_complex(d_ntimeslots,false,1);
      d_sc_ifft_in = d_sc_ifft->get_inbuf();
      d_sc_ifft_out = d_sc_ifft->get_outbuf();

    }


    /*
     * Our virtual destructor.
     */
    receiver_cc_impl::~receiver_cc_impl()
    {
    }

    int
    receiver_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = d_nsubcarrier*d_ntimeslots;
      if (ninput_items[0] != d_fft_len)
      {
        throw std::runtime_error("frame_len must be equal to fft_len");
      }
      return noutput_items;
    }

    int
    receiver_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];
        std::memset(&out[0],0x00,sizeof(gr_complex)*d_N);

        // 1. FFT on input
        std::memcpy(&d_in_fft_in[0],&in[0],sizeof(gr_complex)*d_fft_len);
        d_in_fft->execute();
        // 2. Extract subcarrier
        // 3. Filter and superposition every subcarrier
        for (int k=0; k<d_nsubcarrier; k++)
        {
          std::vector<gr_complex> sc_postfft(d_ntimeslots*d_filter_width);
          std::vector<gr_complex> sc_postfilter(d_ntimeslots*d_filter_width);
          std::vector<gr_complex> sc_postifft(d_ntimeslots);
          //Subcarrier-Offset = d_fft_len/2 + (d_fft_len-d_N)/2 - ((d_filter_width-1)*(d_ntimeslots))/2 + k*d_ntimeslots ) % d_fft_len
          int sc_offset = (d_fft_len/2 + (d_fft_len - d_N)/2 - ((d_filter_width-1)*(d_ntimeslots))/2 + k*d_ntimeslots) % d_fft_len;
          gr_complex * sc = &d_in_fft_out[sc_offset];
          // Only valid for d_filter_width = 2
          for (int n=0;n<d_filter_width*d_ntimeslots;n++)
          {
            sc_postfft.push_back(sc[n + d_ntimeslots % d_ntimeslots*d_filter_width]);

          }

          ::volk_32fc_x2_multiply_32fc(&sc_postfilter[0],
              &sc_postfft[0],&d_filter_taps[0],d_ntimeslots*d_filter_width);
          // Only valid for d_filter_width = 2
          ::volk_32f_x2_add_32f((float*)&d_sc_ifft_in[0],
              (float*)&sc_postfilter[0],(float*)&sc_postfilter[d_ntimeslots],2*d_ntimeslots);
        // 4. apply ifft on every filtered and superpositioned subcarrier
          d_sc_ifft->execute();
          ::volk_32fc_s32fc_multiply_32fc(&sc_postifft[0],&d_sc_ifft_out[0],static_cast<gr_complex>(1.0/d_ntimeslots),d_ntimeslots);

          for (int m=0; m<d_ntimeslots; m++)
          {
            out[k+m*d_nsubcarrier] = sc_postifft[m];
          }
                              
        }
        // (5.) Provide hooks for advanced IC-Receiver

        return d_N;
    }

  } /* namespace gfdm */
} /* namespace gr */

