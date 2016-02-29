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
#include <gfdm/preamble_generator.h>

namespace gr {
  namespace gfdm {

    preamble_generator::preamble_generator(int nsubcarrier, double filter_alpha, int sync_fft_len)
      :d_sync_fft_len(sync_fft_len)
    {
      float qam_energy = (1/::sqrt(2.0));
      std::vector<float> symbol_choices(2);
      symbol_choices.push_back(-qam_energy);
      symbol_choices.push_back(qam_energy);
      d_symbols.resize(nsubcarrier);
      d_samp_preamble.resize(sync_fft_len);

      for (int sc=0; sc<nsubcarrier;sc++)
      {
        d_symbols[sc] = gr_complex(symbol_choices[std::rand() % 2],symbol_choices[std::rand() %2]);
      }
      //Initialize stuff
      std::vector<gr_complex> filter_taps(2*2);
      //Initialize IFFT
      fft::fft_complex* ifft = new fft::fft_complex(sync_fft_len,false,1);
      gr_complex* ifft_in = ifft->get_inbuf();
      gr_complex* ifft_out = ifft->get_outbuf();
      //Initialize SC_FFT
      fft::fft_complex* sc_fft = new fft::fft_complex(2,true,1);
      gr_complex* sc_fft_in = sc_fft->get_inbuf();
      gr_complex* sc_fft_out = sc_fft->get_outbuf();
      //Initialize SC_Filter
      rrc_filter_sparse* sc_filter = new gfdm::rrc_filter_sparse(2*nsubcarrier, filter_alpha, 2, nsubcarrier, 2);
      // Get sc filtertaps
      sc_filter->get_taps(filter_taps);
      for (int sc=0; sc<nsubcarrier;sc++)
      {
        std::vector<gr_complex> sc_tmp(2*2,0j);
        sc_fft_in[0] = d_symbols[sc];
        sc_fft_in[1] = d_symbols[sc];
        sc_fft->execute();
                
        for (int l=0; l<2;l++)
        {
          ::volk_32fc_x2_multiply_32fc(&sc_tmp[l*2],&filter_taps[l*2],&sc_fft_out[0],2);
        }

        int ifft_offset = ( sync_fft_len/2 + (sync_fft_len-2*nsubcarrier)/2 -(sc*2)) % sync_fft_len;
        for (int n=0; n< 2*2; n++)
        {
          ifft_in[(ifft_offset+n) % sync_fft_len] += sc_tmp[(n+(2*2)/2) % (2*2)];
        }      
      
      }
      ifft->execute();
      ::volk_32fc_s32fc_multiply_32fc(&d_samp_preamble[0],&ifft_out[0], static_cast<gr_complex>(1.0/sync_fft_len),sync_fft_len);
      delete ifft;
      delete sc_fft;
      delete sc_filter;

    }

    preamble_generator::~preamble_generator()
    {
    }
    
    preamble_generator::sptr
    preamble_generator::make(int nsubcarrier, double filter_alpha, int sync_fft_len)
    {
      return preamble_generator::sptr(
          new preamble_generator(nsubcarrier, filter_alpha, sync_fft_len));
      
    }

  } /* namespace gfdm */
} /* namespace gr */

