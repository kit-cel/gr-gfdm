/* -*- c++ -*- */
/*
 * Copyright 2016 Andrej Rode.
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

#include <gfdm/gfdm_receiver.h>

namespace gr {
  namespace gfdm {
    namespace kernel {
          
      gfdm_receiver::gfdm_receiver(int nsubcarrier,
                                   int ntimeslots,
                                   double filter_alpha,
                                   int fft_len)
        : 
        d_nsubcarrier(nsubcarrier),
        d_ntimeslots(ntimeslots),
        d_N(ntimeslots*nsubcarrier),
        d_fft_len(fft_len)

      {
        d_filter_width = 2;
        d_filter_taps.resize(d_ntimeslots*d_filter_width);
        rrc_filter_sparse *filter_gen = new rrc_filter_sparse(d_N,filter_alpha,d_filter_width,nsubcarrier,ntimeslots);
        filter_gen->get_taps(d_filter_taps);
        delete filter_gen;
      
        //Initialize input FFT
        d_in_fft = new fft::fft_complex(d_fft_len,true,1);
        d_in_fft_in = d_in_fft->get_inbuf();
        d_in_fft_out = d_in_fft->get_outbuf();
      
        //Initialize IFFT per subcarrier
        d_sc_ifft = new fft::fft_complex(d_ntimeslots,false,1);
        d_sc_ifft_in = d_sc_ifft->get_inbuf();
        d_sc_ifft_out = d_sc_ifft->get_outbuf();
        //Initialize vector of vectors for temporary subcarrier data
        d_sc_fdomain.resize(nsubcarrier);
        for (std::vector< std::vector<gr_complex> >::iterator it = d_sc_fdomain.begin();it != d_sc_fdomain.end();++it)
        {
          it->resize(ntimeslots);
        }
        d_sc_symbols.resize(nsubcarrier);
        for (std::vector< std::vector<gr_complex> >::iterator it = d_sc_symbols.begin();it != d_sc_symbols.end();++it)
        {
          it->resize(ntimeslots);
        }
      }
     
      gfdm_receiver::~gfdm_receiver()
      {
        delete d_in_fft;
        delete d_sc_ifft;
      }
      
      void
      gfdm_receiver::filter_superposition(std::vector< std::vector<gr_complex> > &out,
          const gr_complex in[] )
      {
        ::volk_32fc_s32fc_multiply_32fc(&d_in_fft_in[0],&in[0],static_cast<gr_complex>(float(d_N)/float(d_fft_len)),d_fft_len);
        //std::memcpy(&d_in_fft_in[0],&in[0],sizeof(gr_complex)*d_fft_len);
        d_in_fft->execute();
        std::vector<gr_complex> fft_out(d_fft_len+d_ntimeslots*d_filter_width);
        std::memcpy(&fft_out[0],&d_in_fft_out[0],sizeof(gr_complex)*d_fft_len);
        std::memcpy(&fft_out[d_fft_len],&d_in_fft_out[0],sizeof(gr_complex)*d_ntimeslots*d_filter_width);
        for (int k=0; k<d_nsubcarrier; k++)
        {
          std::vector<gr_complex> sc_postfft(d_ntimeslots*d_filter_width);
          std::vector<gr_complex> sc_postfilter(d_ntimeslots*d_filter_width);
          //FFT output is not centered:
          //Subcarrier-Offset = d_fft_len/2 + (d_fft_len-d_N)/2 - ((d_filter_width-1)*(d_ntimeslots))/2 + k*d_ntimeslots ) modulo d_fft_len
          int sc_offset = (d_fft_len/2 + (d_fft_len - d_N)/2 - ((d_filter_width-1)*(d_ntimeslots))/2 + k*d_ntimeslots) % d_fft_len;
          gr_complex * sc = &fft_out[sc_offset];
          for (int n=0;n<d_filter_width*d_ntimeslots;n++)
          {
            sc_postfft[n] = sc[(n + (d_ntimeslots*d_filter_width)/2) % (d_ntimeslots*d_filter_width)];
          }
          ::volk_32fc_x2_multiply_32fc(&sc_postfilter[0],
              &sc_postfft[0],&d_filter_taps[0],d_ntimeslots*d_filter_width);
          // Only valid for d_filter_width = 2, write own kernel for complex addition
          ::volk_32f_x2_add_32f((float*)&out[k][0],
              (float*)(&sc_postfilter[0]),(float*)(&sc_postfilter[d_ntimeslots]),2*d_ntimeslots);
          }

      }

      void
      gfdm_receiver::demodulate_subcarrier(std::vector< std::vector<gr_complex> > &out,
          std::vector< std::vector<gr_complex> > &sc_fdomain)
      {
        // 4. apply ifft on every filtered and superpositioned subcarrier
        for (int k=0; k<d_nsubcarrier; k++)
        {
          std::vector<gr_complex> sc_postifft(d_ntimeslots);
          std::memcpy(&d_sc_ifft_in[0],&sc_fdomain[k][0],sizeof(gr_complex)*d_ntimeslots);
          d_sc_ifft->execute();
          ::volk_32fc_s32fc_multiply_32fc(&out[k][0],&d_sc_ifft_out[0],static_cast<gr_complex>(1.0/(float)d_ntimeslots),d_ntimeslots);
        }

      }

      void
      gfdm_receiver::serialize_output(gr_complex out[],
          std::vector< std::vector<gr_complex> > &sc_symbols)
      {
        for (int k=0; k<d_nsubcarrier; k++)
        {
          for (int m=0; m<d_ntimeslots; m++)
          {
            out[k+m*d_nsubcarrier] = sc_symbols[k][m];
          }
        }
      }

      void
      gfdm_receiver::gfdm_work(gr_complex out[],const gr_complex in[], int ninput_items, int noutputitems)
      {
       filter_superposition(d_sc_fdomain,in);
       demodulate_subcarrier(d_sc_symbols,d_sc_fdomain);
       serialize_output(out,d_sc_symbols);
      }

    } /* namespace kernel */
  } /* namespace filter */
} /* namespace gr */


