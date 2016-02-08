/* -*- c++ -*- */
/* 
 * Copyright 2016 <+YOU OR YOUR COMPANY+>.
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
#include "sync_cc_impl.h"

namespace gr {
  namespace gfdm {

    sync_cc::sptr
    sync_cc::make(int sync_fft_len, int cp_length, int fft_len)
    {
      return gnuradio::get_initial_sptr
        (new sync_cc_impl(sync_fft_len, cp_length, fft_len));
    }

    /*
     * The private constructor
     */
    sync_cc_impl::sync_cc_impl(int sync_fft_len, int cp_length, int fft_len)
      : gr::block("sync_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_initialized(false),
      d_sync_fft_len(sync_fft_len),
      d_cp_length(cp_length),
      d_block_len(2*cp_length+fft_len+sync_fft_len),
      d_L(sync_fft_len/2)
    {
      set_history(2);
      // Make sure to have only multiple of( one GFDM Block + Sync) in input
      gr::block::set_output_multiple( d_block_len+d_sync_fft_len );
    }

    /*
     * Our virtual destructor.
     */
    sync_cc_impl::~sync_cc_impl()
    {
    }

    void
    sync_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    int
    sync_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      std::vector<gr_complex> P_d(d_block_len);
      
      // Main task: autocorrelate L samples with previos L samples -> need to know whats the upsampling factor
      // Flow:
      // 1. autocorrelate L samples with next L samples and save value.
      // 2. multiply/conjugate L+1 with (2*L+1) and add to autocorrelation value
      // 3. multiply/conjugate 0 and L and subtract from autocorrelation value
      // (4. detect plateau{ Integrate (length CP) along previous autocorrelation values })
      // (5. match against threshold (detect max/first peak -> multipath))
      //
      // We have sync_fft_len + H*(2*cp_length+fft_len+sync_fft_len+sync_fft_len/2) items
      if (!d_initialized)
      {
        initialize(in); //Initial Value in 
        d_initialized = true;
        P_d[0] = d_autocorr_value;
        //P_d[0] corresponds to in[0]
        iterate(&P_d[1], &in[0], d_block_len-1);
        std::memcpy(&out[0],&in[0],sizeof(gr_complex)*d_block_len);
      } else
      {
        //P_d[0] corresponds to in[1] due to history?
        iterate(&P_d[0], &in[0], d_block_len);
        std::memcpy(&out[0],&in[1],sizeof(gr_complex)*d_block_len);
      }

      consume_each(d_block_len);
      
      return d_block_len;
    }

    void
    sync_cc_impl::initialize (const gr_complex in[])
    {
      d_autocorr_value = 0j;
      std::vector<gr_complex> corr_out(d_sync_fft_len);
      ::volk_32fc_x2_multiply_conjugate_32fc(&corr_out[0],&in[d_L],&in[0],d_L);
      for (int i=0; i<d_sync_fft_len; i++)
      {
        d_autocorr_value += corr_out[i];
      }
      
    }

    void
    sync_cc_impl::iterate( gr_complex out[], const gr_complex start[], int num_items)
    {
      for (int i=0; i<num_items;i++)
      {
        d_autocorr_value -= conj(start[i])*start[i+d_L];
        d_autocorr_value += conj(start[d_L+i])*start[2*d_L+i];
        out[i] = d_autocorr_value;
       
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

