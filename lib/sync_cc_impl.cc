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
    sync_cc::make(int sync_fft_len, int cp_length)
    {
      return gnuradio::get_initial_sptr
        (new sync_cc_impl(sync_fft_len, cp_length));
    }

    /*
     * The private constructor
     */
    sync_cc_impl::sync_cc_impl(int sync_fft_len, int cp_length)
      : gr::block("sync_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_initialized(false),
      d_sync_fft_len(sync_fft_len),
      d_cp_length(cp_length)
    {
      set_history( 2*sync_fft_len );
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
      
      // Main task: autocorrelate L samples with previos L samples -> need to know whats the upsampling factor
      // Flow:
      // 1. autocorrelate L samples with next L samples and save value.
      // 2. multiply/conjugate L+1 with (2*L+1) and add to autocorrelation value
      // 3. multiply/conjugate 0 and L and subtract from autocorrelation value
      // (4. detect plateau{ Integrate (length CP) along previous autocorrelation values })
      for (int item=0; item<ninput_items[0]; item++)
      {
      
        if (!d_initialized)
        {
          initialize(in);
          d_initialized = true;
          out[0] = d_autocorr_value;
          item++;
          consume_each(1);
        }
        d_autocorr_value += conj(in[item+d_sync_fft_len])*in[item+2*d_sync_fft_len];
        d_autocorr_value -= conj(in[item])*in[item+d_sync_fft_len];
        out[item] = d_autocorr_value;

        consume_each (1);
      }

      // Tell runtime system how many output items we produced.
      return ninput_items[0];
    }

    void
    sync_cc_impl::initialize (const gr_complex in[])
    {
      d_autocorr_value = 0j;
      std::vector<gr_complex> corr_out(d_sync_fft_len);
      ::volk_32fc_x2_multiply_conjugate_32fc(&corr_out[0],&in[d_sync_fft_len],&in[0],d_sync_fft_len);
      for (int i=0; i<d_sync_fft_len; i++)
      {
        d_autocorr_value += corr_out[i];
      }
      
    }

  } /* namespace gfdm */
} /* namespace gr */

