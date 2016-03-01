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
#include "sync_cc_impl.h"

namespace gr {
  namespace gfdm {

    sync_cc::sptr
    sync_cc::make(int sync_fft_len, int cp_length, int fft_len, gr::gfdm::preamble_generator_sptr preamble_generator, const std::string& gfdm_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new sync_cc_impl(sync_fft_len, cp_length, fft_len, preamble_generator, gfdm_tag_key));
    }

    /*
     * The private constructor
     */
    sync_cc_impl::sync_cc_impl(int sync_fft_len, int cp_length, int fft_len, gr::gfdm::preamble_generator_sptr preamble_generator, const std::string& gfdm_tag_key)
      : gr::block("sync_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make3(1, 4, sizeof(gr_complex),sizeof(gr_complex),sizeof(float))),
      d_initialized(false),
      d_fft_len(fft_len),
      d_sync_fft_len(sync_fft_len),
      d_cp_length(cp_length),
      d_block_len(2*cp_length+fft_len+sync_fft_len),
      d_L(sync_fft_len/2),
      d_preamble_generator(preamble_generator),
      d_gfdm_tag_key(gfdm_tag_key)
    {
      set_tag_propagation_policy(TPP_DONT);
      set_history(2);
      // Make sure to have only multiple of( one GFDM Block + Sync) in input
      gr::block::set_output_multiple( d_block_len+d_sync_fft_len );
      d_P_d_abs_prev.resize(cp_length,0);
      d_known_preamble = d_preamble_generator->get_preamble();
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
      gr_complex *corr_out;
      float *corr_i_out;
      float *res_out;

      //Initialize some vectors to hold Correlation_data
      //P_d: (complex) autocorrelation of length sync_fft_len/2 to detect signal with two identical halves length sync_fft_len
      //P_d_abs: absolute value
      //P_d_i: integrated P_d_abs -cp_length:0 to eliminate CP Plateau
      std::vector<gr_complex> P_d(d_block_len);
      //Add last (length cp_length) samples from P_d_abs to start (to integrate cp_length samples)
      std::vector<float> P_d_abs(d_block_len+d_cp_length);
      std::vector<float> P_d_i(d_block_len);
      
      // Main task: autocorrelate L samples with following L samples -> need to know whats the upsampling factor
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
      }
      iterate(&P_d[0], &in[0], d_block_len);
      
      //Now integrate along time axis (length cp)
      if (d_cp_length)
      {
        ::volk_32fc_magnitude_32f(&P_d_abs[d_cp_length],&P_d[0],d_block_len);
        std::memcpy(&P_d_abs[0],&d_P_d_abs_prev[0],sizeof(float)*d_cp_length);
        std::memcpy(&d_P_d_abs_prev[0],&P_d_abs[d_block_len],sizeof(float)*d_cp_length);
        for (int i=0;i<(d_block_len);i++)
        {
          for (int k=0;k<d_cp_length;k++)
          {
            P_d_i[i] += P_d_abs[i+k];
          }
        }
      }else
      {
        ::volk_32fc_magnitude_32f(&P_d_abs[0], &P_d[0],d_block_len);
        std::memcpy(&P_d_i[0],&P_d_abs[0],sizeof(float)*d_block_len);
      }

      
      //Find max abs
      std::vector<float>::iterator max;
      max = std::max_element(P_d_i.begin(),P_d_i.end());
      int max_index1 = std::distance(P_d_i.begin(),max);
      //Calculate angle <P_d
      float angle =(float) std::atan2((double) std::imag(P_d[max_index1]), (double) std::real(P_d[max_index1]));
      float cfo = angle/M_PI;

      //std::cout << "Carrier Frequency Offset: " <<cfo<<std::endl;
      
      //Correct CFO with epsilon = angle/pi
      std::vector<gr_complex> cfo_correction(d_block_len+d_sync_fft_len);
      for (int i=0;i<(d_block_len+d_sync_fft_len);i++)
      {
        cfo_correction[i] = std::exp( gr_complex(2j*M_PI*(cfo/(d_sync_fft_len/2))*i));
      }
      std::vector<gr_complex> corrected_input_sequence(d_block_len+d_sync_fft_len);
      ::volk_32fc_x2_multiply_32fc(&corrected_input_sequence[0],&cfo_correction[0],&in[1],d_block_len+d_sync_fft_len);
      //Crosscorrelate known preamble and sync_preamble
      //Known Preamble must have length d_sync_fft_len
      std::vector<gr_complex> cross_correlation(d_block_len);
      for (int i=0; i<(d_block_len);i++)
      {
        std::vector<gr_complex> cc_tmp(d_sync_fft_len);
        ::volk_32fc_x2_multiply_conjugate_32fc(&cc_tmp[0],&d_known_preamble[0],&corrected_input_sequence[i],d_sync_fft_len);
        for (int k=0; k<d_sync_fft_len;k++)
        {
          cross_correlation[i] += (gr_complex) (((gr_complex)1.0) / ((gr_complex) d_sync_fft_len)) * cc_tmp[k];
        }
      }
      //multiply crosscorelation with P_d (autocorrelation)
      std::vector<float> cc_abs(d_block_len);
      ::volk_32fc_magnitude_32f(&cc_abs[0],&cross_correlation[0],d_block_len);
      std::vector<float> P_d_res(d_block_len);
      ::volk_32f_x2_multiply_32f(&P_d_res[0],&cc_abs[0],&P_d_abs[0],d_block_len);

      max = std::max_element(P_d_res.begin(),P_d_res.end());
      int max_index2 = std::distance(P_d_res.begin(),max);

      //Add evaluation with threshold and multipath detection (argfirst)
      //Add Stream tags on max
      
      //std::cout << "First Autocorrelation maximum: " << max_index1 <<std::endl;
      //std::cout << "CC maximum: " << max_index2 <<std::endl;
      
      add_item_tag(0, nitems_written(0)+max_index2,
          pmt::string_to_symbol(d_gfdm_tag_key),
          pmt::from_long(d_sync_fft_len));

      //in[0] is last item of previous block
      std::memcpy(&out[0],&in[1],sizeof(gr_complex)*d_block_len);
      // Copy P_d into second (float) port
      if (output_items.size()>1)
      {
        corr_out = (gr_complex *) output_items[1];
        std::memcpy(&corr_out[0],&corrected_input_sequence[0],sizeof(gr_complex)*d_block_len);
        add_item_tag(1, nitems_written(0)+max_index2,
          pmt::string_to_symbol(d_gfdm_tag_key),
          pmt::from_long(d_sync_fft_len));
      }
      if (output_items.size()>2)
      {
        corr_i_out = (float *) output_items[2];
        std::memcpy(&corr_i_out[0],&P_d_i[0],sizeof(float)*d_block_len);
        add_item_tag(2, nitems_written(0)+max_index1,
            pmt::string_to_symbol(d_gfdm_tag_key),
            pmt::from_long(d_sync_fft_len));       
      }
      if (output_items.size()>3)
      {
        res_out = (float *) output_items[3];
        std::memcpy(&res_out[0], &P_d_res[0], sizeof(float)*d_block_len);
      }
      consume_each(d_block_len);
      
      return d_block_len;
    }

    void
    sync_cc_impl::initialize (const gr_complex in[])
    {
      d_autocorr_value = 0j;
      std::vector<gr_complex> corr_out(d_L);
      ::volk_32fc_x2_multiply_conjugate_32fc(&corr_out[0],&in[d_L],&in[0],d_L);
      for (int i=0; i<d_L; i++)
      {
        d_autocorr_value += corr_out[i];
      }
      d_autocorr_value *= (gr_complex) (1.0/d_L);
      
    }

    void
    sync_cc_impl::iterate( gr_complex out[], const gr_complex start[], int num_items)
    {
      for (int i=0; i<num_items;i++)
      {
        d_autocorr_value -=((gr_complex) (1.0/d_L))*conj(start[i])*start[i+d_L];
        d_autocorr_value +=((gr_complex) (1.0/d_L))*conj(start[d_L+i])*start[2*d_L+i];
        out[i] = d_autocorr_value;
       
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

