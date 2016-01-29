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
#include "advanced_receiver_cc_impl.h"

namespace gr {
  namespace gfdm {

    advanced_receiver_cc::sptr
    advanced_receiver_cc::make(int nsubcarrier, int ntimeslots, double filter_alpha, int fft_len, int ic_iter, gr::digital::constellation_sptr constellation, const std::string& len_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new advanced_receiver_cc_impl(nsubcarrier, ntimeslots, filter_alpha, fft_len, ic_iter, constellation, len_tag_key));
    }

    /*
     * The private constructor
     */
    advanced_receiver_cc_impl::advanced_receiver_cc_impl(int nsubcarrier, int ntimeslots, double filter_alpha, int fft_len, int ic_iter, gr::digital::constellation_sptr constellation, const std::string& len_tag_key)
      : gr::tagged_stream_block("advanced_receiver_cc",
              gr::io_signature::make(1,1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              len_tag_key),
      gfdm_receiver(nsubcarrier, ntimeslots, filter_alpha, fft_len),
      d_constellation(constellation)
    {
      d_ic_filter_taps.resize(ntimeslots);
      // Only works for d_filter_width = 2
      ::volk_32fc_x2_multiply_32fc(&d_ic_filter_taps[0],&d_filter_taps[0],&d_filter_taps[d_ntimeslots],d_ntimeslots);
    
    }

    /*
     * Our virtual destructor.
     */
    advanced_receiver_cc_impl::~advanced_receiver_cc_impl()
    {
    }

    int
    advanced_receiver_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0];
      return noutput_items ;
    }

    int
    advanced_receiver_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      filter_superposition(d_sc_fdomain,&in[0]);
      demodulate_subcarrier(d_sc_symbols,d_sc_fdomain);
      for (int j=0;j<d_ic_iter;j++)
      {
        map_sc_symbols(d_sc_symbols);
        remove_sc_interference(d_sc_symbols,d_sc_fdomain);
        //Should work since output is assigned after operation and no volk calls are in demodulate_subcarrier
        demodulate_subcarrier(d_sc_symbols,d_sc_symbols);

      }



      return noutput_items;
    }

    void
    advanced_receiver_cc_impl::map_sc_symbols( std::vector< std::vector<gr_complex> > &sc_symbols)
    {
      unsigned int symbol_tmp = 0;
      std::vector<gr_complex> const_poinst = d_constellation->points();
      for (int k=0;k<d_nsubcarrier;k++)
      {
        for (int m=0;m<d_ntimeslots;m++)
        {
          symbol_tmp =d_constellation->decision_maker(&sc_symbols[k][m]);
          sc_symbols[k][m] = const_poinst[symbol_tmp];
        }
      }
    }

    void
    advanced_receiver_cc_impl::remove_sc_interference(std::vector< std::vector<gr_complex> > &sc_symbols, std::vector< std::vector<gr_complex> > &sc_fdomain)
    {
      std::vector< std::vector<gr_complex> > prev_sc_symbols = sc_symbols;
      std::vector<gr_complex> sc_tmp(d_ntimeslots);
      for (int k=0; k<d_nsubcarrier; k++)
      {
        ::volk_32f_x2_add_32f((float*)&sc_tmp[0],(float*)&prev_sc_symbols[(k-1 % d_nsubcarrier + d_nsubcarrier) % d_nsubcarrier][0],(float*)&prev_sc_symbols[(k+1% d_nsubcarrier)][0],2*d_ntimeslots);
        ::volk_32fc_x2_multiply_32fc(&sc_symbols[k][0],&d_ic_filter_taps[0],&sc_tmp[0],d_ntimeslots);
        ::volk_32f_x2_subtract_32f((float*)&sc_tmp[0],(float*)&sc_fdomain[k][0],(float*)&sc_symbols[k][0],2*d_ntimeslots);
        ::std::memcpy(&sc_symbols[k][0],&sc_tmp[0],sizeof(gr_complex)*d_ntimeslots);

      }

    }

  } /* namespace gfdm */
} /* namespace gr */

