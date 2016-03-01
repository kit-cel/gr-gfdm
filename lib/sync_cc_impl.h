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

#ifndef INCLUDED_GFDM_SYNC_CC_IMPL_H
#define INCLUDED_GFDM_SYNC_CC_IMPL_H

#include <gfdm/sync_cc.h>
#include <volk/volk.h>
#include <pmt/pmt.h>

namespace gr {
  namespace gfdm {

    class sync_cc_impl : public sync_cc
    {
     private:
       int d_fft_len;
       int d_sync_fft_len;
       int d_cp_length;
       bool d_initialized;
       int d_block_len;
       int d_L;
       gr_complex d_autocorr_value;
       gr::gfdm::preamble_generator_sptr d_preamble_generator;
       std::vector<gr_complex> d_known_preamble;
       std::vector<float> d_P_d_abs_prev;
       std::string d_gfdm_tag_key;
       
       void initialize( const gr_complex in[] );
       void iterate( gr_complex out[], const gr_complex start[], int num_items);



     public:
      sync_cc_impl(int sync_fft_len, int cp_length, int fft_len, gr::gfdm::preamble_generator_sptr preamble_generator, const std::string& gfdm_tag_key);
      ~sync_cc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_SYNC_CC_IMPL_H */

