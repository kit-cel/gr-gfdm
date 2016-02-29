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

#ifndef INCLUDED_GFDM_REMOVE_PREFIX_CC_IMPL_H
#define INCLUDED_GFDM_REMOVE_PREFIX_CC_IMPL_H

#include <gfdm/remove_prefix_cc.h>

namespace gr {
  namespace gfdm {

    class remove_prefix_cc_impl : public remove_prefix_cc
    {
     private:
       int d_fft_len;
       int d_sync_fft_len;
       int d_cp_length;
       int d_block_len;
       int d_block_left;
       const std::string& d_gfdm_sync_tag_key;
       const std::string& d_gfdm_len_tag_key;

     public:
      remove_prefix_cc_impl(int sync_fft_len, int fft_len, int cp_length, const std::string& gfdm_sync_tag_key, const std::string& gfdm_len_tag_key);
      ~remove_prefix_cc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_REMOVE_PREFIX_CC_IMPL_H */

