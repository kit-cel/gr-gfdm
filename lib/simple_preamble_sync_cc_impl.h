/* -*- c++ -*- */
/* 
 * Copyright 2016 Johannes Demel.
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

#ifndef INCLUDED_GFDM_SIMPLE_PREAMBLE_SYNC_CC_IMPL_H
#define INCLUDED_GFDM_SIMPLE_PREAMBLE_SYNC_CC_IMPL_H

#include <gfdm/simple_preamble_sync_cc.h>
#include <gfdm/auto_cross_corr_multicarrier_sync_cc.h>

namespace gr {
  namespace gfdm {

    class simple_preamble_sync_cc_impl : public simple_preamble_sync_cc
    {
     private:
      int d_frame_len;
      pmt::pmt_t d_tag_in_key;
      pmt::pmt_t d_tag_out_key;
      pmt::pmt_t d_tag_srcid;
      pmt::pmt_t d_tag_value;

      auto_cross_corr_multicarrier_sync_cc::sptr d_kernel;

      int d_remaining_items;

      int d_call_to_work_count;

      int get_offset_from_tag(const gr::tag_t& t);
      int get_window_size_from_tag(const gr::tag_t& t);

     public:
      simple_preamble_sync_cc_impl(int frame_len, int subcarriers, int cp_len, std::vector<gr_complex> preamble, const std::string& in_key, const std::string& out_key);
      ~simple_preamble_sync_cc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_SIMPLE_PREAMBLE_SYNC_CC_IMPL_H */

