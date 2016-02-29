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
#include <pmt/pmt.h>
#include "remove_prefix_cc_impl.h"

namespace gr {
  namespace gfdm {

    remove_prefix_cc::sptr
    remove_prefix_cc::make(int sync_fft_len, int fft_len, int cp_length,const std::string& gfdm_sync_tag_key, const std::string& gfdm_len_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new remove_prefix_cc_impl(sync_fft_len, fft_len, cp_length,gfdm_sync_tag_key, gfdm_len_tag_key));
    }

    /*
     * The private constructor
     */
    remove_prefix_cc_impl::remove_prefix_cc_impl(int sync_fft_len, int fft_len, int cp_length,const std::string& gfdm_sync_tag_key, const std::string& gfdm_len_tag_key)
      : gr::block("remove_prefix_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_sync_fft_len(sync_fft_len),
      d_cp_length(cp_length),
      d_fft_len(fft_len),
      d_block_len(fft_len+sync_fft_len+2*cp_length),
      d_gfdm_sync_tag_key(gfdm_sync_tag_key),
      d_gfdm_len_tag_key(gfdm_len_tag_key),
      d_block_left(0)
    {
      set_output_multiple(d_fft_len);
      set_tag_propagation_policy(TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    remove_prefix_cc_impl::~remove_prefix_cc_impl()
    {
    }

    void
    remove_prefix_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items+d_sync_fft_len+d_cp_length*2;
    }

    int
    remove_prefix_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];
      std::vector< tag_t > sync_tags(1);
      noutput_items = 0;
      if (d_block_left>0){
        std::memcpy(&out[0],&in[0],sizeof(gr_complex)*d_block_left);
        noutput_items += d_block_left;
      }
      get_tags_in_range(sync_tags,0,nitems_read(0),nitems_read(0)+ninput_items[0]-1, pmt::string_to_symbol(d_gfdm_sync_tag_key));
      
      if (sync_tags.size() > 0)
      {
        for (std::vector<tag_t>::iterator it = sync_tags.begin();it != sync_tags.end();++it)
        {
          int sync_start = it->offset - nitems_read(0);
          int block_start = sync_start + d_sync_fft_len + d_cp_length;
          if (ninput_items[0] - block_start >= d_fft_len){
            std::memcpy(&out[noutput_items],&in[block_start],sizeof(gr_complex)*(d_fft_len));
            add_item_tag(0, nitems_written(0)+noutput_items,
                pmt::string_to_symbol(d_gfdm_len_tag_key),
                pmt::from_long(d_fft_len));
            noutput_items += d_fft_len;
            d_block_left = 0;
          }else{
            std::memcpy(&out[noutput_items],&in[block_start],sizeof(gr_complex)*(ninput_items[0]-block_start));
            add_item_tag(0, nitems_written(0)+noutput_items,
                pmt::string_to_symbol(d_gfdm_len_tag_key),
                pmt::from_long(d_fft_len));
            noutput_items += ninput_items[0] - sync_start - d_sync_fft_len - d_cp_length;
            d_block_left = d_fft_len - (ninput_items[0] - block_start); 
          }
        }
      }

      consume_each (ninput_items[0]);

      return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

