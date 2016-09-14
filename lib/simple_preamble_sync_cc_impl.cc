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
#include "simple_preamble_sync_cc_impl.h"

namespace gr
{
  namespace gfdm
  {

    simple_preamble_sync_cc::sptr
    simple_preamble_sync_cc::make(int frame_len, int subcarriers, int cp_len, std::vector<gr_complex> preamble,
                                  const std::string &in_key, const std::string &out_key)
    {
      return gnuradio::get_initial_sptr
              (new simple_preamble_sync_cc_impl(frame_len, subcarriers, cp_len, preamble, in_key, out_key));
    }

    /*
     * The private constructor
     */
    simple_preamble_sync_cc_impl::simple_preamble_sync_cc_impl(int frame_len, int subcarriers, int cp_len,
                                                               std::vector<gr_complex> preamble,
                                                               const std::string &in_key, const std::string &out_key)
            : gr::block("simple_preamble_sync_cc",
                        gr::io_signature::make(1, 1, sizeof(gr_complex)),
                        gr::io_signature::make(1, 1, sizeof(gr_complex))), d_frame_len(frame_len), d_remaining_items(0)
    {
      d_tag_in_key = pmt::string_to_symbol(in_key);
      d_tag_out_key = pmt::string_to_symbol(out_key);
      d_tag_srcid = pmt::string_to_symbol(name());
      d_tag_value = pmt::from_long(frame_len);
      d_kernel = auto_cross_corr_multicarrier_sync_cc::sptr(
              new auto_cross_corr_multicarrier_sync_cc(subcarriers, cp_len, preamble));
      set_output_multiple(frame_len);
    }

    /*
     * Our virtual destructor.
     */
    simple_preamble_sync_cc_impl::~simple_preamble_sync_cc_impl()
    {
    }

    void
    simple_preamble_sync_cc_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    int
    simple_preamble_sync_cc_impl::general_work(int noutput_items,
                                               gr_vector_int &ninput_items,
                                               gr_vector_const_void_star &input_items,
                                               gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      std::vector<gr::tag_t> tags;
      get_tags_in_window(tags, 0, 0, ninput_items[0], d_tag_in_key);
      int consumed_items = noutput_items;

      if(tags.size() > 0){ // assume we mostly only receive one tag at most!
        gr::tag_t tag = tags[0];
        const int fixed_search_window = d_kernel->cp_len() + 3 * d_kernel->subcarriers();
        int backoff = pmt::to_long(tag.value) - d_frame_len;
        int search_window = fixed_search_window + backoff;
        int search_offset = tag.offset - nitems_read(0);
        if(search_offset + search_window < ninput_items[0]){
          int frame_start = d_kernel->detect_frame_start(in + search_offset, search_window);
//          std::cout << backoff << ", " << search_window <<  ", " << tag.offset << ", call offset: "
//          << search_offset << ", " << noutput_items << ", NOW SEARCH " << frame_start << std::endl;
          add_item_tag(0, nitems_written(0) + search_offset + frame_start, d_tag_out_key, d_tag_value, d_tag_srcid);

        }
        else{
//          std::cout << backoff << ", " << search_window <<  ", " << tag.offset << ", call offset: "
//          << search_offset << ", " << noutput_items << ", NEXT" << std::endl;
          consumed_items = search_offset;
          noutput_items = search_offset; // change to 0 later!
        }



      }
      else{
        // set noutput items to 0 and consume all input items!
      }


      memcpy(out, in, sizeof(gr_complex) * noutput_items);

      // Do <+signal processing+>
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(consumed_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

