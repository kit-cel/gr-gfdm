/* -*- c++ -*- */
/* 
 * Copyright 2017 <+YOU OR YOUR COMPANY+>.
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
#include "extract_burst_cc_impl.h"
#include <cmath>

namespace gr {
  namespace gfdm {

    extract_burst_cc::sptr
    extract_burst_cc::make(int burst_len, std::string burst_start_tag)
    {
      return gnuradio::get_initial_sptr
        (new extract_burst_cc_impl(burst_len, burst_start_tag));
    }

    /*
     * The private constructor
     */
    extract_burst_cc_impl::extract_burst_cc_impl(int burst_len, std::string burst_start_tag)
      : gr::block("extract_burst_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
              d_burst_len(burst_len)
    {
      set_output_multiple(burst_len);
      d_burst_start_tag = pmt::string_to_symbol(burst_start_tag);
      set_tag_propagation_policy(TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    extract_burst_cc_impl::~extract_burst_cc_impl()
    {
    }

    void
    extract_burst_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
        ninput_items_required[0] = noutput_items;
    }

    int
    extract_burst_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      int n_out_bursts = noutput_items / d_burst_len;
      int avail_items = ninput_items[0];
      int consumed_items = avail_items;
      int produced_items = 0;

      std::vector<tag_t> tags;
      get_tags_in_window(tags, 0, 0, avail_items, d_burst_start_tag);
      const int n_max_bursts = std::min(int(tags.size()), n_out_bursts);
      for (int i = 0; i < n_max_bursts; ++i) {
        int burst_start = tags[i].offset - nitems_read(0);
        if(avail_items - burst_start >= d_burst_len){
          memcpy(out, in + burst_start, sizeof(gr_complex) * d_burst_len);

          add_item_tag(0, nitems_written(0) + produced_items, d_burst_start_tag,
                       pmt::from_long(d_burst_len), pmt::string_to_symbol(name()));

          produced_items += d_burst_len;
          consumed_items = burst_start + d_burst_len;
          out += d_burst_len;
        }
        else{
          consumed_items = burst_start;
          break;
        }
      }

      consume_each (consumed_items);
      return produced_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

