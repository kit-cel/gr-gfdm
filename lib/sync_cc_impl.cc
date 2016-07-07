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
    sync_cc::make(int n_subcarriers, int cp_length, int frame_len, std::vector<gr_complex> preamble,
                  const std::string &gfdm_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new sync_cc_impl(n_subcarriers, cp_length, frame_len, preamble, gfdm_tag_key));
    }

    /*
     * The private constructor
     */
    sync_cc_impl::sync_cc_impl(int n_subcarriers, int cp_length, int frame_len, std::vector<gr_complex> preamble,
                               const std::string &gfdm_tag_key)
      : gr::block("sync_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make3(1, 4, sizeof(gr_complex),sizeof(gr_complex),sizeof(float))),
      d_block_len(frame_len),
      d_gfdm_tag_key(gfdm_tag_key)

    {
      // Make sure to have only multiple of( one GFDM Block + Sync) in input
      gr::block::set_output_multiple(d_block_len + 2 * n_subcarriers);

      std::cout << "preamble: " << d_known_preamble.size() << std::endl;
      d_kernel = new improved_sync_algorithm_kernel_cc(n_subcarriers, cp_length, preamble);
    }

    /*
     * Our virtual destructor.
     */
    sync_cc_impl::~sync_cc_impl()
    {
      delete d_kernel;
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

      const int frame_pos = d_kernel->detect_frame_start(in, ninput_items[0]);
      const int avail_items = ninput_items[0] - frame_pos;

      if(avail_items < d_block_len){
        consume_each(frame_pos);
        return 0;
      }
      else{
        memcpy(out, in + frame_pos, sizeof(gr_complex) * d_block_len);
        consume_each(frame_pos + d_block_len);
        add_item_tag(0, nitems_written(0),
                     pmt::string_to_symbol(d_gfdm_tag_key),
                     pmt::from_long(d_sync_fft_len));
        return d_block_len;
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

