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
      d_gfdm_tag_key(gfdm_tag_key),
      d_is_at_frame_start(false)

    {
      // Make sure to have only multiple of( one GFDM Block + Sync) in input
      set_output_multiple(d_block_len);
      set_history(n_subcarriers + 1);

      std::cout << "preamble: " << d_known_preamble.size() << std::endl;
      d_kernel = new improved_sync_algorithm_kernel_cc(n_subcarriers, cp_length, preamble, output_multiple() + output_multiple() / 2);
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
      for (int i = 0; i < ninput_items_required.size(); ++i) {
        ninput_items_required[i] = noutput_items;
      }
    }

    int
    sync_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const int new_items_pos = history() - 1;
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      if(d_is_at_frame_start){ // this is easier and avoids errors! also, no "resync"!
        std::cout << "frame aligned call to work! produce and return!\n";
        produce_output_frame(out, in + new_items_pos);
        consume_each(d_block_len);
        d_is_at_frame_start = false;
        return d_block_len;
      }

      // if no frame is detected, kernel returns -2 * ninput_items!
      const int ninput_item_size = std::min(ninput_items[0] - new_items_pos, d_kernel->max_ninput_size());
      const int frame_pos = d_kernel->detect_frame_start(in + new_items_pos, ninput_item_size);
      const int avail_items = ninput_item_size - frame_pos;
      std::cout << "frame_pos: " << frame_pos << ", avail_items: " << avail_items << ", ninput_items: " << ninput_items[0] << ", ninput_item_size: " << ninput_item_size << ", " << ((frame_pos == -2 * ninput_item_size) ? "true" : "false") << std::endl;
      if(frame_pos == -2 * ninput_item_size){
        std::cout << "NO frame found!\n";
        consume_each(ninput_item_size);
        return 0;
      }
      else if(avail_items < d_block_len){ // align buffered items for next call to work!
        std::cout << "Not enough items! try again @next call to work!\n";
        d_is_at_frame_start = true;
        consume_each(frame_pos);
        return 0;
      }
      else{
        std::cout << "Boya! Frame FOUND! noutput_items: " << noutput_items << std::endl;
        produce_output_frame(out, in + frame_pos + new_items_pos);
        consume_each(frame_pos + d_block_len);
        return d_block_len;
      }
    }

    void sync_cc_impl::produce_output_frame(gr_complex* outbuf, const gr_complex*inbuf)
    {
      memcpy(outbuf, inbuf, sizeof(gr_complex) * d_block_len);
      add_item_tag(0, nitems_written(0),
                   pmt::string_to_symbol(d_gfdm_tag_key),
                   pmt::from_long(d_sync_fft_len));
    }

  } /* namespace gfdm */
} /* namespace gr */

