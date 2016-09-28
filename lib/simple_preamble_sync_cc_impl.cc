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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "simple_preamble_sync_cc_impl.h"
#include <volk/volk.h>

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
      d_correct_cfo = false;
      d_tag_in_key = pmt::string_to_symbol(in_key);
      d_tag_out_key = pmt::string_to_symbol(out_key);
      d_tag_srcid = pmt::string_to_symbol(name());
      d_tag_value = pmt::from_long(frame_len);
      d_kernel = auto_cross_corr_multicarrier_sync_cc::sptr(
              new auto_cross_corr_multicarrier_sync_cc(subcarriers, cp_len, preamble));
      set_output_multiple(frame_len);
      set_tag_propagation_policy(TPP_DONT);
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
      int avail_items = ninput_items[0];
      int consumed_items = 0;
      noutput_items = 0;

      if(d_remaining_items > 0){ // frame from last call still present.

        if(d_correct_cfo){
          std::cout << "correct CFO: " << d_kernel->last_cfo() << std::endl;
//          remove_cfo(out, in, d_kernel->last_cfo(), d_frame_len);
          memcpy(out, in, sizeof(gr_complex) * d_frame_len);
        }
        else{
          memcpy(out, in, sizeof(gr_complex) * d_frame_len);
        }

        d_remaining_items = 0;
        noutput_items = d_frame_len;
        consumed_items = d_frame_len;
      }

      std::vector<gr::tag_t> tags;
      get_tags_in_window(tags, 0, consumed_items, avail_items, d_tag_in_key);
      if(tags.size() > 0){ // assume only one tag per call to work.
        int search_window = get_window_size_from_tag(tags[0]);
        int search_offset = get_offset_from_tag(tags[0]);

        if(search_offset + search_window < avail_items){
          int frame_start = d_kernel->detect_frame_start(in + search_offset, search_window) - d_kernel->cp_len();

          add_item_tag(0, nitems_written(0) + noutput_items, d_tag_out_key, d_tag_value, d_tag_srcid);
          d_remaining_items = d_frame_len;
          consumed_items = search_offset + frame_start;
        }
        else{
          consumed_items = search_offset;
        }
      }
      else{ // no tags, no remaining items!
        consumed_items = avail_items;
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(consumed_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    int simple_preamble_sync_cc_impl::get_offset_from_tag(const gr::tag_t& t)
    {
      return t.offset - nitems_read(0);
    }

    int simple_preamble_sync_cc_impl::get_window_size_from_tag(const gr::tag_t& t)
    {
      const int fixed_search_window = d_kernel->cp_len() + 3 * d_kernel->subcarriers();
      const int backoff = pmt::to_long(t.value) - d_frame_len;
      return fixed_search_window + backoff;
    }

    void
    simple_preamble_sync_cc_impl::remove_cfo(gr_complex* p_out, const gr_complex* p_in, const float cfo, const int ninput_size)
    {
      gr_complex initial_phase = gr_complex(1.0, 0.0);
      const float cfo_incr = -1.0f * M_PI * cfo / d_kernel->subcarriers();
      gr_complex phase_increment = gr_complex(std::cos(cfo_incr), std::sin(cfo_incr));
      volk_32fc_s32fc_x2_rotator_32fc(p_out, p_in, phase_increment, &initial_phase, ninput_size);
    }

  } /* namespace gfdm */
} /* namespace gr */

