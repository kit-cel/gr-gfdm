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
#include "frame_energy_detector_cc_impl.h"

namespace gr
{
  namespace gfdm
  {

    frame_energy_detector_cc::sptr
    frame_energy_detector_cc::make(float alpha, int average_len, int frame_len, int backoff_len, const std::string& tag_key)
    {
      return gnuradio::get_initial_sptr
              (new frame_energy_detector_cc_impl(alpha, average_len, frame_len, backoff_len, tag_key));
    }

    /*
     * The private constructor
     */
    frame_energy_detector_cc_impl::frame_energy_detector_cc_impl(float alpha, int average_len, int frame_len,
                                                                 int backoff_len, const std::string& tag_key)
            : gr::block("frame_energy_detector_cc",
                        gr::io_signature::make(1, 1, sizeof(gr_complex)),
                        gr::io_signature::make(1, 1, sizeof(gr_complex))),
              d_frame_len(frame_len), d_backoff_len(backoff_len), d_remaining_frame(0), d_detector_stage(SEARCHING)
    {
      d_tag_key = pmt::string_to_symbol(tag_key);
      d_tag_srcid = pmt::string_to_symbol(name());
      d_tag_value = pmt::from_long(d_frame_len + (d_backoff_len > -1 ? 2 * d_backoff_len : 0));
      d_kernel = detect_frame_energy_kernel_cl::sptr(new detect_frame_energy_kernel_cl(alpha, average_len));
    }

    /*
     * Our virtual destructor.
     */
    frame_energy_detector_cc_impl::~frame_energy_detector_cc_impl()
    {
    }

    void
    frame_energy_detector_cc_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = (noutput_items < d_kernel->average_len()) ? d_kernel->average_len() : noutput_items;
    }

    int
    frame_energy_detector_cc_impl::general_work(int noutput_items,
                                                gr_vector_int &ninput_items,
                                                gr_vector_const_void_star &input_items,
                                                gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

//      std::cout << "noutputitems: " << noutput_items << ", ninput_items: " << ninput_items[0] << std::endl;

      int available_items = d_kernel->average_len() * (ninput_items[0] / d_kernel->average_len());
      int consumed_items = 0;
      int produced_items = 0;


      if(d_remaining_frame > 0){
        int out_items = std::min(d_remaining_frame, noutput_items);
        std::cout << "frame_stage! " << out_items << std::endl;
        memcpy(out, in, sizeof(gr_complex) * out_items);
        out += out_items;
        in += out_items;
        d_remaining_frame -= out_items;
        consumed_items += out_items;
        produced_items += out_items;
        available_items = d_kernel->average_len() * ((available_items - out_items) / d_kernel->average_len());
      }

      if(available_items > 0){
        int frame_pos = d_kernel->detect_frame(in, available_items);
        if(frame_pos > -1){
          int backoff_start = 0;
          if(d_backoff_len > -1){ // only output detected frames!
            backoff_start = std::max(frame_pos - d_backoff_len, 0);
            add_item_tag(0, nitems_written(0) + produced_items, d_tag_key, d_tag_value, d_tag_srcid);
          }
          else{
            add_item_tag(0, nitems_written(0) + produced_items + frame_pos, d_tag_key, d_tag_value, d_tag_srcid);
          }

          d_remaining_frame = d_frame_len + 2 * d_backoff_len;
          int avail_items = std::min(d_remaining_frame, available_items - backoff_start);
          memcpy(out, in + backoff_start, sizeof(gr_complex) * avail_items);



          d_remaining_frame -= avail_items;
          produced_items += avail_items;
        }
        if(d_backoff_len < 0){
          memcpy(out, in, sizeof(gr_complex) * noutput_items);
        }
        consumed_items += available_items;
      }



//      if(d_remaining_frame > 0){
//        int avail_items = std::min(d_remaining_frame, noutput_items);
//        std::cout << "frame_stage! " << avail_items << std::endl;
//        memcpy(out, in, sizeof(gr_complex) * avail_items);
//        d_remaining_frame -= avail_items;
//        consumed_items = avail_items;
//        noutput_items = avail_items;
//      }
//      else{
//        int frame_pos = d_kernel->detect_frame(in, noutput_items);
//
//        if(d_backoff_len < 0){
//          memcpy(out, in, sizeof(gr_complex) * noutput_items);
//          if(frame_pos > -1){
//            add_item_tag(0, nitems_written(0) + frame_pos, d_tag_key, d_tag_value, d_tag_srcid);
//          }
//        }
//        else{
//          if(frame_pos > -1){
//            int backoff_start = std::max(frame_pos - d_backoff_len, 0);
//            int avail_items = std::min(d_frame_len + 2 * d_backoff_len, noutput_items - backoff_start);
//            memcpy(out, in + backoff_start, sizeof(gr_complex) * avail_items);
//            add_item_tag(0, nitems_written(0), d_tag_key, d_tag_value, d_tag_srcid);
//            d_remaining_frame = d_frame_len + 2 * d_backoff_len - avail_items;
//            consumed_items = backoff_start + avail_items;
//            noutput_items = avail_items;
//          }
//          else{
//            noutput_items = 0;
//          }
//        }
//      }


      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(consumed_items);

      // Tell runtime system how many output items we produced.
      return produced_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

