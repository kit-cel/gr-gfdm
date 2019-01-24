/* -*- c++ -*- */
/*
 * Copyright 2019 Johannes Demel.
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

#include <exception>
#include <cstring>
#include <gnuradio/io_signature.h>
#include "short_burst_shaper_impl.h"

namespace gr {
  namespace gfdm {

    short_burst_shaper::sptr
    short_burst_shaper::make(int pre_padding, int post_padding,
                             const std::string &length_tag_name)
    {
      return gnuradio::get_initial_sptr
        (new short_burst_shaper_impl(pre_padding, post_padding, length_tag_name));
    }

    /*
     * The private constructor
     */
    short_burst_shaper_impl::short_burst_shaper_impl(int pre_padding, int post_padding,
                                                     const std::string &length_tag_name)
      : gr::tagged_stream_block("short_burst_shaper",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), length_tag_name),
      d_pre_padding(pre_padding), d_post_padding(post_padding)
    {
      if(d_pre_padding < 0){
        throw std::invalid_argument("Pre-padding length MUST be >= 0!");
      }
      if(d_post_padding < 0){
        throw std::invalid_argument("Post-padding length MUST be >= 0!");
      }
    }

    /*
     * Our virtual destructor.
     */
    short_burst_shaper_impl::~short_burst_shaper_impl()
    {
    }

    int
    short_burst_shaper_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0] + d_pre_padding + d_post_padding;
      return noutput_items;
    }

    int
    short_burst_shaper_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      // std::cout << "noutput_items = " << noutput_items << std::endl;
      // std::cout << "ninput_items = " << ninput_items[0] << std::endl;
      // std::cout << "pre_padding = " << d_pre_padding << std::endl;
      // std::cout << "post_padding = " << d_post_padding << std::endl;

      std::memset(out, 0, sizeof(gr_complex) * d_pre_padding);
      std::memcpy(out + d_pre_padding, in, sizeof(gr_complex) * ninput_items[0]);
      std::memset(out + d_pre_padding + ninput_items[0], 0,
                  sizeof(gr_complex) * d_post_padding);

      // Tell runtime system how many output items we produced.
      return ninput_items[0] + d_pre_padding + d_post_padding;
    }

  } /* namespace gfdm */
} /* namespace gr */

