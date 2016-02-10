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
#include "cyclic_prefixer_cc_impl.h"

namespace gr {
  namespace gfdm {

    cyclic_prefixer_cc::sptr
    cyclic_prefixer_cc::make(int cp_length,const std::string& len_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new cyclic_prefixer_cc_impl(cp_length, len_tag_key));
    }

    /*
     * The private constructor
     */
    cyclic_prefixer_cc_impl::cyclic_prefixer_cc_impl(int cp_length, const std::string& len_tag_key)
      : gr::tagged_stream_block("cyclic_prefixer_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), len_tag_key),
        d_cp_length(cp_length)
    {}

    /*
     * Our virtual destructor.
     */
    cyclic_prefixer_cc_impl::~cyclic_prefixer_cc_impl()
    {
    }

    int
    cyclic_prefixer_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0] + d_cp_length;
      return noutput_items ;
    }

    int
    cyclic_prefixer_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      std::memcpy(&out[d_cp_length], &in[0], ninput_items[0]*sizeof(gr_complex));
      std::memcpy(&out[0], in+ninput_items[0]-d_cp_length, d_cp_length*sizeof(gr_complex));

      return d_cp_length+ninput_items[0];
    }

  } /* namespace gfdm */
} /* namespace gr */

