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
#include "modulator_cc_impl.h"

namespace gr {
  namespace gfdm {

    modulator_cc::sptr
    modulator_cc::make(
        const std::string& len_tag_key,
        int nsubcarrier,
        int ntimeslots,
        int filter_width,
        double filter_alpha)
    {
      return gnuradio::get_initial_sptr
        (new modulator_cc_impl(len_tag_key,
                               nsubcarrier,
                               ntimeslots,
                               filter_width,
                               filter_alpha)
         );

    }

    /*
     * The private constructor
     */
    modulator_cc_impl::modulator_cc_impl(
        const std::string& len_tag_key,
        int nsubcarrier,
        int ntimeslots,
        int filter_width,
        double filter_alpha)
      : gr::tagged_stream_block("modulator_cc",
              gr::io_signature::make(nsubcarrier*ntimeslots, nsubcarrier*ntimeslots, sizeof(gr_complex)),
              gr::io_signature::make(nsubcarrier*ntimeslots, nsubcarrier*ntimeslots, sizeof(gr_complex)),
              len_tag_key),
      d_ntimeslots(ntimeslots),
      d_nsubcarrier(nsubcarrier),
      d_filter_width(filter_width)
    {}

    /*
     * Our virtual destructor.
     */
    modulator_cc_impl::~modulator_cc_impl()
    {
    }

    int
    modulator_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0];
      return noutput_items ;
    }

    int
    modulator_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];

        // Do signal processing!


        return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

