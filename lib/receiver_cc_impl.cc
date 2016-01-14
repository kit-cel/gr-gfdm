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
#include "receiver_cc_impl.h"

namespace gr {
  namespace gfdm {

    receiver_cc::sptr
    receiver_cc::make(int nsubcarrier,
                      int ntimeslots,
                      int filter_width,
                      double filter_alpha)
    {
      return gnuradio::get_initial_sptr
        (new receiver_cc_impl(int nsubcarrier,int ntimeslots,int filter_width, double filter_alpha));
    }

    /*
     * The private constructor
     */
    receiver_cc_impl::receiver_cc_impl(int nsubcarrier, int ntimeslots, int filter_width, double filter_alpha)
      : gr::block("receiver_cc",
              gr::io_signature::make(nsubcarrier*ntimeslots, nsubcarrier*ntimeslots , sizeof(gr_complex)),
              gr::io_signature::make(nsubcarrier*ntimeslots, nsubcarrier*ntimeslots, sizeof(gr_complex))),
      d_nsubcarrier(nsubcarrier),
      d_ntimeslots(ntimeslots)
    {}

    /*
     * Our virtual destructor.
     */
    receiver_cc_impl::~receiver_cc_impl()
    {
    }

    void
    receiver_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
    }

    int
    receiver_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];

        // Do <+signal processing+>
        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

