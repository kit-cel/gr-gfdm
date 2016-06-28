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
#include "simple_receiver_cc_impl.h"

namespace gr {
  namespace gfdm {

    simple_receiver_cc::sptr
    simple_receiver_cc::make(int n_timeslots, int n_subcarriers, int overlap, std::vector<gr_complex> frequency_taps)
    {
      return gnuradio::get_initial_sptr
        (new simple_receiver_cc_impl(n_timeslots, n_subcarriers, overlap, frequency_taps));
    }

    /*
     * The private constructor
     */
    simple_receiver_cc_impl::simple_receiver_cc_impl(int n_timeslots, int n_subcarriers, int overlap, std::vector<gr_complex> frequency_taps)
      : gr::sync_block("simple_receiver_cc",
              gr::io_signature::make(<+MIN_IN+>, <+MAX_IN+>, sizeof(<+ITYPE+>)),
              gr::io_signature::make(<+MIN_OUT+>, <+MAX_OUT+>, sizeof(<+OTYPE+>)))
    {}

    /*
     * Our virtual destructor.
     */
    simple_receiver_cc_impl::~simple_receiver_cc_impl()
    {
    }

    int
    simple_receiver_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const <+ITYPE+> *in = (const <+ITYPE+> *) input_items[0];
      <+OTYPE+> *out = (<+OTYPE+> *) output_items[0];

      // Do <+signal processing+>

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

