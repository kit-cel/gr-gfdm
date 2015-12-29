/* -*- c++ -*- */
/* 
 * Copyright 2015 Andrej Rode.
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
#include "transmitter_cvc_impl.h"

namespace gr {
  namespace gfdm {

    transmitter_cvc::sptr
    transmitter_cvc::make(
		    int nsubcarrier,
		    int ntimeslots,
		    int filter_width,
		    int filter_type,
		    float filter_alpha)
    {
      return gnuradio::get_initial_sptr
        (new transmitter_cvc_impl(nsubcarrier,
				 ntimeslots,
				 filter_width,
				 filter_type,
				 filter_alpha));
    }

    /*
     * The private constructor
     */
    transmitter_cvc_impl::transmitter_cvc_impl()
      : gr::block("transmitter_cvc",
              gr::io_signature::make(1,1, sizeof(gr_complex)),
              gr::io_signature::make(1,1 , sizeof(gr_complex)))
    {}

    /*
     * Our virtual destructor.
     */
    transmitter_cvc_impl::~transmitter_cvc_impl()
    {
    }

    void
    transmitter_cvc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
    }

    int
    transmitter_cvc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const <+ITYPE*> *in = (const <+ITYPE*> *) input_items[0];
        <+OTYPE*> *out = (<+OTYPE*> *) output_items[0];

        // Do <+signal processing+>
        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

