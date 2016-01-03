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
#include <gnuradio/fft/fft.h>
#include <gnuradio/filter/firdes.h>

namespace gr {
  namespace gfdm {

    transmitter_cvc::sptr
    transmitter_cvc::make(
		    int nsubcarrier,
		    int ntimeslots,
		    int filter_width,
		    double filter_alpha)
    {
      return gnuradio::get_initial_sptr
        (new transmitter_cvc_impl(nsubcarrier,
				 ntimeslots,
				 filter_width,
				 filter_alpha));
    }

    /*
     * The private constructor
     */
    transmitter_cvc_impl::transmitter_cvc_impl(
		    int nsubcarrier,
		    int ntimeslots,
		    int filter_width,
		    double filter_alpha)
      : gr::block("transmitter_cvc",
              gr::io_signature::make(1,-1, sizeof(gr_complex)),
              gr::io_signature::make(1,1 , sizeof(gr_complex) * nsubcarrier * ntimeslots)),
        d_nsubcarrier(nsubcarrier),
        d_ntimeslots(ntimeslots),
        d_filter_width(filter_width)

    {
	    double sampling_freq = 1.0/nsubcarrier;
	    std::vector<float> filtertaps = gr::filter::firdes::root_raised_cosine(
			    1.0,
			    1.0,
			    sampling_freq,
			    filter_alpha,
			    d_nsubcarrier * d_ntimeslots);
             
            
            fft::fft_real_fwd *filter_fft = new fft::fft_real_fwd(d_nsubcarrier*d_ntimeslots,1);
            float *in = filter_fft->get_inbuf();
            gr_complex *out = filter_fft->get_outbuf();
            for (int i=0;i<d_nsubcarrier*d_ntimeslots-1;i++)
            {
              in[i] = filtertaps[i];
            }
            filter_fft->execute();
            // filter_width*d_ntimeslots must be power of 2
            int filter_length = (filter_width*d_ntimeslots)/2;
            for (int i=0;i<filter_length;i++)
            {
              d_filtertaps.push_back(out[i]);
            }
            for (int i=0;i<filter_width;i++)
            {
              d_filtertaps.push_back(out[(d_nsubcarrier*d_ntimeslots-filter_length)+i]);
            }
            delete filter_fft;

    }

    /*
     * Our virtual destructor.
     */
    transmitter_cvc_impl::~transmitter_cvc_impl()
    {
    }

    void
    transmitter_cvc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = d_nsubcarrier * d_ntimeslots * noutput_items;
    }

    int
    transmitter_cvc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];

        memset((void *) out, 0x00, sizeof(gr_complex) * d_nsubcarrier * d_ntimeslots * noutput_items);
        fft::fft_complex *subcarrier_fft = new fft::fft_complex(d_ntimeslots,1,1);
        gr_complex *fft_in = subcarrier_fft->get_inbuf();
        gr_complex *fft_out = subcarrier_fft->get_outbuf();
        for (int c = 0; c < d_nsubcarrier;c++)
        {
          for (int t = 0; t < d_ntimeslots ;t++)
          {
            fft_in[t] = in[c+d_nsubcarrier*t];
          }
          subcarrier_fft->execute();
          // Filter output (tile with factor filter_width) with filter_width * ntimeslots length filtersamples
          for (int t =0 ; t < d_filter_width*d_ntimeslots;t++)
          {
            fft_out[t]*d_filtertaps[t];
          }
        }

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

