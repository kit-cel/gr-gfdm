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
    cyclic_prefixer_cc::make(int block_len, int cp_len, int cs_len, int ramp_len, std::vector<gr_complex> window_taps)
    {
      return gnuradio::get_initial_sptr
              (new cyclic_prefixer_cc_impl(block_len, cp_len, cs_len, ramp_len, window_taps));
    }

    /*
     * The private constructor
     */
    cyclic_prefixer_cc_impl::cyclic_prefixer_cc_impl(int block_len, int cp_len, int cs_len, int ramp_len,
                                                     std::vector<gr_complex> window_taps)
            : gr::block("cyclic_prefixer_cc",
                                      gr::io_signature::make(1, 1, sizeof(gr_complex)),
                                      gr::io_signature::make(1, 1, sizeof(gr_complex))),
              d_cp_len(cp_len)
    {
      // all the work is done in the kernel!
      d_kernel = add_cyclic_prefix_cc::sptr(
              new add_cyclic_prefix_cc(block_len, cp_len, cs_len, ramp_len, window_taps));

      // set block properties!
      set_relative_rate(1.0 * d_kernel->frame_size() / d_kernel->block_size());
      set_fixed_rate(true);
      set_output_multiple(d_kernel->frame_size());
    }

    /*
     * Our virtual destructor.
     */
    cyclic_prefixer_cc_impl::~cyclic_prefixer_cc_impl()
    {
    }

    void
    cyclic_prefixer_cc_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
    {
      for (int i = 0; i < ninput_items_required.size(); ++i) {
        ninput_items_required[i] = fixed_rate_noutput_to_ninput(noutput_items);
      }
    }

    int
    cyclic_prefixer_cc_impl::fixed_rate_ninput_to_noutput(int ninput)
    {
      return (ninput / d_kernel->block_size()) * d_kernel->frame_size();;
    }

    int
    cyclic_prefixer_cc_impl::fixed_rate_noutput_to_ninput(int noutput)
    {
      return (noutput / d_kernel->frame_size()) * d_kernel->block_size();;
    }

    int
    cyclic_prefixer_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];
      const int n_frames = noutput_items / d_kernel->frame_size();

      for (int i = 0; i < n_frames; ++i) {
        d_kernel->generic_work(out, in);
        in += d_kernel->block_size();
        out += d_kernel->frame_size();
      }

      consume_each(n_frames * d_kernel->block_size());
      return n_frames * d_kernel->frame_size();
    }

  } /* namespace gfdm */
} /* namespace gr */

