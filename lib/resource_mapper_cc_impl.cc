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
#include "resource_mapper_cc_impl.h"

namespace gr {
  namespace gfdm {

    resource_mapper_cc::sptr
    resource_mapper_cc::make(int active_subcarriers, int fft_len, int timeslots, std::vector<int> subcarrier_map, bool per_timeslot)
    {
      return gnuradio::get_initial_sptr
        (new resource_mapper_cc_impl(active_subcarriers, fft_len, timeslots, subcarrier_map, per_timeslot));
    }

    /*
     * The private constructor
     */
    resource_mapper_cc_impl::resource_mapper_cc_impl(int timeslots, int subcarriers, int active_subcarriers, std::vector<int> subcarrier_map, bool per_timeslot)
      : gr::block("resource_mapper_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
      d_kernel = resource_mapper_kernel_cc::sptr(new resource_mapper_kernel_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, per_timeslot));

      set_relative_rate(1.0 * d_kernel->input_vector_size() / d_kernel->output_vector_size());
      set_fixed_rate(true);
      set_output_multiple(d_kernel->output_vector_size());
      }

    /*
     * Our virtual destructor.
     */
    resource_mapper_cc_impl::~resource_mapper_cc_impl()
    {
    }

    void
    resource_mapper_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
      for (int i = 0; i < ninput_items_required.size(); ++i) {
        ninput_items_required[i] = fixed_rate_noutput_to_ninput(noutput_items);
      }
    }

    int
    resource_mapper_cc_impl::fixed_rate_ninput_to_noutput(int ninput)
    {
      return (ninput / d_kernel->input_vector_size()) * d_kernel->output_vector_size();
    }

    int
    resource_mapper_cc_impl::fixed_rate_noutput_to_ninput(int noutput)
    {
      return (noutput / d_kernel->output_vector_size()) * d_kernel->input_vector_size();
    }

    int
    resource_mapper_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];
      // The gr block uses a less flexible version of the kernel. Zero-padding is ignored.
      int n_frames = std::min(noutput_items / d_kernel->output_vector_size(), ninput_items[0] / d_kernel->input_vector_size());
      for (int i = 0; i < n_frames; ++i) {
        d_kernel->generic_work(out, in, d_kernel->input_vector_size());
        out += d_kernel->output_vector_size();
        in += d_kernel->input_vector_size();
      }

      consume_each (n_frames * d_kernel->input_vector_size());

      // Tell runtime system how many output items we produced.
      return n_frames * d_kernel->output_vector_size();
    }

  } /* namespace gfdm */
} /* namespace gr */

