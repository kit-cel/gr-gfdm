/* -*- c++ -*- */
/*
 * Copyright 2016 Johannes Demel.
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

#include "simple_modulator_cc_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace gfdm {

simple_modulator_cc::sptr
simple_modulator_cc::make(int n_timeslots,
                          int n_subcarriers,
                          int overlap,
                          std::vector<gr_complex> frequency_taps)
{
    return gnuradio::get_initial_sptr(new simple_modulator_cc_impl(
        n_timeslots, n_subcarriers, overlap, frequency_taps));
}

/*
 * The private constructor
 */
simple_modulator_cc_impl::simple_modulator_cc_impl(int n_timeslots,
                                                   int n_subcarriers,
                                                   int overlap,
                                                   std::vector<gr_complex> frequency_taps)
    : gr::sync_block("simple_modulator_cc",
                     gr::io_signature::make(1, 1, sizeof(gr_complex)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_kernel(std::unique_ptr<modulator_kernel_cc>(
          new modulator_kernel_cc(n_timeslots, n_subcarriers, overlap, frequency_taps)))
{
    // d_kernel = modulator_kernel_cc::sptr(new modulator_kernel_cc(n_timeslots,
    // n_subcarriers, overlap, frequency_taps));
    set_output_multiple(d_kernel->block_size());
}

/*
 * Our virtual destructor.
 */
simple_modulator_cc_impl::~simple_modulator_cc_impl() {}

int simple_modulator_cc_impl::work(int noutput_items,
                                   gr_vector_const_void_star& input_items,
                                   gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    const int n_blocks = noutput_items / d_kernel->block_size();
    //      std::cout << "noutput_items = " << noutput_items << ", block_size = " <<
    //      d_kernel->block_size() << ", #blocks = " << n_blocks << std::endl;
    for (int i = 0; i < n_blocks; ++i) {
        d_kernel->generic_work(out, in);
        in += d_kernel->block_size();
        out += d_kernel->block_size();
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace gfdm */
} /* namespace gr */
