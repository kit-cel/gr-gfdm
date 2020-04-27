/* -*- c++ -*- */
/*
 * Copyright 2017 <+YOU OR YOUR COMPANY+>.
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

#include "channel_estimator_cc_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace gfdm {

channel_estimator_cc::sptr channel_estimator_cc::make(int timeslots,
                                                      int fft_len,
                                                      int active_subcarriers,
                                                      bool is_dc_free,
                                                      int which_estimator,
                                                      std::vector<gr_complex> preamble)
{
    return gnuradio::get_initial_sptr(new channel_estimator_cc_impl(
        timeslots, fft_len, active_subcarriers, is_dc_free, which_estimator, preamble));
}

/*
 * The private constructor
 */
channel_estimator_cc_impl::channel_estimator_cc_impl(int timeslots,
                                                     int fft_len,
                                                     int active_subcarriers,
                                                     bool is_dc_free,
                                                     int which_estimator,
                                                     std::vector<gr_complex> preamble)
    : gr::block("channel_estimator_cc",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex)))
{
    d_estimator_kernel = std::unique_ptr<preamble_channel_estimator_cc>(
        new preamble_channel_estimator_cc(timeslots,
                                          fft_len,
                                          active_subcarriers,
                                          is_dc_free,
                                          which_estimator,
                                          preamble));

    // set block properties!
    set_relative_rate(timeslots / 2.0);
    set_fixed_rate(true);
    set_output_multiple(fft_len * timeslots);
}

/*
 * Our virtual destructor.
 */
channel_estimator_cc_impl::~channel_estimator_cc_impl() {}

void channel_estimator_cc_impl::forecast(int noutput_items,
                                         gr_vector_int& ninput_items_required)
{
    for (int i = 0; i < ninput_items_required.size(); ++i) {
        ninput_items_required[i] = fixed_rate_noutput_to_ninput(noutput_items);
    }
}

int channel_estimator_cc_impl::fixed_rate_ninput_to_noutput(int ninput)
{
    return ninput * d_estimator_kernel->timeslots() / 2;
}

int channel_estimator_cc_impl::fixed_rate_noutput_to_ninput(int noutput)
{
    return 2 * noutput / d_estimator_kernel->timeslots();
}

int channel_estimator_cc_impl::general_work(int noutput_items,
                                            gr_vector_int& ninput_items,
                                            gr_vector_const_void_star& input_items,
                                            gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    const int invec_len = 2 * d_estimator_kernel->fft_len();
    const int frame_len = d_estimator_kernel->frame_len();
    const int n_frames = noutput_items / frame_len;

    for (int i = 0; i < n_frames; ++i) {
        d_estimator_kernel->estimate_frame(out, in);
        /*
        const float snr_lin = d_estimator_kernel->estimate_snr(in);
        add_item_tag(0,
                     nitems_written(0) + i * frame_len,
                     pmt::intern("snr_lin"),
                     pmt::from_float(snr_lin));
        */
        in += invec_len;
        out += frame_len;
    }

    consume_each(n_frames * invec_len);
    return n_frames * frame_len;
}

} /* namespace gfdm */
} /* namespace gr */
