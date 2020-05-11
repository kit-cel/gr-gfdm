/* -*- c++ -*- */
/*
 * Copyright 2018 Johannes Demel.
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

#include "transmitter_cc_impl.h"
#include <gnuradio/io_signature.h>
#include <chrono>

namespace gr {
namespace gfdm {

transmitter_cc::sptr transmitter_cc::make(int timeslots,
                                          int subcarriers,
                                          int active_subcarriers,
                                          int cp_len,
                                          int cs_len,
                                          int ramp_len,
                                          std::vector<int> subcarrier_map,
                                          bool per_timeslot,
                                          int overlap,
                                          std::vector<gr_complex> frequency_taps,
                                          std::vector<gr_complex> window_taps,
                                          std::vector<gr_complex> preamble,
                                          const std::string& tsb_tag_key)
{
    return gnuradio::get_initial_sptr(new transmitter_cc_impl(timeslots,
                                                              subcarriers,
                                                              active_subcarriers,
                                                              cp_len,
                                                              cs_len,
                                                              ramp_len,
                                                              subcarrier_map,
                                                              per_timeslot,
                                                              overlap,
                                                              frequency_taps,
                                                              window_taps,
                                                              preamble,
                                                              tsb_tag_key));
}

/*
 * The private constructor
 */
transmitter_cc_impl::transmitter_cc_impl(int timeslots,
                                         int subcarriers,
                                         int active_subcarriers,
                                         int cp_len,
                                         int cs_len,
                                         int ramp_len,
                                         std::vector<int> subcarrier_map,
                                         bool per_timeslot,
                                         int overlap,
                                         std::vector<gr_complex> frequency_taps,
                                         std::vector<gr_complex> window_taps,
                                         std::vector<gr_complex> preamble,
                                         const std::string& tsb_tag_key)
    : gr::block("transmitter_cc",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_length_tag_key_str(tsb_tag_key),
      d_length_tag_key(pmt::string_to_symbol(tsb_tag_key))
{
    d_kernel =
        std::unique_ptr<transmitter_kernel>(new transmitter_kernel(timeslots,
                                                                   subcarriers,
                                                                   active_subcarriers,
                                                                   cp_len,
                                                                   cs_len,
                                                                   ramp_len,
                                                                   subcarrier_map,
                                                                   per_timeslot,
                                                                   overlap,
                                                                   frequency_taps,
                                                                   window_taps,
                                                                   preamble));
    set_relative_rate(1.0 * d_kernel->output_vector_size() /
                      d_kernel->input_vector_size());
    set_fixed_rate(true);
    set_output_multiple(d_kernel->output_vector_size());
}

/*
 * Our virtual destructor.
 */
transmitter_cc_impl::~transmitter_cc_impl() {}

void transmitter_cc_impl::forecast(int noutput_items,
                                   gr_vector_int& ninput_items_required)
{
    for (int i = 0; i < ninput_items_required.size(); ++i) {
        ninput_items_required[i] = fixed_rate_noutput_to_ninput(noutput_items);
    }
}

int transmitter_cc_impl::fixed_rate_ninput_to_noutput(int ninput)
{
    return (ninput / d_kernel->input_vector_size()) * d_kernel->output_vector_size();
}

int transmitter_cc_impl::fixed_rate_noutput_to_ninput(int noutput)
{
    return (noutput / d_kernel->output_vector_size()) * d_kernel->input_vector_size();
}

int transmitter_cc_impl::general_work(int noutput_items,
                                      gr_vector_int& ninput_items,
                                      gr_vector_const_void_star& input_items,
                                      gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    // std::vector<tag_t> tags;
    // get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + noutput_items,
    // pmt::string_to_symbol("time")); for(auto t: tags){
    //   auto cn = std::chrono::high_resolution_clock::now().time_since_epoch();
    //   auto s = std::chrono::nanoseconds(pmt::to_long(t.value));
    //   auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(cn - s);
    //   std::cout << "Transmitter before duration: " << d.count() << "ns" << std::endl;
    // }

    int n_frames = std::min(noutput_items / d_kernel->output_vector_size(),
                            ninput_items[0] / d_kernel->input_vector_size());

    if (not d_length_tag_key_str.empty()) {
        std::vector<tag_t> tags;
        get_tags_in_range(tags,
                          0,
                          nitems_read(0),
                          nitems_read(0) + n_frames * d_kernel->input_vector_size(),
                          d_length_tag_key);
        for (auto tag : tags) {
            if (tag.key == d_length_tag_key) {
                //GR_LOG_INFO(this->d_logger, "length: " + std::to_string(pmt::to_long(tag.value)) +
                //        "\tkey: " + pmt::symbol_to_string(tag.key));	
    //                             std::to_string(header_duration.count()) +
    //                             "ns");
                assert(pmt::to_long(tag.value) == d_kernel->input_vector_size());
                remove_item_tag(0, tag);
            }
        }
    }


    for (int i = 0; i < n_frames; ++i) {
        d_kernel->generic_work(out, in, d_kernel->input_vector_size());
        out += d_kernel->output_vector_size();
        in += d_kernel->input_vector_size();
    }

    consume_each(n_frames * d_kernel->input_vector_size());

    if (not d_length_tag_key_str.empty()) {
        for (int i = 0; i < n_frames; ++i) {
            add_item_tag(0,
                         nitems_written(0) + i * d_kernel->output_vector_size(),
                         d_length_tag_key,
                         pmt::from_long(d_kernel->output_vector_size()));
        }
    }
    // for(auto t: tags){
    //   auto cn = std::chrono::high_resolution_clock::now().time_since_epoch();
    //   auto s = std::chrono::nanoseconds(pmt::to_long(t.value));
    //   auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(cn - s);
    //   std::cout << "Transmitter after duration: " << d.count() << "ns" << std::endl;
    // }

    // Tell runtime system how many output items we produced.
    return n_frames * d_kernel->output_vector_size();
}

} /* namespace gfdm */
} /* namespace gr */
