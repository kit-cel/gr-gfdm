/* -*- c++ -*- */
/*
 * Copyright 2016 Andrej Rode, Johannes Demel.
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

#include "remove_prefix_cc_impl.h"
#include <gnuradio/io_signature.h>
#include <pmt/pmt.h>

namespace gr {
namespace gfdm {

remove_prefix_cc::sptr remove_prefix_cc::make(int frame_len,
                                              int block_len,
                                              int offset,
                                              const std::string& gfdm_sync_tag_key)
{
    return gnuradio::make_block_sptr<remove_prefix_cc_impl>(
        frame_len, block_len, offset, gfdm_sync_tag_key);
}

/*
 * The private constructor
 */
remove_prefix_cc_impl::remove_prefix_cc_impl(int frame_len,
                                             int block_len,
                                             int offset,
                                             const std::string& gfdm_sync_tag_key)
    : gr::block("remove_prefix_cc",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_frame_len(frame_len),
      d_block_len(block_len),
      d_offset(offset),
      d_block_left(0)
{
    d_tag_key = pmt::string_to_symbol(gfdm_sync_tag_key);
    set_output_multiple(block_len);
    set_fixed_rate(true);
    set_relative_rate(1. * block_len / frame_len);
    set_tag_propagation_policy(TPP_DONT);
}

/*
 * Our virtual destructor.
 */
remove_prefix_cc_impl::~remove_prefix_cc_impl() {}

void remove_prefix_cc_impl::forecast(int noutput_items,
                                     gr_vector_int& ninput_items_required)
{
    ninput_items_required[0] = fixed_rate_noutput_to_ninput(noutput_items);
}

int remove_prefix_cc_impl::fixed_rate_ninput_to_noutput(int ninput)
{
    return (ninput / d_frame_len) * d_block_len;
}

int remove_prefix_cc_impl::fixed_rate_noutput_to_ninput(int noutput)
{
    return (noutput / d_block_len) * d_frame_len;
}

int remove_prefix_cc_impl::general_work(int noutput_items,
                                        gr_vector_int& ninput_items,
                                        gr_vector_const_void_star& input_items,
                                        gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];
    // This general_work relies on the scheduler in several ways!
    // due to fixed rate and correct forecast return values, ninput_items is a multiple of
    // d_frame_len also, noutput_items is a multiple of d_block_len. Just crash if that
    // doesn't hold true. (hopefully)

    int frames = ninput_items[0] / d_frame_len;
    int blocks = noutput_items / d_block_len;
    int avail_items = std::min(blocks, frames) * d_frame_len;
    int consumed_items = 0;
    int produced_items = 0;

    std::vector<tag_t> tags;
    get_tags_in_window(tags, 0, 0, avail_items, d_tag_key);
    for (int i = 0; i < tags.size(); ++i) {
        memcpy(out, in + d_offset, sizeof(gr_complex) * d_block_len);
        add_item_tag(0, nitems_written(0) + produced_items, d_tag_key, tags[i].value);
        consumed_items += d_frame_len;
        produced_items += d_block_len;
        in += d_frame_len;
        out += d_block_len;
    }

    consume_each(consumed_items);
    return produced_items;
}

} /* namespace gfdm */
} /* namespace gr */
