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

#include "advanced_receiver_sb_cc_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>

namespace gr {
namespace gfdm {

advanced_receiver_sb_cc::sptr
advanced_receiver_sb_cc::make(int n_timeslots,
                              int n_subcarriers,
                              int overlap,
                              int ic_iter,
                              std::vector<gr_complex> frequency_taps,
                              gr::digital::constellation_sptr constellation,
                              std::vector<int> subcarrier_map,
                              int do_phase_compensation)
{
    return gnuradio::get_initial_sptr(
        new advanced_receiver_sb_cc_impl(n_timeslots,
                                         n_subcarriers,
                                         overlap,
                                         ic_iter,
                                         frequency_taps,
                                         constellation,
                                         subcarrier_map,
                                         do_phase_compensation));
}

/*
 * The private constructor
 */
advanced_receiver_sb_cc_impl::advanced_receiver_sb_cc_impl(
    int n_timeslots,
    int n_subcarriers,
    int overlap,
    int ic_iter,
    std::vector<gr_complex> frequency_taps,
    gr::digital::constellation_sptr constellation,
    std::vector<int> subcarrier_map,
    int do_phase_compensation)
    : gr::sync_block("advanced_receiver_sb_cc",
                     gr::io_signature::make(1, 2, sizeof(gr_complex)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex)))
{
    // int do_phase_compensation = 1;
    d_adv_kernel = std::unique_ptr<advanced_receiver_kernel_cc>(
        new advanced_receiver_kernel_cc(n_timeslots,
                                        n_subcarriers,
                                        overlap,
                                        frequency_taps,
                                        subcarrier_map,
                                        ic_iter,
                                        constellation,
                                        do_phase_compensation));
    set_output_multiple(d_adv_kernel->block_size());
    set_tag_propagation_policy(TPP_DONT);
}

/*
 * Our virtual destructor.
 */
advanced_receiver_sb_cc_impl::~advanced_receiver_sb_cc_impl() {}

int advanced_receiver_sb_cc_impl::work(int noutput_items,
                                       gr_vector_const_void_star& input_items,
                                       gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    const int n_blocks = noutput_items / d_adv_kernel->block_size();

    std::vector<tag_t> tags;
    if (input_items.size() > 1) {
        const gr_complex* in_eq = (const gr_complex*)input_items[1];
        for (int i = 0; i < n_blocks; ++i) {
            d_adv_kernel->generic_work_equalize(out, in, in_eq);

            in += d_adv_kernel->block_size();
            in_eq += d_adv_kernel->block_size();
            out += d_adv_kernel->block_size();
        }

        get_tags_in_window(tags, 1, 0, n_blocks * d_adv_kernel->block_size());
    } else {
        for (int i = 0; i < n_blocks; ++i) {
            d_adv_kernel->generic_work(out, in);

            in += d_adv_kernel->block_size();
            out += d_adv_kernel->block_size();
        }

        get_tags_in_window(tags, 0, 0, n_blocks * d_adv_kernel->block_size());
    }

    for (auto t : tags) {
        add_item_tag(0, t);
    }

    return n_blocks * d_adv_kernel->block_size();
}

} /* namespace gfdm */
} /* namespace gr */
