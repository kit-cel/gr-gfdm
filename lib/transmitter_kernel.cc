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

#include <gnuradio/io_signature.h>
#include <gfdm/transmitter_kernel.h>
#include <volk/volk.h>
#include <cstring>
#include <exception>

namespace gr {
namespace gfdm {

transmitter_kernel::transmitter_kernel(int timeslots,
                                       int subcarriers,
                                       int active_subcarriers,
                                       int cp_len,
                                       int cs_len,
                                       int ramp_len,
                                       std::vector<int> subcarrier_map,
                                       bool per_timeslot,
                                       int overlap,
                                       std::vector<gfdm_complex> frequency_taps,
                                       std::vector<gfdm_complex> window_taps,
                                       std::vector<int> cyclic_shifts,
                                       std::vector<std::vector<gfdm_complex>> preambles)
    : d_mapper(std::make_unique<resource_mapper_kernel_cc>(
          timeslots, subcarriers, active_subcarriers, subcarrier_map, per_timeslot)),
      d_modulator(std::make_unique<modulator_kernel_cc>(
          timeslots, subcarriers, overlap, frequency_taps)),
      d_prefixer(std::make_unique<add_cyclic_prefix_cc>(
          timeslots * subcarriers, cp_len, cs_len, ramp_len, window_taps)),
      d_cyclic_shifts(cyclic_shifts.begin(), cyclic_shifts.end()),
      d_preamble_size(preambles[0].size())
{
    if (cyclic_shifts.size() != preambles.size()) {
        throw std::invalid_argument(
            "Number of cyclic shifts and number of preambles do not match!");
    }

    for (const auto& p : preambles) {
        if (d_preamble_size != p.size()) {
            throw std::invalid_argument("All preambles must have equal size!");
        }
    }

    for (unsigned i = 0; i < cyclic_shifts.size(); ++i) {
        d_preambles.emplace(
            cyclic_shifts[i],
            volk::vector<gfdm_complex>(preambles[i].begin(), preambles[i].end()));
    }
    d_mapped.resize(d_mapper->output_vector_size());
    d_frame.resize(d_modulator->block_size());
}

transmitter_kernel::~transmitter_kernel() {}

void transmitter_kernel::modulate(gfdm_complex* out,
                                  const gfdm_complex* in,
                                  const int ninput_size)
{
    d_mapper->map_to_resources(d_mapped.data(), in, ninput_size);
    d_modulator->generic_work(out, d_mapped.data());
}

void transmitter_kernel::insert_preamble(gfdm_complex* out, const int cyclic_shift)
{
    std::memcpy(
        out, d_preambles[cyclic_shift].data(), sizeof(gfdm_complex) * d_preamble_size);
}

void transmitter_kernel::add_frame(gfdm_complex* out,
                                   const gfdm_complex* in,
                                   const int cyclic_shift)
{
    insert_preamble(out, cyclic_shift);
    d_prefixer->add_cyclic_prefix(out + d_preamble_size, in, cyclic_shift);
}


void transmitter_kernel::generic_work(gfdm_complex* p_out,
                                      const gfdm_complex* p_in,
                                      const int ninput_size)
{
    modulate(d_frame.data(), p_in, ninput_size);
    add_frame(p_out, d_frame.data(), d_cyclic_shifts[0]);
}


} /* namespace gfdm */
} /* namespace gr */
