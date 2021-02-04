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
                                       std::vector<gfdm_complex> preamble)
    : d_preamble(preamble.begin(), preamble.end()),
      d_mapper(std::make_unique<resource_mapper_kernel_cc>(
          timeslots, subcarriers, active_subcarriers, subcarrier_map, per_timeslot)),
      d_modulator(std::make_unique<modulator_kernel_cc>(
          timeslots, subcarriers, overlap, frequency_taps)),
      d_prefixer(std::make_unique<add_cyclic_prefix_cc>(
          timeslots * subcarriers, cp_len, cs_len, ramp_len, window_taps))
{

    d_mapped.resize(d_mapper->output_vector_size());
    d_frame.resize(d_modulator->block_size());
}

transmitter_kernel::~transmitter_kernel() {}

void transmitter_kernel::generic_work(gfdm_complex* p_out,
                                      const gfdm_complex* p_in,
                                      const int ninput_size)
{
    d_mapper->map_to_resources(d_mapped.data(), p_in, ninput_size);
    d_modulator->generic_work(d_frame.data(), d_mapped.data());
    d_prefixer->generic_work(p_out + d_preamble.size(), d_frame.data());
    std::memcpy(p_out, d_preamble.data(), sizeof(gfdm_complex) * d_preamble.size());
}


} /* namespace gfdm */
} /* namespace gr */
