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

    transmitter_kernel::transmitter_kernel(int timeslots, int subcarriers,
                                           int active_subcarriers, int cp_len,
                                           int cs_len, int ramp_len,
                                           std::vector<int> subcarrier_map,
                                           bool per_timeslot, int overlap,
                                           std::vector<gfdm_complex> frequency_taps,
                                           std::vector<gfdm_complex> window_taps,
                                           std::vector<gfdm_complex> preamble)
        : d_preamble(preamble),
          d_mapper(std::unique_ptr<resource_mapper_kernel_cc>(new resource_mapper_kernel_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, per_timeslot))),
          d_modulator(std::unique_ptr<modulator_kernel_cc>(new modulator_kernel_cc(timeslots, subcarriers, overlap, frequency_taps))),
          d_prefixer(std::unique_ptr<add_cyclic_prefix_cc>(new add_cyclic_prefix_cc(timeslots * subcarriers, cp_len, cs_len, ramp_len, window_taps)))
    {
        // d_mapper = resource_mapper_kernel_cc::sptr(
        //     new resource_mapper_kernel_cc(timeslots, subcarriers, active_subcarriers,
        //                                   subcarrier_map, per_timeslot));

        d_mapped = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * d_mapper->output_vector_size(),
                                                volk_get_alignment());
        d_frame = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * d_modulator->block_size(),
                                               volk_get_alignment());

    }

    transmitter_kernel::~transmitter_kernel()
    {
    }

    void
    transmitter_kernel::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in,
                                     const int ninput_size)
    {
        d_mapper->map_to_resources(d_mapped, p_in, ninput_size);
        d_modulator->generic_work(d_frame, d_mapped);
        d_prefixer->generic_work(p_out + d_preamble.size(), d_frame);
        std::memcpy(p_out, d_preamble.data(), sizeof(gfdm_complex) * d_preamble.size());
    }


  } /* namespace gfdm */
} /* namespace gr */

