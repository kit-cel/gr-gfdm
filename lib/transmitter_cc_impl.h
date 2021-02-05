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

#ifndef INCLUDED_GFDM_TRANSMITTER_CC_IMPL_H
#define INCLUDED_GFDM_TRANSMITTER_CC_IMPL_H

#include <gfdm/transmitter_cc.h>
#include <gfdm/transmitter_kernel.h>
#include <volk/volk_alloc.hh>

namespace gr {
namespace gfdm {

class transmitter_cc_impl : public transmitter_cc
{
private:
    std::unique_ptr<transmitter_kernel> d_kernel;
    std::string d_length_tag_key_str;
    pmt::pmt_t d_length_tag_key;
    volk::vector<gr_complex> d_modulated;

public:
    transmitter_cc_impl(int timeslots,
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
                        std::vector<int> cyclic_shifts,
                        std::vector<std::vector<gr_complex>> preambles,
                        const std::string& tsb_tag_key = "");
    ~transmitter_cc_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);
    int fixed_rate_ninput_to_noutput(int ninput);
    int fixed_rate_noutput_to_ninput(int noutput);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_TRANSMITTER_CC_IMPL_H */
