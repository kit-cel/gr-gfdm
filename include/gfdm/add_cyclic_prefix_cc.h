/* -*- c++ -*- */
/*
 * Copyright 2016, 2021 Johannes Demel.
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


#ifndef INCLUDED_GFDM_ADD_CYCLIC_PREFIX_CC_H
#define INCLUDED_GFDM_ADD_CYCLIC_PREFIX_CC_H

#include <volk/volk_alloc.hh>
#include <complex>
#include <stdexcept>
#include <vector>

namespace gr {
namespace gfdm {

/*!
 * \brief Kernel adds cyclic prefix to GFDM frame and applies block pinching window.
 *
 */

class add_cyclic_prefix_cc
{
public:
    typedef std::complex<float> gfdm_complex;

    add_cyclic_prefix_cc(int block_len,
                         int cp_len,
                         int cs_len,
                         int ramp_len,
                         std::vector<gfdm_complex> window_taps,
                         int cyclic_shift = 0);
    ~add_cyclic_prefix_cc();
    void generic_work(gfdm_complex* p_out, const gfdm_complex* p_in);
    void add_cyclic_prefix(gfdm_complex* p_out,
                           const gfdm_complex* p_in,
                           const int cyclic_prefix);
    void remove_cyclic_prefix(gfdm_complex* p_out, const gfdm_complex* p_in);
    int block_size() { return d_block_len; };
    int frame_size() { return block_size() + d_cp_len + d_cs_len; };
    int cyclic_shift() const { return d_cyclic_shift; };

private:
    const int d_block_len;
    const int d_cp_len;
    const int d_cs_len;
    const int d_ramp_len;
    const int d_cyclic_shift;

    void add_cyclic_extension(gfdm_complex* out,
                              const gfdm_complex* in,
                              const int cyclic_shift);

    volk::vector<gfdm_complex> d_front_ramp;
    volk::vector<gfdm_complex> d_back_ramp;
    void apply_ramp(gfdm_complex* out, const gfdm_complex* in);
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_ADD_CYCLIC_PREFIX_CC_H */
