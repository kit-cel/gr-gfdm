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


#ifndef INCLUDED_GFDM_MODULATOR_KERNEL_CC_H
#define INCLUDED_GFDM_MODULATOR_KERNEL_CC_H

#include "gfdm_kernel_utils.h"
#include <fftw3.h>
#include <volk/volk_alloc.hh>
#include <complex>
#include <stdexcept>
#include <vector>

namespace gr {
namespace gfdm {

void foobar();
/*!
 * \brief modulate a GFDM block.
 *  This class initializes and performs all operations necessary to modulate a GFDM block.
 *
 */
class modulator_kernel_cc : public gfdm_kernel_utils
{
public:
    modulator_kernel_cc(int n_timeslots,
                        int n_subcarriers,
                        int overlap,
                        std::vector<gfdm_complex> frequency_taps);
    ~modulator_kernel_cc();
    void generic_work(gfdm_complex* p_out, const gfdm_complex* p_in);
    int block_size() { return d_n_subcarriers * d_n_timeslots; };
    std::vector<gfdm_complex> filter_taps();

private:
    int d_n_timeslots;
    int d_n_subcarriers;
    int d_ifft_len;
    int d_overlap;
    volk::vector<gfdm_complex> d_filter_taps;

    void initialize_taps_vector(gfdm_complex* filter_taps,
                                std::vector<gfdm_complex> frequency_taps,
                                const int n_timeslots);

    volk::vector<gfdm_complex> d_sub_fft_in;
    volk::vector<gfdm_complex> d_sub_fft_out;
    fftwf_plan d_sub_fft_plan;
    volk::vector<gfdm_complex> d_filtered;
    volk::vector<gfdm_complex> d_ifft_in;
    volk::vector<gfdm_complex> d_ifft_out;
    fftwf_plan d_ifft_plan;

    // DEBUG function
    const void print_vector(const gfdm_complex* v, const int size);
    static bool complex_compare(gfdm_complex i, gfdm_complex j)
    {
        return std::abs(i) < std::abs(j);
    };
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_MODULATOR_KERNEL_CC_H */
