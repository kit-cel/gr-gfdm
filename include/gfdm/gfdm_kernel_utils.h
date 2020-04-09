/* -*- c++ -*- */
/*
 * Copyright 2017 Johannes Demel.
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


#ifndef INCLUDED_GFDM_GFDM_KERNEL_UTILS_H
#define INCLUDED_GFDM_GFDM_KERNEL_UTILS_H

//#include <gfdm/api.h>
#include <fftw3.h>
#include <complex>
#include <stdexcept>
#include <vector>

namespace gr {
namespace gfdm {

/*!
 * \brief <+description+>
 *
 */
class gfdm_kernel_utils
{
public:
    typedef std::complex<float> gfdm_complex;

    gfdm_kernel_utils();
    ~gfdm_kernel_utils();

    fftwf_plan initialize_fft(gfdm_complex* out_buf,
                              gfdm_complex* in_buf,
                              const int fft_size,
                              bool forward);
    float calculate_signal_energy(const gfdm_complex* p_in, const int ninput_size);

private:
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_GFDM_KERNEL_UTILS_H */
