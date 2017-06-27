/* -*- c++ -*- */
/* 
 * Copyright 2017 <+YOU OR YOUR COMPANY+>.
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

//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif
#include <volk/volk.h>
#include <gnuradio/io_signature.h>
#include <gfdm/gfdm_kernel_utils.h>

namespace gr {
  namespace gfdm {

    gfdm_kernel_utils::gfdm_kernel_utils()
    {
    }

    gfdm_kernel_utils::~gfdm_kernel_utils()
    {
    }

    fftwf_plan
    gfdm_kernel_utils::initialize_fft(gfdm_complex *out_buf, gfdm_complex *in_buf, const int fft_size, bool forward)
    {
      std::string filename(getenv("HOME"));
      filename += "/.gr_fftw_wisdom";
      FILE *fpr = fopen(filename.c_str(), "r");
      if (fpr != 0) {
        fftwf_import_wisdom_from_file(fpr);
        fclose(fpr);
      }

      fftwf_plan plan = fftwf_plan_dft_1d(fft_size,
                                          reinterpret_cast<fftwf_complex *>(in_buf),
                                          reinterpret_cast<fftwf_complex *>(out_buf),
                                          forward ? FFTW_FORWARD : FFTW_BACKWARD,
                                          FFTW_MEASURE);

      FILE *fpw = fopen(filename.c_str(), "w");
      if (fpw != 0) {
        fftwf_export_wisdom_to_file(fpw);
        fclose(fpw);
      }
      return plan;
    }

    float
    gfdm_kernel_utils::calculate_signal_energy(const gfdm_complex* p_in, const int ninput_size)
    {
      gfdm_complex energy = gfdm_complex(0.0, 0.0);
      volk_32fc_x2_conjugate_dot_prod_32fc(&energy, p_in, p_in, ninput_size);
      return energy.real();
    }

  } /* namespace gfdm */
} /* namespace gr */

