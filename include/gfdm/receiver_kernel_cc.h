/* -*- c++ -*- */
/*
 * Copyright 2016 Andrej Rode
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_GFDM_RECEIVER_KERNEL_CC_H
#define INCLUDED_GFDM_RECEIVER_KERNEL_CC_H

#include <complex>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <fftw3.h>
#include <stdexcept>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Demodulate a GFDM block
     *  This class initializes and performs all operations necessary to demodulate a GFDM block.
     *
     */
    /*!
     * \details
     * The GFDM receiver kernel class provides all necessary operations to blocks which instantiate it.
     * Further functions and methods not depending on GNU Radio should be implemented here.
     * This receiver implementation is based on [Gas+13].
     * It is recommended to use overlap = 2 to make use of sparse frequency domain processing.
     *
     * [Gas+13] I.S. Gaspar et al. "Low Complexity GFDM Receiver Based on Sparse Frequency Domain Processing"
     *
     */
    class  receiver_kernel_cc
    {
    public:
      typedef std::complex<float> gfdm_complex;
      typedef boost::shared_ptr<receiver_kernel_cc> sptr;

      receiver_kernel_cc(int n_timeslots, int n_subcarriers, int overlap, std::vector<gfdm_complex> frequency_taps);
      ~receiver_kernel_cc();

      void generic_work(gfdm_complex* out, const gfdm_complex* in);
      void generic_work_equalize(gfdm_complex *out, const gfdm_complex *in, const gfdm_complex* f_eq_in);
      void fft_filter_downsample(gfdm_complex* p_out, const gfdm_complex* p_in);
      void transform_subcarriers_to_td(gfdm_complex *p_out, const gfdm_complex *p_in);
      void cancel_sc_interference(gfdm_complex* p_out, const gfdm_complex* p_td_in, const gfdm_complex* p_fd_in);

      void filter_superposition(std::vector< std::vector<gfdm_complex> > &out, const gfdm_complex* in);
      void demodulate_subcarrier(std::vector< std::vector<gfdm_complex> > &out, std::vector< std::vector<gfdm_complex> > &sc_fdomain);
      void serialize_output(gfdm_complex out[], std::vector< std::vector<gfdm_complex> > &sc_symbols);
      void vectorize_2d(std::vector< std::vector<gfdm_complex> > &out_vector, const gfdm_complex* p_in);
      void remove_sc_interference(std::vector< std::vector<gfdm_complex> > &sc_symbols, std::vector< std::vector<gfdm_complex> > &sc_fdomain);
      int block_size()
      {
        return d_block_len;
      };
      std::vector<gfdm_complex> filter_taps();
      std::vector<gfdm_complex> ic_filter_taps();
      int timeslots() {return d_n_timeslots;};

    private:
      int d_n_subcarriers;
      int d_n_timeslots;
      int d_block_len;
      int d_overlap;
      gfdm_complex* d_filter_taps;
      gfdm_complex* d_ic_filter_taps;

      fftwf_plan initialize_fft(gfdm_complex* out_buf, gfdm_complex* in_buf, const int fft_size, bool forward);
      void initialize_taps_vector(gfdm_complex* filter_taps, std::vector<gfdm_complex> frequency_taps, const int n_timeslots);

      fftwf_plan d_in_fft_plan;
      gfdm_complex* d_in_fft_in;
      gfdm_complex* d_in_fft_out;

      gfdm_complex* d_equalized;

      fftwf_plan d_sc_ifft_plan;
      gfdm_complex* d_sc_ifft_in;
      gfdm_complex* d_sc_ifft_out;

      fftwf_plan d_sc_fft_plan;
      gfdm_complex* d_sc_fft_in;
      gfdm_complex* d_sc_fft_out;

      gfdm_complex* d_sc_postfilter;
      gfdm_complex* d_sc_filtered;

      void filter_subcarriers_and_downsample_fd(gfdm_complex *p_out, const gfdm_complex *p_in);

    };
  } /* namespace gfdm */
} /* namespace gr */

#endif /* INCLUDED_GFDM_RECEIVER_KERNEL_CC_H */




