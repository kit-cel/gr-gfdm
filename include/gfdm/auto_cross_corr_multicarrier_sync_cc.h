/* -*- c++ -*- */
/* 
 * Copyright 2016 <+YOU OR YOUR COMPANY+>.
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


#ifndef INCLUDED_GFDM_AUTO_CROSS_CORR_MULTICARRIER_SYNC_CC_H
#define INCLUDED_GFDM_AUTO_CROSS_CORR_MULTICARRIER_SYNC_CC_H

#include <complex>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Simplified version of "Improved Preamble-Aided Timing Estimation for OFDM Systems"
     *
     */
    class auto_cross_corr_multicarrier_sync_cc
    {
    public:
      typedef std::complex<float> gfdm_complex;
      typedef boost::shared_ptr<auto_cross_corr_multicarrier_sync_cc> sptr;

      auto_cross_corr_multicarrier_sync_cc(int subcarriers, int cp_len, std::vector<gfdm_complex> preamble);
      ~auto_cross_corr_multicarrier_sync_cc();

      int detect_frame_start(const gfdm_complex *p_in, int ninput_size);
    private:
      int d_subcarriers;
      int d_cp_len;
      int d_buffer_len;

      gfdm_complex* d_preamble;
      gfdm_complex* d_auto_corr;
      float* d_abs_auto_corr;
      gfdm_complex* d_xcorr;
      float* d_abs_xcorr;

      int find_peak(float* vals, const int ninput_size);
      float calculate_normalized_cfo(const gfdm_complex corr_val);
      void cross_correlate_preamble(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size);
      void adjust_buffer_size(const int ninput_size);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_AUTO_CROSS_CORR_MULTICARRIER_SYNC_CC_H */

