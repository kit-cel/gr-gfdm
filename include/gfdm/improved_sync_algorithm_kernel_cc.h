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


#ifndef INCLUDED_GFDM_IMPROVED_SYNC_ALGORITHM_KERNEL_CC_H
#define INCLUDED_GFDM_IMPROVED_SYNC_ALGORITHM_KERNEL_CC_H

#include <gfdm/api.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Perform STO/CFO synchronization for multicarrier systems
     * Compare: Awoseyila et. al. "Improved Preamble-Aided Timing Estimation for OFDM Systems"
     *
     * \param n_subcarriers: number of subcarriers of the multicarrier system.
     * \param cp_len: Cyclic Prefix Length of the system
     * \param preamble: actual used preamble. preamble.size MUST be equal to 2*n_subcarriers.
     *
     */
    class GFDM_API improved_sync_algorithm_kernel_cc
    {
    public:
      improved_sync_algorithm_kernel_cc(int n_subcarriers, int cp_len, std::vector<gr_complex> preamble);
      ~improved_sync_algorithm_kernel_cc();

      int generic_work(gr_complex* p_out, const gr_complex* p_in, int ninput_size);
      std::vector<gr_complex> find_preamble(std::vector<gr_complex> in_vec);
    private:
      struct auto_correlation_result_t{
        int nm;
        float cfo;
      };
      int d_n_subcarriers;
      int d_cp_len;
      gr_complex* d_preamble;
      float* d_abs_auto_corr_vals;

      float calculate_normalized_cfo(const gr_complex corr_val);

      auto_correlation_result_t find_auto_correlation_peak(float* abs_auto_corr_vals, const gr_complex* p_in, int ninput_size);
      void remove_cfo(gr_complex* p_out, const gr_complex* p_in, const float cfo, const float ninput_size);

      float* d_cc_float;
      gr_complex* d_cc_complex;
      int find_exact_cross_correlation_peak(const gr_complex* p_in, const float* abs_auto_corr_vals);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_IMPROVED_SYNC_ALGORITHM_KERNEL_CC_H */

