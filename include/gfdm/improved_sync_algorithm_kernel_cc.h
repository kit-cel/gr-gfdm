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
#include <deque>

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

      int detect_frame_start(const gr_complex *p_in, int ninput_size);

      // The following public functions are mainly a debugging interface to Python!
      int find_preamble(std::vector<gr_complex> in_vec);
      std::vector<gr_complex> auto_correlate_preamble(std::vector<gr_complex> in_vec);
      std::vector<float> abs_integrate_preamble(std::vector<gr_complex> in_vec);
      int find_peak_preamble(std::vector<float> in_vec);
      float calculate_normalized_cfo_preamble(const gr_complex val){ return calculate_normalized_cfo(val);};
      std::vector<gr_complex> remove_cfo_preamble(std::vector<gr_complex> in_vec, const float cfo);
      std::vector<gr_complex> cross_correlate_preamble(std::vector<gr_complex> in_vec);
      std::vector<gr_complex> preamble();
      void set_false_alarm_probability(float false_alarm_prob);
      std::vector<gr_complex> input_buffer();
      std::vector<gr_complex> auto_corr_buffer();
      std::vector<float> integration_buffer();
      std::vector<float> auto_corr_integrate(std::vector<gr_complex> in_vec);

    private:
      int d_buffer_len;
      float d_false_alarm_prob_factor;
      int d_n_subcarriers;
      int d_cp_len;
      gr_complex* d_preamble;
      gr_complex* d_p_in_buffer;
      gr_complex* d_auto_corr_vals;
      float* d_abs_auto_corr_vals;

      void auto_correlate(gr_complex* corr_vals, const gr_complex* p_in, const int ninput_size);

      // following functions take care of absolute value integration over CP length.
      std::deque<float> d_fifo;
      float integrate_fifo(float next_val);
      void abs_integrate(float* vals, const gr_complex* p_in, const int ninput_size);

      // calculate results from auto correlation.
      int find_peak(float* vals, const int ninput_size);

      // perform the auto correlation stage and write all results to the provided buffers!
      void perform_auto_correlation_stage(float *abs_corr_vals, gr_complex *corr_vals,
                                          const gr_complex *p_in, const int window_size);

      // derive subcarrier CFO from correlation value peak.
      float calculate_normalized_cfo(const gr_complex corr_val);

      // function assumes enough samples are available. Just lile xcorr stage does!
      void prepare_xcorr_input_array(gr_complex *xcorr_in, const gr_complex *p_in,
                                     const int offset);

      // following lines hold arrays and functions for xcorr fine STO peak detection.
      gr_complex* d_xcorr_vals;
      float* d_abs_xcorr_vals;
      int find_cross_correlation_peak(const gr_complex* p_in, const float* abs_int_vals, const float cfo);
      void remove_cfo(gr_complex* p_out, const gr_complex* p_in, const float cfo, const int ninput_size);
      void cross_correlate(gr_complex* p_out, const gr_complex* p_in, const int ninput_size);
      void combine_abs_auto_and_cross_correlation(float* p_out, const float* p_auto, const float* p_cross, const int ninput_size);


    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_IMPROVED_SYNC_ALGORITHM_KERNEL_CC_H */

