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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gfdm/improved_sync_algorithm_kernel_cc.h>
#include <stdexcept>
#include <iostream>
#include <string.h>
#include <complex>
#include <queue>

#include <volk/volk.h>

namespace gr {
  namespace gfdm {

    improved_sync_algorithm_kernel_cc::improved_sync_algorithm_kernel_cc(int n_subcarriers, int cp_len, std::vector<gr_complex> preamble):
      d_n_subcarriers(n_subcarriers), d_cp_len(cp_len)
    {
      // perform initial checks!
      if(preamble.size() != 2 * n_subcarriers){
        throw std::runtime_error("ERROR: preamble.size() MUST be equal to 2 * n_subcarriers!");
      }

      // calculate energy of preamble
      gr_complex energy = gr_complex(0.0, 0.0);
      volk_32fc_x2_conjugate_dot_prod_32fc(&energy, &preamble[0], &preamble[0], preamble.size());

      // now calculate amplitude, assume Q part == 0.0
      float amplitude = std::sqrt(energy.real());
      float scaling_factor = 1.0 / amplitude;

      std::cout << "preamble_energy: " << energy << ", amplitude: " << amplitude << std::endl;

      // malloc array for preamble and copy scaled version to array.
      d_preamble = (gr_complex*) volk_malloc(sizeof(gr_complex) * 2 * n_subcarriers, volk_get_alignment());
      volk_32f_s32f_multiply_32f((float*) d_preamble, (float*) &preamble[0], scaling_factor, 2 * 2 * n_subcarriers);

      // DEBUG ONLY: check if values match!
//      volk_32fc_x2_conjugate_dot_prod_32fc(&energy, d_preamble, d_preamble, 2 * n_subcarriers);
//      std::cout << "preamble_energy: " << energy << std::endl;
      d_fifo.resize(d_cp_len);
//      for (int i = 0; i < d_fifo.size(); ++i) {
//        std::cout << ", " << i << ":" << d_fifo.at(i);
//      }
//      std::cout << std::endl;

      d_abs_auto_corr_vals = (float*) volk_malloc(sizeof(float) * 2 * n_subcarriers, volk_get_alignment());

      d_xcorr_vals = (gr_complex*) volk_malloc(sizeof(gr_complex) * 4 * n_subcarriers, volk_get_alignment());
      d_abs_xcorr_vals = (float*) volk_malloc(sizeof(float) * 2 * n_subcarriers, volk_get_alignment());
    }

    improved_sync_algorithm_kernel_cc::~improved_sync_algorithm_kernel_cc()
    {
      volk_free(d_preamble);
      volk_free(d_abs_auto_corr_vals);
      volk_free(d_xcorr_vals);
      volk_free(d_abs_xcorr_vals);
    }

    std::vector<gr_complex> improved_sync_algorithm_kernel_cc::preamble()
    {
      std::vector<gr_complex> preamble(2 * d_n_subcarriers);
      memcpy(&preamble[0], d_preamble, sizeof(gr_complex) * 2 * d_n_subcarriers);
      return preamble;
    }

    int
    improved_sync_algorithm_kernel_cc::find_preamble(std::vector<gr_complex> in_vec)
    {
      int res = detect_frame_start(&in_vec[0], in_vec.size());
      return res;
    }

    int
    improved_sync_algorithm_kernel_cc::detect_frame_start(const gr_complex *p_in, int ninput_size)
    {
      const int buf_len = (ninput_size - 2 * d_n_subcarriers);

      gr_complex* corr_vals = (gr_complex*) volk_malloc(sizeof(gr_complex) * buf_len, volk_get_alignment());
      float* abs_int_vals = (float*) volk_malloc(sizeof(float) * buf_len, volk_get_alignment());

      const float threshold = 0.5;
      const int step_size = 1070;
      const int max_part_size = step_size + 2 * d_n_subcarriers;
      for (int i = 0; i < buf_len; i += step_size) {
        const int part_step_size = std::min(step_size, buf_len - i);
        auto_correlate(corr_vals + i, p_in + i, part_step_size + 2 * d_n_subcarriers);
        abs_integrate(abs_int_vals + i, corr_vals + i, part_step_size);
        int nm = find_peak(abs_int_vals + i, std::min(step_size, buf_len - i));
        if (*(abs_int_vals + i + nm) > threshold && std::min(nm, step_size - nm) < d_n_subcarriers){
          std::cout << "threshold detected! " << nm << std::endl;
          i -= d_n_subcarriers;
        }
      }
//      auto_correlate(corr_vals, p_in, ninput_size);
//      abs_integrate(abs_int_vals, corr_vals, buf_len);

      int nm = find_peak(abs_int_vals, buf_len);
      float cfo = calculate_normalized_cfo(corr_vals[nm]);

      const int offset = nm - d_n_subcarriers;
      int rnc = offset + find_cross_correlation_peak(p_in + offset, abs_int_vals + offset, cfo);

      volk_free(corr_vals);
      volk_free(abs_int_vals);

      return rnc;
    }

    std::vector<gr_complex>
    improved_sync_algorithm_kernel_cc::auto_correlate_preamble(std::vector<gr_complex> in_vec)
    {
      const int p_len = 2 * d_n_subcarriers;
      const int buf_len = in_vec.size() - p_len;
      std::vector<gr_complex> res(buf_len, 0);
      auto_correlate(&res[0], &in_vec[0], in_vec.size());
      return res;
    }

    void
    improved_sync_algorithm_kernel_cc::auto_correlate(gr_complex* corr_vals, const gr_complex* p_in, const int ninput_size)
    {
      const int p_len = 2 * d_n_subcarriers;
      const int buf_len = ninput_size - p_len;
      gr_complex energy = gr_complex(0.0, 0.0);
      gr_complex val = gr_complex(0.0, 0.0);
      for (int i = 0; i < buf_len; ++i) {
        // calculate symbol energy for normalization
        volk_32fc_x2_conjugate_dot_prod_32fc(&energy, p_in, p_in, p_len);

        // correlate over half preamble length
        // ATTENTION: second array is conjugated! Not first!
        volk_32fc_x2_conjugate_dot_prod_32fc(&val, p_in + p_len / 2, p_in, p_len / 2);

        // normalize result!
        float abs_energy = 0.5 * energy.real();
        *corr_vals++ = val / abs_energy;

        ++p_in;
      }
    }

    std::vector<float>
    improved_sync_algorithm_kernel_cc::abs_integrate_preamble(std::vector<gr_complex> in_vec)
    {
      std::vector<float> res(in_vec.size(), 0);
      abs_integrate(&res[0], &in_vec[0], in_vec.size());
      return res;
    }

    float improved_sync_algorithm_kernel_cc::integrate_fifo(float next_val)
    {
      d_fifo.push_back(next_val);

      float fifo_sum = 0.0;
      for (int j = 0; j < d_fifo.size(); ++j) {
        fifo_sum += d_fifo.at(j);
      }
      d_fifo.pop_front();
      return fifo_sum;
    }

    void
    improved_sync_algorithm_kernel_cc::abs_integrate(float* vals, const gr_complex* p_in, const int ninput_size)
    {
      const float norm_factor = 1.0 / (d_cp_len + 1.0);
      for (int i = 0; i < ninput_size; ++i) {
        float fifo_sum = integrate_fifo(std::abs(*p_in++));
        *vals++ = fifo_sum * norm_factor;
      }
    }

    int improved_sync_algorithm_kernel_cc::find_peak_preamble(std::vector<float> in_vec)
    {
      return find_peak(&in_vec[0], in_vec.size());
    }

    int
    improved_sync_algorithm_kernel_cc::find_peak(float* vals, const int ninput_size)
    {
      unsigned int nm = 0;
      volk_32f_index_max_32u(&nm, vals, ninput_size);
//      unsigned int nm = std::distance(vals, std::max_element(vals, vals + ninput_size));
      return (int) nm;
    }

    float improved_sync_algorithm_kernel_cc::calculate_normalized_cfo(const gr_complex corr_val)
    {
      return std::arg(corr_val) / M_PI;
    }

    int
    improved_sync_algorithm_kernel_cc::find_cross_correlation_peak(const gr_complex* p_in, const float* abs_int_vals, const float cfo)
    {
      const int n_corr_vals = 2 * d_n_subcarriers;
      const int buf_len = 2 * n_corr_vals;

      remove_cfo(d_xcorr_vals, p_in, cfo, buf_len);
      cross_correlate(d_xcorr_vals, d_xcorr_vals, buf_len);
      volk_32fc_magnitude_32f(d_abs_xcorr_vals, d_xcorr_vals, n_corr_vals);
      combine_abs_auto_and_cross_correlation(d_abs_xcorr_vals, abs_int_vals, d_abs_xcorr_vals, n_corr_vals);

      return find_peak(d_abs_xcorr_vals, n_corr_vals);
    }

    std::vector<gr_complex>
    improved_sync_algorithm_kernel_cc::remove_cfo_preamble(std::vector<gr_complex> in_vec, const float cfo)
    {
      std::vector<gr_complex> res(in_vec.size(), 0);
      remove_cfo(&res[0], &in_vec[0], cfo, in_vec.size());
      return res;
    }

    void
    improved_sync_algorithm_kernel_cc::remove_cfo(gr_complex* p_out, const gr_complex* p_in, const float cfo, const int ninput_size)
    {
      gr_complex initial_phase = gr_complex(1.0, 0.0);
      const float cfo_incr = -1.0f * M_PI * cfo / d_n_subcarriers;
      gr_complex phase_increment = gr_complex(std::cos(cfo_incr), std::sin(cfo_incr));
      volk_32fc_s32fc_x2_rotator_32fc(p_out, p_in, phase_increment, &initial_phase, ninput_size);
    }

    std::vector<gr_complex> improved_sync_algorithm_kernel_cc::cross_correlate_preamble(std::vector<gr_complex> in_vec)
    {
      std::vector<gr_complex> res(in_vec.size() - 2 * d_n_subcarriers, 0);
      cross_correlate(&res[0], &in_vec[0], in_vec.size());
      return res;
    }


    void
    improved_sync_algorithm_kernel_cc::cross_correlate(gr_complex* p_out, const gr_complex* p_in, const int ninput_size)
    {
      const int p_len = 2 * d_n_subcarriers;
      const int buf_len = ninput_size - p_len;
      for(int i = 0; i < buf_len; ++i){
        volk_32fc_x2_conjugate_dot_prod_32fc(p_out++, p_in++, d_preamble, p_len);
      }
    }

    void
    improved_sync_algorithm_kernel_cc::combine_abs_auto_and_cross_correlation(float* p_out, const float* p_auto, const float* p_cross, const int ninput_size)
    {
      volk_32f_x2_multiply_32f(p_out, p_auto, p_cross, ninput_size);
    }

  } /* namespace gfdm */
} /* namespace gr */

