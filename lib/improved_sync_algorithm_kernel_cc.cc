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
#include <deque>
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
      volk_32fc_x2_conjugate_dot_prod_32fc(&energy, d_preamble, d_preamble, 2 * n_subcarriers);
      std::cout << "preamble_energy: " << energy << std::endl;

      d_abs_auto_corr_vals = (float*) volk_malloc(sizeof(float) * 2 * n_subcarriers, volk_get_alignment());
      d_cc_float = (float*) volk_malloc(sizeof(float) * 2 * n_subcarriers, volk_get_alignment());
      d_cc_complex = (gr_complex*) volk_malloc(sizeof(gr_complex) * 2 * n_subcarriers, volk_get_alignment());
    }

    improved_sync_algorithm_kernel_cc::~improved_sync_algorithm_kernel_cc()
    {
      volk_free(d_preamble);
      volk_free(d_abs_auto_corr_vals);
    }

    std::vector<gr_complex> improved_sync_algorithm_kernel_cc::find_preamble(std::vector<gr_complex> in_vec)
    {
      std::vector<gr_complex> res(in_vec.size() - 2 * d_n_subcarriers);
      generic_work(&res[0], &in_vec[0], in_vec.size());
      return res;
    }

    int
    improved_sync_algorithm_kernel_cc::generic_work(gr_complex* p_out, const gr_complex* p_in, int ninput_size)
    {
      const int buf_len = (ninput_size - 2 * d_n_subcarriers);
      float* abs_corr_vals = (float*) volk_malloc(sizeof(float) * buf_len, volk_get_alignment());
      auto_correlation_result_t res = find_auto_correlation_peak(abs_corr_vals, p_in, ninput_size);
      std::cout << "max element position: " << res.nm << ", corrval: " << abs_corr_vals[res.nm] << ", CFO:" << res.cfo << std::endl;
      memcpy(p_out, abs_corr_vals, sizeof(float) * buf_len);
      remove_cfo(p_out, p_in, res.cfo, buf_len);
      const int offset = res.nm - d_n_subcarriers;
      int nc = find_exact_cross_correlation_peak(p_out + offset, abs_corr_vals + offset);
      std::cout << "exact_preamble_start: " << nc << ", result: " << res.nm - d_n_subcarriers + nc << std::endl;

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

    void
    improved_sync_algorithm_kernel_cc::abs_integrate(float* vals, const gr_complex* p_in, const int ninput_size)
    {
      std::deque<float> fifo(d_cp_len, 0.0);
      const float norm_factor = 1.0 / (d_cp_len + 1.0);

      for (int i = 0; i < ninput_size; ++i) {
        fifo.push_back(std::abs(*p_in++));

        float fifo_sum = 0.0;
        for (int j = 0; j < fifo.size(); ++j) {
          fifo_sum += fifo.at(j);
        }
        fifo.pop_front();

        *vals++ = fifo_sum * norm_factor;
      }
    }

    void
    improved_sync_algorithm_kernel_cc::auto_correlate_integrate(float* abs_corr_vals, gr_complex* corr_vals, const gr_complex* p_in, const int ninput_size)
    {
      const int p_len = 2 * d_n_subcarriers;
      const int buf_len = ninput_size - p_len;
      std::deque<float> fifo_abs_corr_vals(d_cp_len, 0.0);
      gr_complex energy = gr_complex(0.0, 0.0);
      for (int i = 0; i < buf_len; ++i) {
        // calculate symbol energy for normalization
        volk_32fc_x2_conjugate_dot_prod_32fc(&energy, p_in, p_in, p_len);
        float abs_energy = energy.real();

        // correlate over half preamble length
        // ATTENTION: second array is conjugated! Not first!
        volk_32fc_x2_conjugate_dot_prod_32fc(corr_vals + i, p_in + p_len / 2, p_in, p_len / 2);

        float cval = 2.0 * std::abs(corr_vals[i]) / abs_energy;
        fifo_abs_corr_vals.push_back(cval);
        float fifo_sum = 0.0;
        for (int j = 0; j < fifo_abs_corr_vals.size(); ++j) {
          fifo_sum += fifo_abs_corr_vals.at(j);
        }
        fifo_abs_corr_vals.pop_front();
        abs_corr_vals[i] = fifo_sum / (d_cp_len + 1.0);

        ++p_in;
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

    improved_sync_algorithm_kernel_cc::auto_correlation_result_t
    improved_sync_algorithm_kernel_cc::find_auto_correlation_peak(float* abs_auto_corr_vals, const gr_complex* p_in, int ninput_size)
    {
      const int p_len = 2 * d_n_subcarriers;
      const int buf_len = ninput_size - p_len;

      gr_complex* corr_vals = (gr_complex*) volk_malloc(sizeof(gr_complex) * buf_len, volk_get_alignment());
      auto_correlate(corr_vals, p_in, ninput_size);
      float* abs_corr_vals = (float*) volk_malloc(sizeof(float) * buf_len, volk_get_alignment());
      abs_integrate(abs_corr_vals, corr_vals, buf_len);

      auto_correlation_result_t res;
      res.nm = find_peak(abs_corr_vals, buf_len);
      res.cfo = calculate_normalized_cfo(corr_vals[res.nm]);

      memcpy(abs_auto_corr_vals, abs_corr_vals, sizeof(float) * buf_len);
      return res;
    }

    float improved_sync_algorithm_kernel_cc::calculate_normalized_cfo(const gr_complex corr_val)
    {
      return std::arg(corr_val) / M_PI;
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


    void improved_sync_algorithm_kernel_cc::cross_correlate(gr_complex* p_out, const gr_complex* p_in, const int ninput_size)
    {
      const int p_len = 2 * d_n_subcarriers;
      const int buf_len = ninput_size - p_len;
      for(int i = 0; i < buf_len; ++i){
//        volk_32fc_x2_dot_prod_32fc(p_out++, p_in++, d_preamble, p_len);
//        volk_32fc_x2_conjugate_dot_prod_32fc(p_out++, d_preamble, p_in++, p_len);
        volk_32fc_x2_conjugate_dot_prod_32fc(p_out++, p_in++, d_preamble, p_len);
      }
    }

    int
    improved_sync_algorithm_kernel_cc::find_exact_cross_correlation_peak(const gr_complex *p_in, const float *abs_auto_corr_vals)
    {
      const int p_len = 2 * d_n_subcarriers;
      for(int i = 0; i < p_len; ++i){
        volk_32fc_x2_dot_prod_32fc(d_cc_complex + i, p_in + i, d_preamble, p_len);
      }
      volk_32fc_magnitude_squared_32f(d_cc_float, d_cc_complex, p_len);
      volk_32f_x2_multiply_32f(d_cc_float, d_cc_float, abs_auto_corr_vals, p_len);
      unsigned int nc = 0;
      volk_32f_index_max_32u(&nc, d_cc_float, p_len);
      return nc;
    }

  } /* namespace gfdm */
} /* namespace gr */

