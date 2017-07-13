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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gfdm/preamble_channel_estimator_cc.h>
#include <volk/volk.h>
#include <cstring>
#include <cmath>
#include <iostream>

namespace gr {
  namespace gfdm {

    preamble_channel_estimator_cc::preamble_channel_estimator_cc(int timeslots, int fft_len, int active_subcarriers, bool is_dc_free, std::vector<gfdm_complex> preamble):
      d_timeslots(timeslots), d_fft_len(fft_len), d_active_subcarriers(active_subcarriers), d_is_dc_free(is_dc_free), d_n_gaussian_taps(9)
    {
      d_preamble_fft_in = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_preamble_fft_out = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_preamble_fft_plan = initialize_fft(d_preamble_fft_out, d_preamble_fft_in, fft_len, true);

      d_inv_freq_preamble0 = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_inv_freq_preamble1 = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      initialize_inv_freq_preamble(d_inv_freq_preamble0, &preamble[0]);
      initialize_inv_freq_preamble(d_inv_freq_preamble1, &preamble[fft_len]);

      d_intermediate_channel_estimate = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());

      d_gaussian_taps = (float *) volk_malloc(sizeof(float) * d_n_gaussian_taps, volk_get_alignment());
      initialize_gaussian_filter(d_gaussian_taps, 1.0f, d_n_gaussian_taps);

      int intermediate_filter_len = active_subcarriers + d_n_gaussian_taps + (is_dc_free ? 1 : 0);
      d_filter_intermediate = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * intermediate_filter_len, volk_get_alignment());

      d_preamble_estimate = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_filtered_estimate = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * active_subcarriers + (is_dc_free ? 1 : 0), volk_get_alignment());

//      self._x_preamble = x_preamble
//      self._fft_len = fft_len
//      self._timeslots = timeslots
//      self._inv_freq_x_preamble0 = 1. / np.fft.fft(x_preamble[0:fft_len])
//      self._inv_freq_x_preamble1 = 1. / np.fft.fft(x_preamble[fft_len:])
//
//      active_sc = np.arange((self._fft_len - active_subcarriers)//2, (self._fft_len + active_subcarriers)//2+1)
//      active_sc = active_sc[3:-3]
//      freqs = np.fft.fftfreq(fft_len)
//      freqs = np.fft.fftshift(freqs)
//      self._active_preamble_freqs = freqs[active_sc]
//      self._active_sc = active_sc
//      fr_freqs = np.fft.fftfreq(self._fft_len * self._timeslots)
//      self._frame_freqs = np.fft.fftshift(fr_freqs)
//
//      g = signal.gaussian(9, 1.0)
//      g_factor = 1.
//      g /= np.sqrt(g_factor * g.dot(g))
//      self._p_filter = g
    }

    preamble_channel_estimator_cc::~preamble_channel_estimator_cc()
    {
    }

    void preamble_channel_estimator_cc::initialize_gaussian_filter(float* taps, const float sigma_sq, const int n_taps)
    {
      float s = 0.0f;
      for(int i = 0; i < n_taps; ++i){
        float val = std::pow(float(i - (n_taps / 2)), 2.0f) / sigma_sq;
        taps[i] = std::exp(-0.5f * val);
        s += taps[i];
      }

      for(int i = 0; i < n_taps; ++i){
        taps[i] = taps[i] / s;
      }
    }

    std::vector<float>
    preamble_channel_estimator_cc::preamble_filter_taps()
    {
      std::vector<float> t = std::vector<float>(d_n_gaussian_taps);
      for(int i = 0; i < d_n_gaussian_taps; ++i){
        t[i] = d_gaussian_taps[i];
      }
      return t;
    }

    void
    preamble_channel_estimator_cc::initialize_inv_freq_preamble(gfdm_complex* p_out, const gfdm_complex* p_preamble_part)
    {
      memcpy(d_preamble_fft_in, p_preamble_part, sizeof(gfdm_complex) * d_fft_len);
      fftwf_execute(d_preamble_fft_plan);
      for(int i = 0; i < d_fft_len; ++i){
        p_out[i] = gfdm_complex(0.5, 0.0) / d_preamble_fft_out[i];
      }
    }

    void
    preamble_channel_estimator_cc::estimate_preamble_channel(gfdm_complex* fd_preamble_channel, const gfdm_complex* rx_preamble)
    {
      estimate_fftlen_preamble_channel(d_intermediate_channel_estimate, rx_preamble, d_inv_freq_preamble0);
      estimate_fftlen_preamble_channel(fd_preamble_channel, rx_preamble + d_fft_len, d_inv_freq_preamble1);
      volk_32f_x2_add_32f((float*) fd_preamble_channel, (float*) fd_preamble_channel, (float*) d_intermediate_channel_estimate, 2 * d_fft_len);
    }

    void
    preamble_channel_estimator_cc::estimate_fftlen_preamble_channel(gfdm_complex* p_out, const gfdm_complex* rx_samples, const gfdm_complex* fd_ref_samples)
    {
      memcpy(d_preamble_fft_in, rx_samples, sizeof(gfdm_complex) * d_fft_len);
      fftwf_execute(d_preamble_fft_plan);
      volk_32fc_x2_multiply_32fc(p_out, d_preamble_fft_out, fd_ref_samples, d_fft_len);
    }

    void
    preamble_channel_estimator_cc::filter_preamble_estimate(gfdm_complex* filtered, const gfdm_complex* estimate)
    {
      for(int i = 0; i < d_n_gaussian_taps / 2; ++i){
        d_filter_intermediate[i] = estimate[d_fft_len - d_active_subcarriers / 2];
      }

      for(int i = 0; i < d_active_subcarriers / 2; ++i){
        d_filter_intermediate[i + d_n_gaussian_taps / 2] = estimate[i + d_fft_len - d_active_subcarriers / 2];
      }

      int offset = 0;
      if(d_is_dc_free){
        d_filter_intermediate[d_n_gaussian_taps / 2 + d_active_subcarriers / 2] = (estimate[d_fft_len - 1] + estimate[1]) / gfdm_complex(2.0f, 0.0f);
        offset = 1;
      }

      for(int i = 0; i < d_active_subcarriers / 2; ++i){
        d_filter_intermediate[i + offset + d_n_gaussian_taps / 2 + d_active_subcarriers / 2] = estimate[offset + i];
      }

      for(int i = d_active_subcarriers / 2; i < d_active_subcarriers / 2 + d_n_gaussian_taps / 2; ++i){
        d_filter_intermediate[i + offset + d_n_gaussian_taps / 2 + d_active_subcarriers / 2] = estimate[offset + d_active_subcarriers / 2 - 1];
      }

      int n_taps = d_active_subcarriers + offset;
      for(int i = 0; i < n_taps; ++i){
        volk_32fc_32f_dot_prod_32fc(filtered + i, d_filter_intermediate + i, d_gaussian_taps, d_n_gaussian_taps);
      }

    }

    void
    preamble_channel_estimator_cc::interpolate_frame(gfdm_complex* frame_estimate, const gfdm_complex* estimate)
    {
      const int n_estimated_taps = d_active_subcarriers + (d_is_dc_free ? 1 : 0);
      const int center = d_fft_len * d_timeslots / 2;
      const int dead_subcarriers = d_fft_len - d_active_subcarriers;

      gfdm_complex step_size = gfdm_complex(1.0f / float(d_timeslots), 0.0f);

      for(int i = center; i < center + d_timeslots * dead_subcarriers / 2; ++i) {
        frame_estimate[i] = estimate[0];
      }

      for(int i = d_timeslots * d_active_subcarriers / 2; i < center; ++i) {
        frame_estimate[i] = estimate[n_estimated_taps - 1];
      }

      for(int i = 0; i < n_estimated_taps / 2; ++i) {
        gfdm_complex inc = (estimate[i + 1] - estimate[i]) * step_size;
        gfdm_complex factor = estimate[i];
        for(int j = 0; j < d_timeslots; ++j) {
          frame_estimate[center + d_timeslots * dead_subcarriers / 2 + i * d_timeslots + j] = factor;
          factor += inc;
        }
      }

      for(int i = n_estimated_taps / 2; i < n_estimated_taps - 1; ++i) {
        int offset = (i - n_estimated_taps / 2) * d_timeslots;
        gfdm_complex inc = (estimate[i + 1] - estimate[i]) * step_size;
        gfdm_complex factor = estimate[i];
        for(int j = 0; j < d_timeslots; ++j) {
          frame_estimate[offset + j] = factor;
          factor += inc;
        }
      }
    }

    void
    preamble_channel_estimator_cc::estimate_frame(gfdm_complex* frame_estimate, const gfdm_complex* rx_preamble)
    {
      estimate_preamble_channel(d_preamble_estimate, rx_preamble);
      filter_preamble_estimate(d_filtered_estimate, d_preamble_estimate);
      interpolate_frame(frame_estimate, d_filtered_estimate);
    }

  } /* namespace gfdm */
} /* namespace gr */

