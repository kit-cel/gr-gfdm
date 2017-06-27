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

namespace gr {
  namespace gfdm {

    preamble_channel_estimator_cc::preamble_channel_estimator_cc(int timeslots, int fft_len, int active_subcarriers, std::vector<gfdm_complex> preamble):
      d_timeslots(timeslots), d_fft_len(fft_len), d_active_subcarriers(active_subcarriers)
    {
      d_preamble_fft_in = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_preamble_fft_out = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_preamble_fft_plan = initialize_fft(d_preamble_fft_out, d_preamble_fft_in, fft_len, true);

      d_inv_freq_preamble0 = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      d_inv_freq_preamble1 = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());
      initialize_inv_freq_preamble(d_inv_freq_preamble0, &preamble[0]);
      initialize_inv_freq_preamble(d_inv_freq_preamble1, &preamble[fft_len]);

      d_intermediate_channel_estimate = (gfdm_complex *) volk_malloc(sizeof(gfdm_complex) * fft_len, volk_get_alignment());

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

  } /* namespace gfdm */
} /* namespace gr */

