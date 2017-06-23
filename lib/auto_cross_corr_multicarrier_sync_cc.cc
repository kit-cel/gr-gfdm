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


#include <gfdm/auto_cross_corr_multicarrier_sync_cc.h>
#include <volk/volk.h>
#include <string.h>  // memset requires this include?!
#include <iostream>

namespace gr {
  namespace gfdm {

    auto_cross_corr_multicarrier_sync_cc::auto_cross_corr_multicarrier_sync_cc(int subcarriers, int cp_len, std::vector<gfdm_complex> preamble):
      d_subcarriers(subcarriers), d_cp_len(cp_len), d_last_cfo(0.0f), d_preamble_attenuation(1.0f)
    {
      if(int(preamble.size()) != 2 * subcarriers){
        throw std::runtime_error("ERROR: preamble.size() MUST be equal to 2 * n_subcarriers!");
      }

      // calculate energy of preamble
//      gfdm_complex energy = gfdm_complex(0.0, 0.0);
//      volk_32fc_x2_conjugate_dot_prod_32fc(&energy, &preamble[0], &preamble[0], preamble.size());
      float energy = calculate_signal_energy(&preamble[0], preamble.size());
      d_reference_preamble_energy = energy;
      // now calculate amplitude, assume Q part == 0.0
      float amplitude = std::sqrt(energy / 2.0);
      float scaling_factor = 1.0 / amplitude;

      // malloc array for preamble and copy scaled version to array.
      d_preamble = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * 2 * subcarriers, volk_get_alignment());
      volk_32f_s32f_multiply_32f((float*) d_preamble, (float*) &preamble[0], scaling_factor, 2 * 2 * subcarriers);

        std::cout << "energy: " << energy << ", my_e: " << calculate_signal_energy(d_preamble, 2 * subcarriers) << ", amplitude: " << amplitude << std::endl;

      d_buffer_len = 4 * subcarriers;
      d_auto_corr = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
      memset(d_auto_corr, 0, sizeof(gfdm_complex) * d_buffer_len);
      d_abs_auto_corr  = (float*) volk_malloc(sizeof(float) * d_buffer_len, volk_get_alignment());

      d_xcorr = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * 2 * subcarriers, volk_get_alignment());
      d_abs_xcorr = (float*) volk_malloc(sizeof(float) * 2 * subcarriers, volk_get_alignment());

      // This part is necessary for heavy sync optimization. e.g. ~120 -> ~45 us latency
      d_fxc_in = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
      d_fxc_out = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
      d_ixc_in = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
      d_ixc_out = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
      d_freq_preamble = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
      d_fxc_plan = initialize_fft(d_fxc_out, d_fxc_in, d_buffer_len, true);
      d_ixc_plan = initialize_fft(d_ixc_out, d_ixc_in, d_buffer_len, false);
      memset(d_fxc_in, 0x0, sizeof(gfdm_complex) * d_buffer_len);
      memcpy(d_fxc_in, d_preamble, sizeof(gfdm_complex) * 2 * subcarriers);
      fftwf_execute(d_fxc_plan);
      memcpy(d_freq_preamble, d_fxc_out, sizeof(gfdm_complex) * d_buffer_len);
    }

    auto_cross_corr_multicarrier_sync_cc::~auto_cross_corr_multicarrier_sync_cc()
    {
      volk_free(d_preamble);
      volk_free(d_auto_corr);
      volk_free(d_abs_auto_corr);
      volk_free(d_xcorr);
      volk_free(d_abs_xcorr);

      fftwf_destroy_plan(d_fxc_plan);
      fftwf_destroy_plan(d_ixc_plan);
      volk_free(d_fxc_in);
      volk_free(d_fxc_out);
      volk_free(d_ixc_in);
      volk_free(d_ixc_out);
      volk_free(d_freq_preamble);

    }

      float
      auto_cross_corr_multicarrier_sync_cc::calculate_signal_energy(const gfdm_complex* p_in, const int ninput_size)
      {
          gfdm_complex energy = gfdm_complex(0.0, 0.0);
          volk_32fc_x2_conjugate_dot_prod_32fc(&energy, p_in, p_in, ninput_size);
          return energy.real();
      }

    fftwf_plan
    auto_cross_corr_multicarrier_sync_cc::initialize_fft(gfdm_complex *out_buf, gfdm_complex *in_buf, const int fft_size, bool forward)
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
    auto_cross_corr_multicarrier_sync_cc::calculate_preamble_attenuation(const gfdm_complex* p_in)
    {
      return std::sqrt(d_reference_preamble_energy / calculate_signal_energy(p_in, 2 * d_subcarriers));
    }

    void
    auto_cross_corr_multicarrier_sync_cc::normalize_power_level(gfdm_complex* p_out, const gfdm_complex* p_in, const float norm_factor, const int ninput_size)
    {
//      volk_32f_s32f_normalize((float*) p_out, const float scalar, unsigned int num_points)
      volk_32f_s32f_multiply_32f((float*) p_out, (const float*) p_in, norm_factor, 2 * ninput_size);
    }


    int
    auto_cross_corr_multicarrier_sync_cc::detect_frame_start(const gfdm_complex *p_in, int ninput_size)
    {
      adjust_buffer_size(ninput_size); // make sure we don't try to write to some random memory.
      const int p_len = 2 * d_subcarriers;
      const int buf_len = ninput_size - p_len;

      fixed_lag_auto_correlate(d_auto_corr, p_in, ninput_size);

      volk_32fc_magnitude_squared_32f(d_abs_auto_corr, d_auto_corr, buf_len);
      const int nm = find_peak(d_abs_auto_corr, buf_len);
      d_last_cfo = calculate_normalized_cfo(d_auto_corr[nm]);

      // don forget! : ac = np.roll(oac, cp_len // 2) -> aka move left
      const int pos_correction_factor = d_cp_len / 2;
      const int xc_start = nm + pos_correction_factor - d_subcarriers;

      cross_correlate_preamble(d_xcorr, p_in + xc_start, 4 * d_subcarriers);
      volk_32fc_magnitude_squared_32f(d_abs_xcorr, d_xcorr, 2 * d_subcarriers);

      volk_32f_x2_multiply_32f(d_abs_xcorr, d_abs_xcorr, d_abs_auto_corr + xc_start, 2 * d_subcarriers);
      const int p_nc = find_peak(d_abs_xcorr, 2 * d_subcarriers);

      const int nc = xc_start + p_nc;
      d_preamble_attenuation = calculate_preamble_attenuation(p_in + nc);
//      std::cout << "preamble attenuation: " << d_preamble_attenuation << ", " << d_abs_xcorr[p_nc]  << ", " << d_abs_auto_corr[nc] << std::endl;
      d_frame_phase = std::arg(d_xcorr[p_nc]);
//      std::cout << "nc: " << nc << "(" << nm + pos_correction_factor << "), cfo: " << d_last_cfo << std::endl;
      return nc;
    }

    void
    auto_cross_corr_multicarrier_sync_cc::fixed_lag_auto_correlate(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size)
    {
      const int p_len = 2 * d_subcarriers;
      const int buf_len = ninput_size - p_len;
      gfdm_complex val = gfdm_complex(0.0, 0.0);
      for (int i = 0; i < buf_len; ++i) {
        // correlate over half preamble length
        // ATTENTION: second array is conjugated! Not first!
        volk_32fc_x2_conjugate_dot_prod_32fc(&val, p_in + p_len / 2, p_in, p_len / 2);
        *p_out++ = val;
        ++p_in;
      }
    }


    int
    auto_cross_corr_multicarrier_sync_cc::find_peak(float* vals, const int ninput_size)
    {
      unsigned int nm = 0;
      volk_32f_index_max_32u(&nm, vals, ninput_size);
      return (int) nm;
    }

    float
    auto_cross_corr_multicarrier_sync_cc::calculate_normalized_cfo(const gfdm_complex corr_val)
    {
      return std::arg(corr_val) / (2.0 * M_PI);
    }

    void
    auto_cross_corr_multicarrier_sync_cc::cross_correlate_preamble(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size)
    {
      const int p_len = 2 * d_subcarriers;
//      const int buf_len = ninput_size - p_len;
//      for(int i = 0; i < buf_len; ++i){
//        volk_32fc_x2_conjugate_dot_prod_32fc(p_out++, p_in++, d_preamble, p_len);
////        std::cout << "res: " << *(p_out - 1) << ",\tin: " << *(p_in - 1) << ",\tp: " << d_preamble[i] << std::endl;
//      }
//      std::cout << "buffer_len: " << d_buffer_len << ", ninput_size: " << ninput_size << std::endl;
      const int fft_len = 4 * d_subcarriers;
      memcpy(d_fxc_in, p_in, sizeof(gfdm_complex) * fft_len);
      fftwf_execute(d_fxc_plan);
      volk_32fc_x2_multiply_conjugate_32fc(d_ixc_in, d_fxc_out, d_freq_preamble, fft_len);
      fftwf_execute(d_ixc_plan);
      memcpy(p_out, d_ixc_out, sizeof(gfdm_complex) * p_len);
    }

    void
    auto_cross_corr_multicarrier_sync_cc::adjust_buffer_size(const int ninput_size)
    {
      if(ninput_size > d_buffer_len){
        volk_free(d_auto_corr);
        volk_free(d_abs_auto_corr);
        d_buffer_len = ninput_size;
        d_auto_corr = (gfdm_complex*) volk_malloc(sizeof(gfdm_complex) * d_buffer_len, volk_get_alignment());
        d_abs_auto_corr  = (float*) volk_malloc(sizeof(float) * d_buffer_len, volk_get_alignment());
      }
    }

  } /* namespace gfdm */
} /* namespace gr */

