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

// These imports are GR related. Removed them to be able to compile a stand-alone kernel.
//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif
//#include <gnuradio/io_signature.h>

// These are the important includes
#include <gfdm/modulator_kernel_cc.h>
#include <iostream>
#include <volk/volk.h>
#include <fftw3.h>
#include <string.h>
#include <algorithm>


namespace gr {
  namespace gfdm {

    modulator_kernel_cc::modulator_kernel_cc(int n_timeslots, int n_subcarriers, int overlap, std::vector<gfdm_complex> frequency_taps):
      d_n_timeslots(n_timeslots), d_n_subcarriers(n_subcarriers), d_overlap(overlap)
    {
      if (int(frequency_taps.size()) != n_timeslots * overlap){
        throw std::invalid_argument("number of frequency taps MUST be equal to n_timeslots * overlap!");
      }
      d_filter_taps = (gfdm_complex *) volk_malloc (sizeof (gfdm_complex) * n_timeslots * overlap, volk_get_alignment ());
      memcpy(d_filter_taps, &frequency_taps[0], sizeof(gfdm_complex) * n_timeslots * overlap);

      // first create input and output buffers for a new FFTW plan.
      d_sub_fft_in = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots, volk_get_alignment ());
      d_sub_fft_out = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots, volk_get_alignment ());
      d_sub_fft_plan = fftwf_plan_dft_1d (n_timeslots,
                                  reinterpret_cast<fftwf_complex *>(d_sub_fft_in),
                                  reinterpret_cast<fftwf_complex *>(d_sub_fft_out),
                                  FFTW_FORWARD,
                                  FFTW_MEASURE);
      d_fd_data = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots * n_subcarriers, volk_get_alignment ());
      d_upfilter = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots * overlap, volk_get_alignment ());
      d_filtered = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_subcarriers * n_timeslots * overlap, volk_get_alignment ());
      d_fd_out = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * (n_subcarriers * n_timeslots + (overlap - 1) * n_timeslots), volk_get_alignment ());

      d_ifft_in = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots * n_subcarriers, volk_get_alignment ());
      d_ifft_out = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots * n_subcarriers, volk_get_alignment ());
      d_ifft_plan = fftwf_plan_dft_1d (n_timeslots * n_subcarriers,
                                       reinterpret_cast<fftwf_complex *>(d_ifft_in),
                                       reinterpret_cast<fftwf_complex *>(d_ifft_out),
                                       FFTW_BACKWARD,
                                       FFTW_MEASURE);
      std::cout << "CTOR kernel finished!\n";
    }

    modulator_kernel_cc::~modulator_kernel_cc()
    {
      volk_free(d_filter_taps);
      fftwf_destroy_plan ((fftwf_plan) d_sub_fft_plan);
      volk_free(d_sub_fft_in);
      volk_free(d_sub_fft_out);

      volk_free(d_fd_data);
      volk_free(d_upfilter);
      volk_free(d_filtered);
      volk_free(d_fd_out);
      fftwf_destroy_plan((fftwf_plan) d_ifft_plan);
      volk_free(d_ifft_in);
      volk_free(d_ifft_out);
    }

    void
    modulator_kernel_cc::gfdm_fftshift(gfdm_complex* p_out, const gfdm_complex* p_in, const int size)
    {
      const unsigned int len = (unsigned int)(ceil(size/2.0));
      memcpy(p_out, p_in + len, sizeof(gfdm_complex) * (size - len));
      memcpy(p_out + (size - len), p_in, sizeof(gfdm_complex) * len);
    }

    void
    modulator_kernel_cc::subcarrier_fft(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      memcpy(d_sub_fft_in, p_in, sizeof(gfdm_complex) * d_n_timeslots);
      fftwf_execute((fftwf_plan) d_sub_fft_plan);
      memcpy(p_out, d_sub_fft_out, sizeof(gfdm_complex) * d_n_timeslots);
    }

    void
    modulator_kernel_cc::block_subcarrier_fft(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      for(int k = 0; k < d_n_subcarriers; ++k){
        subcarrier_fft(p_out, p_in);
        p_in += d_n_timeslots;
        p_out += d_n_timeslots;
      }
    }

    void
    modulator_kernel_cc::upsample_and_filter(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
//      gfdm_complex* upfilter = d_upfilter;
//      for(int p = 0; p < d_overlap; ++p){
//        memcpy(upfilter, p_in, sizeof(gfdm_complex) * d_n_timeslots);
//        upfilter += d_n_timeslots;
//      }
//      volk_32fc_x2_multiply_32fc(p_out, d_upfilter, d_filter_taps, d_n_timeslots * d_overlap);

      gfdm_complex* taps = d_filter_taps;
      for(int p = 0; p < d_overlap; ++p){
        volk_32fc_x2_multiply_32fc(p_out, p_in, taps, d_n_timeslots);
        p_out += d_n_timeslots;
        taps += d_n_timeslots;
      }
    }

    void
    modulator_kernel_cc::block_upsample_and_filter(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      for(int k = 0; k < d_n_subcarriers; ++k){
        upsample_and_filter(p_out, p_in);
        p_in += d_n_timeslots;
        p_out += d_n_timeslots * d_overlap;
      }
    }

    const void
    modulator_kernel_cc::print_vector(const gfdm_complex* v, const int size)
    {
      for (int i = 0; i < size; ++i) {
        std::cout << v[i] << ", ";
        if(i % 8 == 0 && i > 0){
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    void
    modulator_kernel_cc::combine_subcarriers(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      const int tail_length = (d_overlap - 1) * d_n_timeslots;
      const int fft_len = d_n_timeslots * d_overlap;
      const int fft_half = fft_len / 2;
      memset(d_fd_out, 0x00, sizeof (gfdm_complex) * (d_n_subcarriers * d_n_timeslots + (d_overlap - 1) * d_n_timeslots));

      gfdm_complex* temp = d_fd_out;
      for(int k = 0; k < d_n_subcarriers; ++k){
        volk_32f_x2_add_32f((float*) temp, (float*) temp, (float*) (p_in + fft_half), fft_len);
        volk_32f_x2_add_32f((float*) (temp + fft_half), (float*) (temp + fft_half), (float*) p_in, fft_len);

        p_in += d_n_timeslots * d_overlap;
        temp += d_n_timeslots;
      }
      volk_32f_x2_add_32f((float*) d_fd_out, (float*) d_fd_out, (float*) (d_fd_out + d_n_timeslots * d_n_subcarriers), 2 * tail_length);
//      std::rotate_copy(d_fd_out, d_fd_out + (d_n_timeslots * d_overlap / 2), d_fd_out + d_n_subcarriers * d_n_timeslots, p_out);
      std::rotate_copy(d_fd_out, d_fd_out + fft_half, d_fd_out + d_n_subcarriers * d_n_timeslots, p_out);
    }

    void
    modulator_kernel_cc::modulate_block(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      const int part_len = std::min(d_n_timeslots * d_overlap / 2, d_n_timeslots);
      memset(p_out, 0x00, sizeof (gfdm_complex) * d_n_subcarriers * d_n_timeslots);

      for(int k = 0; k < d_n_subcarriers; ++k){
        memcpy(d_sub_fft_in, p_in, sizeof(gfdm_complex) * d_n_timeslots);
        fftwf_execute((fftwf_plan) d_sub_fft_plan);
        upsample_and_filter(d_filtered, d_sub_fft_out);

        for (int i = 0; i < d_overlap; ++i) {
          int src_part_pos = ((i + d_overlap / 2) % d_overlap) * d_n_timeslots;
          int target_part_pos = ((k + i + d_n_subcarriers - (d_overlap / 2)) % d_n_subcarriers) * d_n_timeslots;
          volk_32f_x2_add_32f((float*) (p_out + target_part_pos), (float*) (p_out + target_part_pos), (float*) (d_filtered + src_part_pos), 2 * part_len);
        }

        p_in += d_n_timeslots;
      }
    }

    void
    modulator_kernel_cc::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      // assume data symbols are stacked subcarrier-wise in 'in'.
//      block_subcarrier_fft(d_fd_data, p_in);
//      block_upsample_and_filter(d_filtered, d_fd_data);
//      combine_subcarriers(d_ifft_in, d_filtered);
      modulate_block(d_ifft_in, p_in);
      fftwf_execute((fftwf_plan) d_ifft_plan);
      memcpy(p_out, d_ifft_out, sizeof(gfdm_complex) * d_n_timeslots * d_n_subcarriers);
    }

  } /* namespace gfdm */
} /* namespace gr */

