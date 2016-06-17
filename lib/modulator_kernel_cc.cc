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

#include <gfdm/modulator_kernel_cc.h>
#include <iostream>
#include <volk/volk.h>
#include <string.h>
//#include <sstream>


namespace gr {
  namespace gfdm {

    modulator_kernel_cc::modulator_kernel_cc(int n_timeslots, int n_subcarriers, int overlap, std::vector<gfdm_complex> frequency_taps):
      d_n_timeslots(n_timeslots), d_n_subcarriers(n_subcarriers), d_ifft_len(n_timeslots * n_subcarriers), d_overlap(overlap)
    {
      if (int(frequency_taps.size()) != n_timeslots * overlap){
        std::stringstream sstm;
        sstm << "number of frequency taps(" << frequency_taps.size() << ") MUST be equal to n_timeslots(";
        sstm << n_timeslots << ") * overlap(" << overlap << ") = " << n_timeslots * overlap << "!";
        std::string err_str = sstm.str();
                //" MUST be equal to n_timeslots * overlap!";
        throw std::invalid_argument(err_str.c_str());
      }
      d_filter_taps = (gfdm_complex *) volk_malloc (sizeof (gfdm_complex) * n_timeslots * overlap, volk_get_alignment ());
      memcpy(d_filter_taps, &frequency_taps[0], sizeof(gfdm_complex) * n_timeslots * overlap);

      // first create input and output buffers for a new FFTW plan.
      d_sub_fft_in = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots, volk_get_alignment ());
      d_sub_fft_out = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots, volk_get_alignment ());
      d_sub_fft_plan = initialize_fft(d_sub_fft_out, d_sub_fft_in, n_timeslots, true);

      d_filtered = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * n_timeslots, volk_get_alignment ());

      d_ifft_in = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_ifft_len, volk_get_alignment ());
      d_ifft_out = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_ifft_len, volk_get_alignment ());
      d_ifft_plan = initialize_fft(d_ifft_out, d_ifft_in, n_timeslots * n_subcarriers, false);
    }

    modulator_kernel_cc::~modulator_kernel_cc()
    {
      volk_free(d_filter_taps);
      fftwf_destroy_plan (d_sub_fft_plan);
      volk_free(d_sub_fft_in);
      volk_free(d_sub_fft_out);

      volk_free(d_filtered);

      fftwf_destroy_plan(d_ifft_plan);
      volk_free(d_ifft_in);
      volk_free(d_ifft_out);
    }


    fftwf_plan
    modulator_kernel_cc::initialize_fft(gfdm_complex *out_buf, gfdm_complex *in_buf, const int fft_size, bool forward)
    {
      std::string filename(getenv("HOME"));
      filename += "/.gr_fftw_wisdom";
      FILE *fpr = fopen (filename.c_str(), "r");
      if (fpr != 0){
        int r = fftwf_import_wisdom_from_file (fpr);
        fclose (fpr);
      }

      fftwf_plan plan = fftwf_plan_dft_1d(fft_size,
                                      reinterpret_cast<fftwf_complex *>(in_buf),
                                      reinterpret_cast<fftwf_complex *>(out_buf),
                                      forward ? FFTW_FORWARD : FFTW_BACKWARD,
                                      FFTW_MEASURE);

      FILE *fpw = fopen (filename.c_str(), "w");
      if (fpw != 0){
        fftwf_export_wisdom_to_file (fpw);
        fclose (fpw);
      }
      return plan;
    }


    void
    modulator_kernel_cc::generic_work(gfdm_complex* p_out, const gfdm_complex* p_in)
    {
      // assume data symbols are stacked subcarrier-wise in 'in'.
      const int part_len = std::min(d_n_timeslots * d_overlap / 2, d_n_timeslots);

      // make sure we don't sum up old results.
      memset(d_ifft_in, 0x00, sizeof (gfdm_complex) * d_ifft_len);

      // perform modulation for each subcarrier separately
      for(int k = 0; k < d_n_subcarriers; ++k){
        // get items into subcarrier FFT
        memcpy(d_sub_fft_in, p_in, sizeof(gfdm_complex) * d_n_timeslots);
        fftwf_execute(d_sub_fft_plan);

        // handle each part separately. The length of a part should always be d_n_timeslots.
        // FIXME: Assumption and algorithm will probably fail for d_overlap = 1 (Should never be used though).
        for (int i = 0; i < d_overlap; ++i) {
          // calculate positions for next part to handle
          int src_part_pos = ((i + d_overlap / 2) % d_overlap) * d_n_timeslots;
          int target_part_pos = ((k + i + d_n_subcarriers - (d_overlap / 2)) % d_n_subcarriers) * d_n_timeslots;
          // perform filtering operation!
          volk_32fc_x2_multiply_32fc(d_filtered, d_sub_fft_out, d_filter_taps + src_part_pos, d_n_timeslots);
          // add generated part at correct position.
          volk_32f_x2_add_32f((float*) (d_ifft_in + target_part_pos), (float*) (d_ifft_in + target_part_pos), (float*) d_filtered, 2 * part_len);
        }
        p_in += d_n_timeslots;
      }

      // Back to time domain!
      fftwf_execute(d_ifft_plan);
//      memcpy(p_out, d_ifft_out, sizeof(gfdm_complex) * d_ifft_len);
      volk_32fc_s32fc_multiply_32fc(p_out, d_ifft_out, gfdm_complex(1.0 / d_ifft_len, 0), d_ifft_len);
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

  } /* namespace gfdm */
} /* namespace gr */

