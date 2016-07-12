/* -*- c++ -*- */
/*
 * Copyright 2016 Andrej Rode.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include <gfdm/receiver_kernel_cc.h>
#include <iostream>
#include <volk/volk.h>
#include <string.h>

namespace gr
{
  namespace gfdm
  {

    receiver_kernel_cc::receiver_kernel_cc(int n_subcarriers,
                                           int n_timeslots,
                                           int overlap,
                                           std::vector<gfdm_complex> frequency_taps):
      d_n_subcarriers(n_subcarriers),
      d_n_timeslots(n_timeslots),
      d_fft_len(n_timeslots * n_subcarriers),
      d_overlap(overlap)
    {
      if (int(frequency_taps.size()) != n_timeslots * overlap)
      {
        std::stringstream sstm;
        sstm << "number of frequency taps(" << frequency_taps.size() << ") MUST be equal to n_timeslots(";
        sstm << n_timeslots << ") * overlap(" << overlap << ") = " << n_timeslots * overlap << "!";
        std::string err_str = sstm.str();
        //" MUST be equal to n_timeslots * overlap!";
        throw std::invalid_argument(err_str.c_str());
      }
      d_filter_taps = (gfdm_complex *) volk_malloc (sizeof (gfdm_complex) * n_timeslots * overlap, volk_get_alignment ());
      memcpy(d_filter_taps, &frequency_taps[0], sizeof(gfdm_complex) * n_timeslots * overlap);
      d_ic_filter_taps = (gfdm_complex *) volk_malloc (sizeof (gfdm_complex) * n_timeslots, volk_get_alignment ());

      //Initialize input FFT
      d_in_fft_in =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_fft_len, volk_get_alignment ());
      d_in_fft_out =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_fft_len, volk_get_alignment ());
      d_in_fft_plan = initialize_fft(d_in_fft_out, d_in_fft_in, d_fft_len, true);

      //Initialize IFFT per subcarrier
      d_sc_ifft_in =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());
      d_sc_ifft_out =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());
      d_sc_ifft_plan = initialize_fft(d_sc_ifft_out, d_sc_ifft_in, d_n_timeslots, false);

      //Initialize FFT per subcarrier (IC)
      d_sc_fft_in =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());
      d_sc_fft_out =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());
      d_sc_fft_plan = initialize_fft(d_sc_ifft_out, d_sc_ifft_in, d_n_timeslots, true);

      d_sc_postfilter = (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());

      //Initialize vector of vectors for temporary subcarrier data
      d_sc_fdomain.resize(d_n_subcarriers);
      for (std::vector< std::vector<gfdm_complex> >::iterator it = d_sc_fdomain.begin(); it != d_sc_fdomain.end(); ++it)
      {
        it->resize(d_n_timeslots);
      }
      d_sc_symbols.resize(d_n_subcarriers);
      for (std::vector< std::vector<gfdm_complex> >::iterator it = d_sc_symbols.begin(); it != d_sc_symbols.end(); ++it)
      {
        it->resize(d_n_timeslots);
      }

      //Calculate IC_Filter_taps d_overlaps MUST be greater or equal to 2
      ::volk_32fc_x2_multiply_32fc(&d_ic_filter_taps[0],&frequency_taps[0],&frequency_taps[n_timeslots],n_timeslots);
    }

    receiver_kernel_cc::~receiver_kernel_cc()
    {
      volk_free(d_filter_taps);

      fftwf_destroy_plan (d_in_fft_plan);
      volk_free(d_in_fft_in);
      volk_free(d_in_fft_out);

      fftwf_destroy_plan (d_sc_ifft_plan);
      volk_free(d_sc_ifft_in);
      volk_free(d_sc_ifft_out);
      volk_free(d_sc_postfilter);
    }

    fftwf_plan
    receiver_kernel_cc::initialize_fft(gfdm_complex *out_buf, gfdm_complex *in_buf, const int fft_size, bool forward)
    {
      std::string filename(getenv("HOME"));
      filename += "/.gr_fftw_wisdom";
      FILE *fpr = fopen (filename.c_str(), "r");
      if (fpr != 0)
      {
        int r = fftwf_import_wisdom_from_file (fpr);
        fclose (fpr);
      }

      fftwf_plan plan = fftwf_plan_dft_1d(fft_size,
                                          reinterpret_cast<fftwf_complex *>(in_buf),
                                          reinterpret_cast<fftwf_complex *>(out_buf),
                                          forward ? FFTW_FORWARD : FFTW_BACKWARD,
                                          FFTW_MEASURE);

      FILE *fpw = fopen (filename.c_str(), "w");
      if (fpw != 0)
      {
        fftwf_export_wisdom_to_file (fpw);
        fclose (fpw);
      }
      return plan;
    }


    void
    receiver_kernel_cc::filter_superposition(std::vector< std::vector<gfdm_complex> > &out,
        const gfdm_complex in[] )
    {
      //No scaling needed -> AGC
      //::volk_32fc_s32fc_multiply_32fc(&d_in_fft_in[0],&in[0],static_cast<gfdm_complex>(float(d_fft_len)/float(d_fft_len)),d_fft_len);
      memset(&d_in_fft_in[0],0x00,sizeof(gfdm_complex)*d_fft_len);
      memcpy(&d_in_fft_in[0],&in[0],sizeof(gfdm_complex)*d_fft_len);
      //To Frequency Domain!
      fftwf_execute(d_in_fft_plan);
      //Append additional d_n_timeslots*d_overlap symbols from beginning to fft_out vector
      //To have 'whole' subcarriers in memory
      for (int k=0; k<d_n_subcarriers; ++k)
      {
        //FFT output is not centered:
        //Subcarrier-Offset = d_fft_len/2 + (d_fft_len-d_fft_len)/2 - ((d_overlap-1)*(d_n_timeslots))/2 + k*d_n_timeslots ) modulo d_fft_len
        // First subcarrier is at the beginning!
        ::memset(&out[k][0],0x00,sizeof(gfdm_complex)*d_n_timeslots);
        for (int i = 0; i<d_overlap; ++i)
        {
          int src_part_pos = ((k + i + d_n_subcarriers - (d_overlap/2)) % d_n_subcarriers) * d_n_timeslots;
          int target_part_pos = ((i+d_overlap / 2)%d_overlap) *d_n_timeslots;
          //Filter every part length d_n_timeslots with appropriate filter_taps
          ::volk_32fc_x2_multiply_32fc(d_sc_postfilter, d_filter_taps+target_part_pos,d_in_fft_out+src_part_pos, d_n_timeslots);
          //Superposition parts in out[k]
          ::volk_32f_x2_add_32f((float*)&out[k][0],
                              (float*)(&out[k][0]),(float*)(&d_sc_postfilter[0]),d_n_timeslots);
        }
      }

    }

    void
    receiver_kernel_cc::demodulate_subcarrier(std::vector< std::vector<gfdm_complex> > &out,
        std::vector< std::vector<gfdm_complex> > &sc_fdomain)
    {
      // 4. apply ifft on every filtered and superpositioned subcarrier
      for (int k=0; k<d_n_subcarriers; ++k)
      {
        memcpy(&d_sc_ifft_in[0],&sc_fdomain[k][0],sizeof(gfdm_complex)*d_n_timeslots);
        fftwf_execute(d_sc_ifft_plan);
        //Scale afterwards if required
        ::volk_32fc_s32fc_multiply_32fc(&out[k][0],&d_sc_ifft_out[0],static_cast<gfdm_complex>(1.0/d_n_subcarriers),d_n_timeslots);
        //memcpy(&out[k][0],&d_sc_ifft_out[0],sizeof(gfdm_complex)*d_n_timeslots);
      }

    }

    void
    receiver_kernel_cc::serialize_output(gfdm_complex out[],
                                         std::vector< std::vector<gfdm_complex> > &sc_symbols)
    {
      for (int k=0; k<d_n_subcarriers; ++k)
      {
        memcpy(out+k*d_n_timeslots,&sc_symbols[k][0],sizeof(gfdm_complex)*d_n_timeslots);
      }
    }
    
    void
    receiver_kernel_cc::remove_sc_interference(std::vector< std::vector<gfdm_complex> > &sc_symbols, std::vector< std::vector<gfdm_complex> > &sc_fdomain)
    {
      std::vector< std::vector<gfdm_complex> > prev_sc_symbols = sc_symbols;
      std::vector<gfdm_complex> sc_tmp(d_n_timeslots);
      //Alle Untertraeger
      for (int k=0; k<d_n_subcarriers; k++)
      {
        if(d_n_timeslots*d_n_subcarriers < (d_n_timeslots*d_n_subcarriers) && ((k==0) || (k==d_n_subcarriers-1)))
        {
        //Sum up neighboring subcarriers and then filter with ic_filter_taps
        ::volk_32f_x2_add_32f((float*)&d_sc_fft_in[0],(float*)&prev_sc_symbols[(((k-1) % d_n_subcarriers) + d_n_subcarriers) % d_n_subcarriers][0],(float*)&prev_sc_symbols[((k+1)% d_n_subcarriers)][0],2*d_n_timeslots);
        fftwf_execute(d_sc_fft_plan);
        ::volk_32fc_x2_multiply_32fc(&sc_symbols[k][0],&d_ic_filter_taps[0],&d_sc_fft_out[0],d_n_timeslots);
        ::volk_32f_x2_subtract_32f((float*)&sc_tmp[0],(float*)&sc_fdomain[k][0],(float*)&sc_symbols[k][0],2*d_n_timeslots);
        memcpy(&sc_symbols[k][0],&sc_tmp[0],sizeof(gfdm_complex)*d_n_timeslots);

      }

    }
    }

    void
    receiver_kernel_cc::generic_work(gfdm_complex* out,const gfdm_complex* in)
    {
      filter_superposition(d_sc_fdomain,in);
      demodulate_subcarrier(d_sc_symbols,d_sc_fdomain);
      serialize_output(out,d_sc_symbols);
    }

  } /* namespace filter */
} /* namespace gr */






