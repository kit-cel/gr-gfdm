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

#ifndef INCLUDED_GFDM_RECEIVER_H
#define INCLUDED_GFDM_RECEIVER_H
#endif

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

  //Initialize input FFT
  d_in_fft_in =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_fft_len, volk_get_alignment ());
  d_in_fft_out =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_fft_len, volk_get_alignment ());
  d_in_fft_plan = initialize_fft(d_in_fft_out, d_in_fft_in, d_fft_len, true);

  //Initialize IFFT per subcarrier
  d_sc_ifft_in =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());
  d_sc_ifft_out =  (gfdm_complex *) volk_malloc(sizeof (gfdm_complex) * d_n_timeslots, volk_get_alignment ());
  d_sc_ifft_plan = initialize_fft(d_sc_ifft_out, d_sc_ifft_in, d_n_timeslots, false);

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
  ::volk_32fc_s32fc_multiply_32fc(&d_in_fft_in[0],&in[0],static_cast<gfdm_complex>(float(d_N)/float(d_fft_len)),d_fft_len);
  //std::memcpy(&d_in_fft_in[0],&in[0],sizeof(gfdm_complex)*d_fft_len);
  d_in_fft->execute();
  std::vector<gfdm_complex> fft_out(d_fft_len+d_ntimeslots*d_filter_width);
  std::memcpy(&fft_out[0],&d_in_fft_out[0],sizeof(gfdm_complex)*d_fft_len);
  std::memcpy(&fft_out[d_fft_len],&d_in_fft_out[0],sizeof(gfdm_complex)*d_ntimeslots*d_filter_width);
  for (int k=0; k<d_nsubcarrier; k++)
  {
    std::vector<gfdm_complex> sc_postfft(d_ntimeslots*d_filter_width);
    std::vector<gfdm_complex> sc_postfilter(d_ntimeslots*d_filter_width);
    //FFT output is not centered:
    //Subcarrier-Offset = d_fft_len/2 + (d_fft_len-d_N)/2 - ((d_filter_width-1)*(d_ntimeslots))/2 + k*d_ntimeslots ) modulo d_fft_len
    int sc_offset = (d_fft_len/2 + (d_fft_len - d_N)/2 - ((d_filter_width-1)*(d_ntimeslots))/2 + k*d_ntimeslots) % d_fft_len;
    gfdm_complex * sc = &fft_out[sc_offset];
    for (int n=0; n<d_filter_width*d_ntimeslots; n++)
    {
      sc_postfft[n] = sc[(n + (d_ntimeslots*d_filter_width)/2) % (d_ntimeslots*d_filter_width)];
    }
    ::volk_32fc_x2_multiply_32fc(&sc_postfilter[0],
                                 &sc_postfft[0],&d_filter_taps[0],d_ntimeslots*d_filter_width);
    // Only valid for d_filter_width = 2, write own kernel for complex addition
    ::volk_32f_x2_add_32f((float*)&out[k][0],
                          (float*)(&sc_postfilter[0]),(float*)(&sc_postfilter[d_ntimeslots]),2*d_ntimeslots);
  }

}

void
receiver_kernel_cc::demodulate_subcarrier(std::vector< std::vector<gfdm_complex> > &out,
    std::vector< std::vector<gfdm_complex> > &sc_fdomain)
{
  // 4. apply ifft on every filtered and superpositioned subcarrier
  for (int k=0; k<d_nsubcarrier; k++)
  {
    std::vector<gfdm_complex> sc_postifft(d_ntimeslots);
    std::memcpy(&d_sc_ifft_in[0],&sc_fdomain[k][0],sizeof(gfdm_complex)*d_ntimeslots);
    d_sc_ifft->execute();
    ::volk_32fc_s32fc_multiply_32fc(&out[k][0],&d_sc_ifft_out[0],static_cast<gfdm_complex>(1.0/(float)d_ntimeslots),d_ntimeslots);
  }

}

void
receiver_kernel_cc::serialize_output(gfdm_complex out[],
                                     std::vector< std::vector<gfdm_complex> > &sc_symbols)
{
  for (int k=0; k<d_nsubcarrier; k++)
  {
    for (int m=0; m<d_ntimeslots; m++)
    {
      out[k+m*d_nsubcarrier] = sc_symbols[k][m];
    }
  }
}

void
receiver_kernel_cc::gfdm_work(gfdm_complex out[],const gfdm_complex in[], int ninput_items, int noutputitems)
{
  filter_superposition(d_sc_fdomain,in);
  demodulate_subcarrier(d_sc_symbols,d_sc_fdomain);
  serialize_output(out,d_sc_symbols);
}

} /* namespace filter */
} /* namespace gr */



