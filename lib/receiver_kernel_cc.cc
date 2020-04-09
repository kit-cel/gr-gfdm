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
#include <string.h>
#include <volk/volk.h>
#include <iostream>

namespace gr {
namespace gfdm {

receiver_kernel_cc::receiver_kernel_cc(int n_timeslots,
                                       int n_subcarriers,
                                       int overlap,
                                       std::vector<gfdm_complex> frequency_taps)
    : d_n_subcarriers(n_subcarriers),
      d_n_timeslots(n_timeslots),
      d_block_len(n_timeslots * n_subcarriers),
      d_overlap(overlap)
{
    if (int(frequency_taps.size()) != n_timeslots * overlap) {
        std::stringstream sstm;
        sstm << "number of frequency taps(" << frequency_taps.size()
             << ") MUST be equal to n_timeslots(";
        sstm << n_timeslots << ") * overlap(" << overlap
             << ") = " << n_timeslots * overlap << "!";
        throw std::invalid_argument(sstm.str().c_str());
    }
    if (d_overlap < 2) {
        std::stringstream sstm;
        sstm << "overlap MUST be greater or equal 2";
        throw std::invalid_argument(sstm.str().c_str());
    }
    d_filter_taps = (gfdm_complex*)volk_malloc(
        sizeof(gfdm_complex) * n_timeslots * overlap, volk_get_alignment());
    initialize_taps_vector(d_filter_taps, frequency_taps, n_timeslots);

    d_ic_filter_taps = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * n_timeslots,
                                                  volk_get_alignment());
    memset(d_ic_filter_taps, 0x00, sizeof(gfdm_complex) * n_timeslots);
    // Calculate IC_Filter_taps d_overlap MUST be greater or equal to 2,only consider
    // interference of neigboring subcarriers
    ::volk_32fc_x2_multiply_32fc(d_ic_filter_taps,
                                 d_filter_taps,
                                 &d_filter_taps[n_timeslots * (overlap - 1)],
                                 n_timeslots);

    // Initialize input FFT
    d_in_fft_in = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_block_len,
                                             volk_get_alignment());
    d_in_fft_out = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_block_len,
                                              volk_get_alignment());
    d_in_fft_plan = initialize_fft(d_in_fft_out, d_in_fft_in, d_block_len, true);

    // OPTIONAL: Equalize!
    d_equalized = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_block_len,
                                             volk_get_alignment());

    // Initialize IFFT per subcarrier
    d_sc_ifft_in = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_n_timeslots,
                                              volk_get_alignment());
    d_sc_ifft_out = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_n_timeslots,
                                               volk_get_alignment());
    d_sc_ifft_plan = initialize_fft(d_sc_ifft_out, d_sc_ifft_in, d_n_timeslots, false);

    // Initialize FFT per subcarrier (IC)
    d_sc_fft_in = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_n_timeslots,
                                             volk_get_alignment());
    d_sc_fft_out = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_n_timeslots,
                                              volk_get_alignment());
    d_sc_fft_plan = initialize_fft(d_sc_fft_out, d_sc_fft_in, d_n_timeslots, true);

    d_sc_postfilter = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_n_timeslots,
                                                 volk_get_alignment());

    d_sc_filtered = (gfdm_complex*)volk_malloc(sizeof(gfdm_complex) * d_block_len,
                                               volk_get_alignment());
}

receiver_kernel_cc::~receiver_kernel_cc()
{
    volk_free(d_filter_taps);
    volk_free(d_ic_filter_taps);

    fftwf_destroy_plan(d_in_fft_plan);
    volk_free(d_in_fft_in);
    volk_free(d_in_fft_out);

    fftwf_destroy_plan(d_sc_ifft_plan);
    volk_free(d_sc_ifft_in);
    volk_free(d_sc_ifft_out);

    fftwf_destroy_plan(d_sc_fft_plan);
    volk_free(d_sc_fft_in);
    volk_free(d_sc_fft_out);

    volk_free(d_sc_postfilter);
    volk_free(d_sc_filtered);
}

void receiver_kernel_cc::initialize_taps_vector(gfdm_complex* filter_taps,
                                                std::vector<gfdm_complex> frequency_taps,
                                                const int n_timeslots)
{
    gfdm_complex res = gfdm_complex(0.0, 0.0);
    volk_32fc_x2_conjugate_dot_prod_32fc(
        &res, &frequency_taps[0], &frequency_taps[0], frequency_taps.size());
    // std::cout << "BEFORE energy of taps: " << std::abs(res) << std::endl;

    const gfdm_complex scaling_factor =
        gfdm_complex(1. / std::sqrt(std::abs(res) / n_timeslots), 0.0f);
    volk_32fc_s32fc_multiply_32fc(
        d_filter_taps, &frequency_taps[0], scaling_factor, frequency_taps.size());

    volk_32fc_x2_conjugate_dot_prod_32fc(
        &res, d_filter_taps, d_filter_taps, frequency_taps.size());
    // std::cout << "AFTER  energy of taps: " << std::abs(res) << std::endl;
}

std::vector<receiver_kernel_cc::gfdm_complex> receiver_kernel_cc::filter_taps() const
{
    std::vector<gfdm_complex> taps(d_n_timeslots * d_overlap);
    memcpy(&taps[0], d_filter_taps, sizeof(gfdm_complex) * d_n_timeslots * d_overlap);
    return taps;
}

std::vector<receiver_kernel_cc::gfdm_complex> receiver_kernel_cc::ic_filter_taps() const
{
    std::vector<gfdm_complex> taps(d_n_timeslots);
    memcpy(&taps[0], d_ic_filter_taps, sizeof(gfdm_complex) * d_n_timeslots);
    return taps;
}

void receiver_kernel_cc::filter_superposition(std::vector<std::vector<gfdm_complex>>& out,
                                              const gfdm_complex* in)
{
    // No scaling needed -> AGC
    //::volk_32fc_s32fc_multiply_32fc(&d_in_fft_in[0],&in[0],static_cast<gfdm_complex>(float(d_block_len)/float(d_block_len)),d_block_len);
    //      memset(&d_in_fft_in[0], 0x00, sizeof(gfdm_complex) * d_block_len);
    memcpy(&d_in_fft_in[0], &in[0], sizeof(gfdm_complex) * d_block_len);
    // To Frequency Domain!
    fftwf_execute(d_in_fft_plan);
    // Append additional d_n_timeslots*d_overlap symbols from beginning to fft_out vector
    // To have 'whole' subcarriers in memory
    filter_subcarriers_and_downsample_fd(d_sc_filtered, d_in_fft_out);
    for (int k = 0; k < d_n_subcarriers; ++k) {
        // FFT output is not centered:
        // Subcarrier-Offset = d_block_len/2 + (d_block_len-d_block_len)/2 -
        // ((d_overlap-1)*(d_n_timeslots))/2 + k*d_n_timeslots ) modulo d_block_len
        // First subcarrier is at the beginning!
        ::memset(&out[k][0], 0x00, sizeof(gfdm_complex) * d_n_timeslots);
        for (int i = 0; i < d_overlap; ++i) {
            int src_part_pos =
                ((k + i + d_n_subcarriers - (d_overlap / 2)) % d_n_subcarriers) *
                d_n_timeslots;
            int target_part_pos = ((i + d_overlap / 2) % d_overlap) * d_n_timeslots;
            // Filter every part length d_n_timeslots with appropriate filter_taps
            ::volk_32fc_x2_multiply_32fc(d_sc_postfilter,
                                         d_filter_taps + target_part_pos,
                                         d_in_fft_out + src_part_pos,
                                         d_n_timeslots);
            // Superposition parts in out[k]
            ::volk_32f_x2_add_32f((float*)&out[k][0],
                                  (float*)(&out[k][0]),
                                  (float*)d_sc_postfilter,
                                  2 * d_n_timeslots);
        }
    }
}

void receiver_kernel_cc::filter_subcarriers_and_downsample_fd(gfdm_complex* p_out,
                                                              const gfdm_complex* p_in)
{
    ::memset(p_out, 0x00, sizeof(gfdm_complex) * d_n_timeslots * d_n_subcarriers);
    for (int k = 0; k < d_n_subcarriers; ++k) {
        // FFT DC carrier on 0th bin!
        // Subcarrier-Offset = d_block_len/2 + (d_block_len-d_block_len)/2 -
        // ((d_overlap-1)*(d_n_timeslots))/2 + k*d_n_timeslots ) modulo d_block_len
        // First subcarrier is at the beginning!
        for (int i = 0; i < d_overlap; ++i) {
            int src_part_pos =
                ((k + i + d_n_subcarriers - (d_overlap / 2)) % d_n_subcarriers) *
                d_n_timeslots;
            int target_part_pos = ((i + d_overlap / 2) % d_overlap) * d_n_timeslots;
            // Filter every part length d_n_timeslots with appropriate filter_taps
            ::volk_32fc_x2_multiply_32fc(d_sc_postfilter,
                                         d_filter_taps + target_part_pos,
                                         p_in + src_part_pos,
                                         d_n_timeslots);

            ::volk_32f_x2_add_32f(
                (float*)p_out, (float*)p_out, (float*)d_sc_postfilter, 2 * d_n_timeslots);
        }
        p_out += d_n_timeslots;
    }
}

void receiver_kernel_cc::demodulate_subcarrier(
    std::vector<std::vector<gfdm_complex>>& out,
    std::vector<std::vector<gfdm_complex>>& sc_fdomain)
{
    // 4. apply ifft on every filtered and superpositioned subcarrier
    for (int k = 0; k < d_n_subcarriers; ++k) {
        memcpy(&d_sc_ifft_in[0], &sc_fdomain[k][0], sizeof(gfdm_complex) * d_n_timeslots);
        fftwf_execute(d_sc_ifft_plan);
        // Scale afterwards if required
        ::volk_32fc_s32fc_multiply_32fc(&out[k][0],
                                        &d_sc_ifft_out[0],
                                        static_cast<gfdm_complex>(1.0 / d_n_timeslots),
                                        d_n_timeslots);
        // memcpy(&out[k][0],&d_sc_ifft_out[0],sizeof(gfdm_complex)*d_n_timeslots);
    }
}

void receiver_kernel_cc::transform_subcarriers_to_td(gfdm_complex* p_out,
                                                     const gfdm_complex* p_in)
{
    // 4. apply ifft on every filtered and superpositioned subcarrier
    const gfdm_complex scaling_factor = static_cast<gfdm_complex>(1.0 / d_n_timeslots);
    for (int k = 0; k < d_n_subcarriers; ++k) {
        memcpy(d_sc_ifft_in, p_in, sizeof(gfdm_complex) * d_n_timeslots);
        fftwf_execute(d_sc_ifft_plan);
        // Scale afterwards if required
        ::volk_32fc_s32fc_multiply_32fc(
            p_out, d_sc_ifft_out, scaling_factor, d_n_timeslots);
        p_in += d_n_timeslots;
        p_out += d_n_timeslots;
    }
}

void receiver_kernel_cc::serialize_output(
    gfdm_complex out[], std::vector<std::vector<gfdm_complex>>& sc_symbols)
{
    for (int k = 0; k < d_n_subcarriers; ++k) {
        memcpy(out + k * d_n_timeslots,
               &sc_symbols[k][0],
               sizeof(gfdm_complex) * d_n_timeslots);
    }
}

void receiver_kernel_cc::vectorize_2d(std::vector<std::vector<gfdm_complex>>& out_vector,
                                      const gfdm_complex* p_in)
{
    for (int k = 0; k < d_n_subcarriers; ++k) {
        memcpy(&out_vector[k][0],
               p_in + k * d_n_timeslots,
               sizeof(gfdm_complex) * d_n_timeslots);
    }
}

void receiver_kernel_cc::remove_sc_interference(
    std::vector<std::vector<gfdm_complex>>& sc_symbols,
    std::vector<std::vector<gfdm_complex>>& sc_fdomain)
{
    std::vector<std::vector<gfdm_complex>> prev_sc_symbols = sc_symbols;
    std::vector<gfdm_complex> sc_tmp(d_n_timeslots);
    // Alle Untertraeger
    for (int k = 0; k < d_n_subcarriers; k++) {
        // Sum up neighboring subcarriers and then filter with ic_filter_taps
        int prev_sc = (((k - 1) % d_n_subcarriers) + d_n_subcarriers) % d_n_subcarriers;
        int next_sc = (((k + 1) % d_n_subcarriers) + d_n_subcarriers) % d_n_subcarriers;
        ::volk_32f_x2_add_32f((float*)&d_sc_fft_in[0],
                              (float*)&prev_sc_symbols[prev_sc][0],
                              (float*)&prev_sc_symbols[next_sc][0],
                              2 * d_n_timeslots);
        fftwf_execute(d_sc_fft_plan);
        // Multiply resulting symbols stream with ic_filter_taps in fd
        ::volk_32fc_x2_multiply_32fc(
            &sc_symbols[k][0], &d_ic_filter_taps[0], &d_sc_fft_out[0], d_n_timeslots);
        // Subtract calculated interference from subcarrier k
        ::volk_32f_x2_subtract_32f((float*)&sc_symbols[k][0],
                                   (float*)&sc_fdomain[k][0],
                                   (float*)&sc_symbols[k][0],
                                   2 * d_n_timeslots);
    }
}

void receiver_kernel_cc::cancel_sc_interference(gfdm_complex* p_out,
                                                const gfdm_complex* p_td_in,
                                                const gfdm_complex* p_fd_in)
{
    // apply to all subcarriers!
    for (int k = 0; k < d_n_subcarriers; ++k) {
        // Sum up neighboring subcarriers and then filter with ic_filter_taps
        int prev_sc = (k - 1 + d_n_subcarriers) % d_n_subcarriers;
        int next_sc = (k + 1 + d_n_subcarriers) % d_n_subcarriers;
        ::volk_32f_x2_add_32f((float*)d_sc_fft_in,
                              (float*)(p_td_in + prev_sc * d_n_timeslots),
                              (float*)(p_td_in + next_sc * d_n_timeslots),
                              2 * d_n_timeslots);
        fftwf_execute(d_sc_fft_plan);
        // Multiply resulting symbols stream with ic_filter_taps in fd
        ::volk_32fc_x2_multiply_32fc(
            d_sc_postfilter, d_ic_filter_taps, d_sc_fft_out, d_n_timeslots);
        // Subtract calculated interference from subcarrier k
        ::volk_32f_x2_subtract_32f((float*)(p_out + k * d_n_timeslots),
                                   (float*)(p_fd_in + k * d_n_timeslots),
                                   (float*)d_sc_postfilter,
                                   2 * d_n_timeslots);
    }
}

void receiver_kernel_cc::fft_filter_downsample(gfdm_complex* p_out,
                                               const gfdm_complex* p_in)
{
    memcpy(d_in_fft_in, p_in, sizeof(gfdm_complex) * d_block_len);
    fftwf_execute(d_in_fft_plan);
    filter_subcarriers_and_downsample_fd(p_out, d_in_fft_out);
}

void receiver_kernel_cc::fft_equalize_filter_downsample(gfdm_complex* p_out,
                                                        const gfdm_complex* p_in,
                                                        const gfdm_complex* f_eq_in)
{
    memcpy(d_in_fft_in, p_in, sizeof(gfdm_complex) * d_block_len);
    fftwf_execute(d_in_fft_plan);
    volk_32fc_x2_divide_32fc(d_equalized, d_in_fft_out, f_eq_in, d_block_len);
    // volk_32fc_x2_multiply_conjugate_32fc(d_equalized, d_in_fft_out, f_eq_in,
    // d_block_len);
    filter_subcarriers_and_downsample_fd(p_out, d_equalized);
}

void receiver_kernel_cc::generic_work(gfdm_complex* out, const gfdm_complex* in)
{
    fft_filter_downsample(d_sc_filtered, in);
    transform_subcarriers_to_td(out, d_sc_filtered);
}

void receiver_kernel_cc::generic_work_equalize(gfdm_complex* out,
                                               const gfdm_complex* in,
                                               const gfdm_complex* f_eq_in)
{
    fft_equalize_filter_downsample(d_sc_filtered, in, f_eq_in);
    transform_subcarriers_to_td(out, d_sc_filtered);
}

} // namespace gfdm
} /* namespace gr */
