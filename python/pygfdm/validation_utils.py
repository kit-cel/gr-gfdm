#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Johannes Demel.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

import numpy as np
import scipy.signal as signal
import utils
import mapping
import preamble
import gfdm_modulation
import cyclic_prefix
import filters
# import synchronization as sync

import matplotlib.pyplot as plt


def generate_reference_frame(timeslots, subcarriers, active_subcarriers, cp_len, cs_len, alpha=.2):
    p_seed = utils.generate_seed('awesome preamble')
    f_seed = utils.generate_seed('awesome frame')
    subcarrier_map = mapping.get_subcarrier_map(subcarriers, active_subcarriers, dc_free=True)
    overlap = 2
    frame_preamble, x_preamble = preamble.mapped_preamble(p_seed, 'rrc', alpha, active_subcarriers, subcarriers, subcarrier_map, overlap, cp_len, cs_len)
    d = utils.get_random_qpsk(timeslots * active_subcarriers, f_seed)
    d_frame = mod_frame = gfdm_modulation.modulate_mapped_gfdm_block(d, timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)
    symbol = cyclic_prefix.add_cyclic_starfix(d_frame, cp_len, cs_len)
    window_ramp = cyclic_prefix.get_raised_cosine_ramp(cs_len, cyclic_prefix.get_window_len(cp_len, timeslots, subcarriers, cs_len))
    d_frame = cyclic_prefix.pinch_block(symbol, window_ramp)

    H = filters.get_frequency_domain_filter('rrc', alpha, timeslots, subcarriers, overlap)
    return np.concatenate((frame_preamble, d_frame)), mod_frame, x_preamble, d, H


def generate_sc_qpsk_frame(timeslots, subcarriers, active_subcarriers, cp_len, cs_len, alpha=.2):
    p_seed = utils.generate_seed('awesome preamble')
    f_seed = utils.generate_seed('awesome frame')
    subcarrier_map = mapping.get_subcarrier_map(subcarriers, active_subcarriers, dc_free=True)
    overlap = 2
    frame_preamble, x_preamble = preamble.mapped_preamble(p_seed, 'rrc', alpha, active_subcarriers, subcarriers, subcarrier_map, overlap, cp_len, cs_len)
    d = .2 * utils.get_random_qpsk(timeslots * subcarriers // 4, f_seed)
    d = signal.resample(d, len(d) * 4)
    # d_frame = mod_frame = gfdm_modulation.modulate_mapped_gfdm_block(d, timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)
    symbol = cyclic_prefix.add_cyclic_starfix(d, cp_len, cs_len)

    window_ramp = cyclic_prefix.get_raised_cosine_ramp(cs_len, cyclic_prefix.get_window_len(cp_len, timeslots, subcarriers, cs_len))
    d_frame = cyclic_prefix.pinch_block(symbol, window_ramp)

    H = filters.get_frequency_domain_filter('rrc', alpha, timeslots, subcarriers, overlap)
    return np.concatenate((frame_preamble, d_frame)), d, x_preamble, d, H


def generate_integrated_frame(timeslots, subcarriers, active_subcarriers, cp_len, cs_len, alpha=.2):
    p_seed = utils.generate_seed('awesome preamble')
    f_seed = utils.generate_seed('awesome frame')
    subcarrier_map = mapping.get_subcarrier_map(subcarriers, active_subcarriers, dc_free=True)
    overlap = 2
    p, p_vals = preamble.symmetric_mapped_preamble(p_seed, 'rrc', alpha, active_subcarriers, subcarriers, subcarrier_map, overlap, cp_len, cs_len)
    frame_preamble, x_preamble = p
    p = gfdm_modulation.modulate_mapped_gfdm_block(np.concatenate((p_vals, p_vals, np.zeros((timeslots - 2) * active_subcarriers))), timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)
    x_preamble = p[0:len(x_preamble)]

    d = utils.get_random_qpsk((timeslots - 4) * active_subcarriers, f_seed)
    d = np.tile(p_vals, timeslots)
    # d = np.concatenate((p_vals, p_vals, d, p_vals, p_vals))
    # d = utils.get_random_qpsk((timeslots - 2) * active_subcarriers, f_seed)
    # d = np.concatenate((p_vals, p_vals, d))


    d_frame = mod_frame = gfdm_modulation.modulate_mapped_gfdm_block(d, timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)

    symbol = cyclic_prefix.add_cyclic_starfix(d_frame, cp_len, cs_len)

    window_ramp = cyclic_prefix.get_raised_cosine_ramp(cs_len, cyclic_prefix.get_window_len(cp_len, timeslots, subcarriers, cs_len))
    # d_frame = cyclic_prefix.pinch_block(symbol, window_ramp)

    H = filters.get_frequency_domain_filter('rrc', alpha, timeslots, subcarriers, overlap)
    return p, mod_frame, x_preamble, d, H


# def synchronize_integrated(frame, ref_frame, x_preamble, fft_len, cp_len):
#     samp_rate = 12.5e6
#     ac = sync.auto_correlate_signal(frame, fft_len)
#
#     nm = np.argmax(np.abs(ac[0:len(ac) // 2]))
#     print('AC start: ', nm)
#     # cfo = 2 * np.angle(ac[nm]) / (2. * np.pi)
#     cfo = np.angle(ac[nm]) / (2. * np.pi)
#     print('CFO:', cfo, cfo * samp_rate / fft_len)
#
#     phase_inc = sync.cfo_to_phase_increment(-cfo, fft_len)
#     wave = sync.complex_sine(phase_inc, len(frame), 0.0)
#     # frame *= wave
#
#     ac = sync.auto_correlate_signal(frame, fft_len)
#     cfo = np.angle(ac[nm]) / (2. * np.pi)
#     print('CFO:', cfo, cfo * samp_rate / fft_len)
#     ac = np.roll(ac, cp_len)
#
#     xc = np.correlate(frame, x_preamble, 'valid')
#     cc = sync.multiply_valid(np.abs(ac), np.abs(xc))
#     nc = np.argmax(np.abs(cc[0:len(cc)//2]))
#     print('correlation frame start:', nc)
#     sample_nc = nc - cp_len
#     print('sample frame start:     ', sample_nc)
#     print('data frame start:       ', nc)
#     phase = np.angle(xc[nc])
#     # phase = 0.0
#     print('phase:', phase)
#     # frame *= np.exp(-1j * phase)
#
#     ref_e = utils.calculate_signal_energy(x_preamble)
#     p = frame[nc:nc + len(x_preamble)]
#     rx_e = utils.calculate_signal_energy(p)
#     agc_factor = np.sqrt(ref_e / rx_e)
#     print('AGC values:', ref_e, rx_e, agc_factor)
#     # frame *= agc_factor
#     sframe = frame[sample_nc:sample_nc + len(ref_frame)]
#     # plt.plot(np.abs(ref_frame))
#     # plt.plot(np.abs(frame))
#     plt.plot(np.abs(ac), label='AC')
#     plt.plot(np.abs(xc), label='XC')
#     plt.plot(cc, label='MM')
#     # # plt.axvline(sample_nc, color='y')
#
#     # print(np.abs(p - x_preamble))
#     # print(np.max(np.abs(p - x_preamble)))
#     # x = p.dot(x_preamble)
#     # print(x, x_preamble.dot(x_preamble), p.dot(p))
#     # plt.scatter(p.real, p.imag)
#     # plt.scatter(x_preamble.real, x_preamble.imag)
#     plt.legend()
#     plt.show()
#     return sframe


def main():
    np.set_printoptions(precision=2, linewidth=150)
    alpha = .2
    active_subcarriers = 52
    timeslots = 9
    fft_len = 64
    cp_len = fft_len // 2
    cs_len = cp_len // 2
    subcarrier_map = mapping.get_subcarrier_map(fft_len, active_subcarriers, dc_free=True)
    ref_frame, modulated_frame, x_preamble, data, freq_filter_taps = generate_integrated_frame(timeslots, fft_len, active_subcarriers, cp_len, cs_len, alpha)
    test_frame = np.concatenate((.001 * utils.get_random_samples(1000), ref_frame, .001 * utils.get_random_samples(1000)))
    # sframe = synchronize_integrated(test_frame, ref_frame, x_preamble, fft_len, cp_len)
    # return


if __name__ == '__main__':
    main()
