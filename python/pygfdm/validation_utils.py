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


class frame_estimator():
    def __init__(self, x_preamble, fft_len, timeslots, active_subcarriers):
        self._x_preamble = x_preamble
        self._fft_len = fft_len
        self._timeslots = timeslots
        self._inv_freq_x_preamble0 = 1. / np.fft.fft(x_preamble[0:fft_len])
        self._inv_freq_x_preamble1 = 1. / np.fft.fft(x_preamble[fft_len:])

        active_sc = np.arange((self._fft_len - active_subcarriers)//2, (self._fft_len + active_subcarriers)//2+1)
        active_sc = active_sc[3:-3]
        freqs = np.fft.fftfreq(fft_len)
        freqs = np.fft.fftshift(freqs)
        self._active_preamble_freqs = freqs[active_sc]
        self._active_sc = active_sc
        fr_freqs = np.fft.fftfreq(self._fft_len * self._timeslots)
        self._frame_freqs = np.fft.fftshift(fr_freqs)

        g = signal.gaussian(9, 1.0)
        # g /= g.dot(np.ones(len(g)))
        g /= np.sum(g)
        self._p_filter = g

    def _estimate_preamble(self, rx_preamble):
        e0 = np.fft.fft(rx_preamble[0:self._fft_len]) * self._inv_freq_x_preamble0
        e1 = np.fft.fft(rx_preamble[self._fft_len:]) * self._inv_freq_x_preamble1
        H = (e0 + e1) / 2
        return H

    def _filter_preamble_estimate(self, H):
        H[0] = (H[1] + H[-1]) / 2.
        H = np.fft.fftshift(H)

        Ha = H[self._active_sc]
        Hb = np.concatenate((np.repeat(Ha[0], 4), Ha, np.repeat(Ha[-1], 4)))
        Hg = np.correlate(Hb, self._p_filter)
        return Hg

    def _interpolate_frame(self, H):
        Hg = self._filter_preamble_estimate(H)

        H_frame = np.interp(self._frame_freqs, self._active_preamble_freqs, Hg.real) + 1j * np.interp(self._frame_freqs, self._active_preamble_freqs, Hg.imag)
        return np.fft.fftshift(H_frame)

    def estimate_frame(self, rx_preamble):
        H = self._estimate_preamble(rx_preamble)
        return self._interpolate_frame(H)


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
