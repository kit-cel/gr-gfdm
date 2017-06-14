#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
import scipy.signal as signal
from fractions import gcd

import synchronization as sync
import preamble
import mapping
import converter
import utils
import gfdm_modulation as gmod
import cyclic_prefix

import matplotlib.pyplot as plt


def generate_reference_frame(timeslots, subcarriers, active_subcarriers, cp_len, cs_len):
    p_seed = utils.generate_seed('awesome preamble')
    f_seed = utils.generate_seed('awesome frame')
    subcarrier_map = mapping.get_subcarrier_map(subcarriers, active_subcarriers, dc_free=True)
    alpha = .5
    overlap = 2
    frame_preamble, x_preamble = preamble.mapped_preamble(p_seed, 'rrc', alpha, active_subcarriers, subcarriers, subcarrier_map, overlap, cp_len, cs_len)
    d = utils.get_random_qpsk(timeslots * active_subcarriers, f_seed)
    d_frame = gmod.modulate_mapped_gfdm_block(d, timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)
    symbol = cyclic_prefix.add_cyclic_starfix(d_frame, cp_len, cs_len)
    window_ramp = cyclic_prefix.get_raised_cosine_ramp(cs_len, cyclic_prefix.get_window_len(cp_len, timeslots, subcarriers, cs_len))
    d_frame = cyclic_prefix.pinch_block(symbol, window_ramp)
    return np.concatenate((frame_preamble, d_frame))


def main():
    np.set_printoptions(precision=2, linewidth=150)
    seed = abs(hash('awesome')) % (2 ** 32)
    active_subcarriers = 52
    fft_len = 64
    subcarrier_map = mapping.get_subcarrier_map(fft_len, active_subcarriers, dc_free=True)
    print(subcarrier_map)
    frame_preamble, x_preamble = preamble.mapped_preamble(seed, 'rrc', .5, active_subcarriers, fft_len, subcarrier_map, 2, 64, 16)
    filename = '/home/demel/iq_samples/gfdm_samples_part0.dat'

    frame = converter.load_gr_iq_file(filename)[25000:225000]
    plt.semilogy(*signal.welch(frame))
    # ref_frame = generate_reference_frame(9, fft_len, active_subcarriers, 64, 16)
    # ref_frame = np.concatenate((np.zeros(len(ref_frame) * 10), ref_frame))
    # ac = sync.auto_correlate_signal(ref_frame, fft_len)
    # plt.plot(np.abs(ref_frame))
    # plt.plot(np.abs(ac))
    # plt.show()
    # ref_frame = converter.convert_to_cf64(ref_frame)
    # ref_frame.tofile('/home/demel/iq_samples/gfdm_reference_frame.dat')


    # print('M * K', 9 * 64, ' + 2 * K', 2 * 64, ' + 2CP + 2CS', 2 * 64 + 2 * 16)
    # print(len(frame))
    # ac = sync.auto_correlate_signal(frame, subcarriers)
    #
    # plt.plot(np.abs(ac))
    # plt.plot(np.abs(frame) * 10)

    # p_frame = np.concatenate((np.zeros(1000), frame_preamble, np.zeros(1000)))
    # ac = sync.auto_correlate_signal(p_frame, fft_len)
    #
    # # plt.plot(np.abs(ac))
    # xc = np.correlate(frame, x_preamble)
    # plt.plot(np.abs(xc))
    plt.show()

if __name__ == '__main__':
    main()
