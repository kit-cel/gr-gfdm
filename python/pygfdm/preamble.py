#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Andrej Rode, Johannes Demel.
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

'''
Relevant paper section
[0] A synchronization technique for generalized frequency division multiplexing
[1] Improved Preamble-Aided Timing Estimation for OFDM Systems

COMMENT
[1] describes the algorithm while [0] additionally explains GFDM preamble generation.

'''
from __future__ import print_function, division, unicode_literals
import numpy as np
import commpy as cp
from .modulation import gfdm_tx, gfdm_tx_fft2
from .utils import get_random_qpsk, generate_seed, calculate_signal_energy
from .mapping import map_to_waveform_resources
from .mapping import get_data_matrix
from .cyclic_prefix import add_cyclic_prefix, get_raised_cosine_ramp, get_window_len, pinch_block, add_cyclic_starfix
from .filters import get_frequency_domain_filter
from .gfdm_modulation import gfdm_modulate_block


def sync_symbol(filtertype, alpha, K, n_mod, N):
    '''
        Generate Schmidls Training Symbols to achieve Receiver Synchronisation
        K: should be an odd number
        Process:
            * First Symbol: Transmit PN-Sequence on all even frequencies while zeros on the odd
            frequencies. Constant signal energy -> Multiply every Symbol with sqrt(2)

            * Second Symbol: Transmit PN-Sequence on all odd frequencies and another PN-Sequence
            on the even frequencies
    '''

    pn_order = 14
    pn_seed = '00101101010010'
    pn_mask = '01001110100111'
    if int(np.floor(K / 2.0)) % 2:
        n_even_freq = int(np.floor(K / 2.0))
    else:
        n_even_freq = int(np.ceil(K / 2.0))
    seq_length = n_even_freq * n_mod
    sym_sequence = np.zeros(K, dtype='complex')
    pn_sequence = cp.pnsequence(pn_order, pn_seed, pn_mask, seq_length)
    qam_mod = cp.modulation.QAMModem(2 ** n_mod)
    qam_sequence = qam_mod.modulate(pn_sequence)
    for i in xrange(len(sym_sequence)):
        if not i % 2:
            sym_sequence[i] = qam_sequence[i / 2]
    ifft_sequence = np.fft.ifftshift(sym_sequence)
    output = gfdm_tx(ifft_sequence, filtertype, alpha, 1, K, N)
    # output = np.fft.ifft(np.sqrt(2)*ifft_sequence)

    return output


def sync_symbol2(filtertype, alpha, K, L, n_mod):
    pn_order = 14
    pn_seed = '01001000111011'
    pn_mask = '01001101001110'
    seq_length = K * n_mod
    pn_sequence = cp.pnsequence(pn_order, pn_seed, pn_mask, seq_length)
    qam_mod = cp.modulation.QAMModem(2 ** n_mod)
    qam_sequence = qam_mod.modulate(pn_sequence)
    output = gfdm_tx_fft2(np.tile(qam_sequence, 2), filtertype, alpha, 2, K, 2, 1)
    return output


def mapped_preamble(seed, filtertype, alpha, active_subcarriers, fft_len, subcarrier_map, overlap, cp_len, ramp_len):
    pn_vals = get_random_qpsk(active_subcarriers, seed)
    pn_sym = map_to_waveform_resources(pn_vals, active_subcarriers, fft_len, subcarrier_map)
    return generate_sync_symbol(pn_sym, filtertype, alpha, fft_len, overlap, cp_len, ramp_len)


def symmetric_mapped_preamble(seed, filtertype, alpha, active_subcarriers, fft_len, subcarrier_map, overlap, cp_len, ramp_len):
    pn_vals = get_random_qpsk(active_subcarriers // 2, seed)
    pn_vals = np.concatenate((pn_vals, np.conj(pn_vals[::-1])))
    pn_sym = map_to_waveform_resources(pn_vals, active_subcarriers, fft_len, subcarrier_map)
    return generate_sync_symbol(pn_sym, filtertype, alpha, fft_len, overlap, cp_len, ramp_len), pn_vals


def get_sync_symbol(pn_symbols, H, K, L, cp_len, ramp_len):
    M = 2  # fixed for preamble
    pn_symbols = np.concatenate((pn_symbols, pn_symbols))
    D = get_data_matrix(pn_symbols, K, group_by_subcarrier=True)  # careful here! group by subcarrier is correct!
    symbol = x_symbol = gfdm_modulate_block(D, H, M, K, L, compat_mode=False)
    symbol = add_cyclic_starfix(symbol, cp_len, ramp_len)
    window_ramp = get_raised_cosine_ramp(ramp_len, get_window_len(cp_len, M, K, ramp_len))
    symbol = pinch_block(symbol, window_ramp)
    return symbol, x_symbol


def generate_sync_symbol(pn_symbols, filtertype, alpha, K, L, cp_len, ramp_len):
    H = get_frequency_domain_filter(filtertype, alpha, 2, K, L)
    filter_energy = calculate_signal_energy(H)
    scaling_factor = 1. / np.sqrt(filter_energy / 2)
    H = H * scaling_factor
    return get_sync_symbol(pn_symbols, H, K, L, cp_len, ramp_len)


def check_preamble_properties(preamble, x_preamble):
    x_1st = x_preamble[0:len(x_preamble) // 2]
    x_2nd = x_preamble[-len(x_preamble) // 2:]
    if not np.all(np.abs(x_1st - x_2nd) < 1e-12):
        print(np.abs(x_1st - x_2nd))
        raise ValueError('preamble timeslots do not repeat!')
    from correlation import cross_correlate_naive, auto_correlate_halfs
    from utils import calculate_signal_energy
    x_ampl = np.sqrt(calculate_signal_energy(x_preamble))
    preamble *= x_ampl
    x_preamble *= x_ampl
    x_energy = calculate_signal_energy(x_preamble)
    if np.abs(2. * auto_correlate_halfs(x_preamble) / x_energy) -1. > 1e-10:
        raise ValueError('auto correlating halfs of preamble fails!')

    print('normalized preamble xcorr val: ', np.correlate(x_preamble, x_preamble) / x_energy)
    print('windowed normalized preamble: ', np.correlate(preamble[-len(x_preamble):], x_preamble) / x_energy)
    fxc = np.correlate(preamble, x_preamble, 'full') / x_energy
    vxc = np.correlate(preamble, x_preamble, 'valid') / x_energy
    nxc = cross_correlate_naive(preamble, x_preamble) / x_energy
    import matplotlib.pyplot as plt

    plt.plot(np.abs(fxc))
    plt.plot(np.abs(vxc))
    plt.plot(np.abs(nxc))
    plt.show()


def main():
    np.set_printoptions(precision=4, suppress=True)
    seed = generate_seed('awesome')
    fft_len = 32
    active_subcarriers = 24
    subcarrier_map = np.arange(fft_len)
    subcarrier_map = np.concatenate((subcarrier_map[0:active_subcarriers//2], subcarrier_map[-active_subcarriers//2:]))
    preamble, x_preamble = mapped_preamble(seed, 'rrc', .1, active_subcarriers, fft_len, subcarrier_map, 2, fft_len // 2, fft_len // 8)
    check_preamble_properties(preamble, x_preamble)


if __name__ == '__main__':
    main()
