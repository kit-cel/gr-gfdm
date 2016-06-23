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

import numpy as np
import commpy as cp
import matplotlib.pyplot as plt
import scipy.signal as signal
from modulation import gfdm_tx, gfdm_tx_fft2, gfdm_modulation_matrix
from filters import get_frequency_domain_filter, gfdm_filter_taps
from gfdm_modulation import gfdm_modulate_block, gfdm_modulate_fft
from cyclic_prefix import add_cyclic_prefix, pinch_block, get_raised_cosine_ramp, get_window_len, get_root_raised_cosine_ramp
from mapping import get_data_matrix
from utils import get_random_qpsk, get_complex_noise_vector, calculate_awgn_noise_variance, calculate_average_signal_energy, calculate_signal_energy


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


def sync_product(x, L):
    '''
    Auto-Korrelation der ersten L Samples mit den naechsten L Samples
    '''
    return np.sum([x[i].conj() * x[i + L] for i in xrange(L)])


def sync_iter(x, L, cp):
    '''
    Schrittweise Iteration ueber alle Samples (zunaechst die ersten 2*L Samples)
    Danach ueber die restlichen len(x)-2*L Samples
    '''
    P_d = np.array([])
    P_d = np.append(P_d, sync_product(x, L))
    for i in xrange(len(x) - 2 * L):
        P_d = np.append(
            P_d, P_d[i] + (x[L + i].conj() * x[2 * L + i]) -
                 (x[i].conj() * x[L + i]))
    P_d_out = P_d
    P_d = np.append(np.zeros(cp, dtype='complex'), P_d)
    # Integrate cp-samples to eliminate cp-plateau
    P_di = np.array([])
    for i in xrange(cp, len(x) - 2 * L):
        P_di = np.append(
            P_di, (1.0 / (cp + 1) * np.sum(np.abs(P_d[i - cp:i]) ** 2)))
    return (P_di, P_d_out)


def sync_energy(x, L):
    '''
    Berechnung der Energie der zweiten Haelfte der Sync-Samples -> normieren
    '''
    R_d = np.array([])
    R_d = np.append(R_d, np.sum([np.abs(x[i + L]) ** 2 for i in xrange(L)]))
    for i in xrange(len(x) - 2 * L):
        R_d = np.append(R_d, R_d[-1] + np.abs(x[2 * L + i]) ** 2 - np.abs(x[L + i]) ** 2)
    return R_d


def sync_perform(x, L, cp):
    (P_di, P_d) = sync_iter(x, L, cp)
    # R_d = sync_energy(x, L)
    # M_d = (np.abs(P_d)**2)/(R_d**2)
    return (P_di, P_d)


def sync_CFO(P_d, P_di):
    '''
    Gewinn von d (Beginn der Pilotsequenz) und d_f Frequenzoffset genormt auf 1/T.
    Kann bei nur einem Pilotsymbol nur +/- 1/T betragen.
    '''
    d = np.argmax(P_di)
    dphi = np.angle(P_d[d])
    d_f = dphi / (np.pi)
    print("P_d:{},df:{})".format(P_d[d], d_f))

    return (d, d_f)


def get_sync_symbol(pn_symbols, H, K, L, cp_len, ramp_len):
    M = 2  # fixed for preamble
    pn_symbols = np.concatenate((pn_symbols, pn_symbols))
    D = get_data_matrix(pn_symbols, K, group_by_subcarrier=True) # careful here! group by subcarrier is correct!
    symbol = x_symbol = gfdm_modulate_block(D, H, M, K, L, compat_mode=False)
    symbol = add_cyclic_prefix(symbol, cp_len)
    window_ramp = get_raised_cosine_ramp(ramp_len, get_window_len(cp_len, M, K))
    symbol = pinch_block(symbol, window_ramp)
    return symbol, x_symbol


def generate_sync_symbol(pn_symbols, filtertype, alpha, K, L, cp_len, ramp_len):
    H = get_frequency_domain_filter(filtertype, alpha, 2, K, L)
    return get_sync_symbol(pn_symbols, H, K, L, cp_len, ramp_len)


def get_gfdm_frame(data, alpha, M, K, L, cp_len, ramp_len):
    window_len = get_window_len(cp_len, M, K)
    b = gfdm_modulate_fft(data, alpha, M, K, L)
    b = add_cyclic_prefix(b, cp_len)
    window_ramp = get_root_raised_cosine_ramp(ramp_len, window_len)
    return pinch_block(b, window_ramp)


def auto_correlate_halfs(s):
    slen = len(s)
    return np.sum(np.conjugate(s[0:slen/2]) * s[slen/2:])


def auto_correlation_sync(s, K, cp_len):
    plen = K * 2
    slen = len(s)
    ac = np.zeros(slen - plen, dtype=np.complex)
    nc = np.zeros(slen - plen, dtype=float)

    le_catched = 0
    for i in range(slen - plen):
        # calc auto-correlation
        c = s[i:i+plen]
        #normalize value
        p = calculate_signal_energy(c)
        if p < 1e-4:
            p = 1.0
            le_catched += 1
        ac[i] = auto_correlate_halfs(c)
        nc[i] = 2 * np.abs(ac[i]) / p
    print 'low energy detected:', le_catched
    ic = np.zeros(slen - plen + cp_len, dtype=float)
    for i in range(slen - plen):
        n = i + cp_len
        p_int = nc[n-cp_len:n + 1]
        ic[n] = (1. / len(p_int)) * np.sum(p_int)
    nm = np.argmax(ic)
    cfo = np.angle(ac[nm]) / np.pi
    print 'max corr val', ic[nm], 'with energy:', calculate_average_signal_energy(s[nm: nm + plen])
    return nm, cfo, ic


def cross_correlation_paper_def(s, preamble):
    # naive implementation of the algorithm described in [0]
    # cc = np.zeros(len(s) - len(preamble), dtype=np.complex)
    # s = np.conjugate(s)
    # for i in range(len(s) - len(preamble)):
    #     cc[i] = np.sum(s[i:i + len(preamble)] * preamble)

    # scipy.signal.correlate way of doing things!
    cc = signal.correlate(s, preamble, 'valid')
    # normalize correlation
    cc /= len(preamble)
    return cc


def correct_frequency_offset(signal, cfo):
    n = np.arange(len(signal))
    return signal * np.exp(-2j * np.pi * cfo * n)


def find_frame_start(s, preamble, K, cp_len):
    # avg_signal_ampl = np.sqrt(calculate_average_signal_energy(s))
    # s /= avg_signal_ampl
    if not len(preamble) == 2 * K:
        raise ValueError('Preamble length must be equal to 2K!')
    if calculate_average_signal_energy(preamble) < .99:
        print 'This preamble is not suited'

    # normalize preamble, maybe it helps
    avg_preamble_ampl = np.sqrt(calculate_average_signal_energy(preamble))
    preamble /= avg_preamble_ampl
    print 'average preamble amplitude:', avg_preamble_ampl, 'normalized preamble:', calculate_average_signal_energy(preamble)

    # first part of the algorithm. rough STO sync and CFO estimation!
    nm, cfo, auto_corr_vals = auto_correlation_sync(s, K, cp_len)
    print 'auto correlation nm:', nm, ', cfo:', cfo, ', val:', auto_corr_vals[nm]

    # Fix CFO!
    s = correct_frequency_offset(s, cfo / (2. * K))

    # Now search for exact preamble position!
    pcc = cross_correlation_paper_def(s, preamble)

    # Find peak and suppress minor peaks at +-N/2
    apcc = np.abs(pcc)

    napcc = apcc * auto_corr_vals[0:len(apcc)]
    nc = np.argmax(napcc)
    print 'cross correlation nc:', nc, np.abs(napcc[nc])

    # norm vectors for plotting
    plt.plot(np.abs(apcc) * (np.abs(napcc[nc] / np.abs(apcc[nc]))))
    plt.plot(np.abs(napcc))
    plt.plot(auto_corr_vals)
    peak = auto_corr_vals > 0.5
    plt.plot(peak)
    print 'peak width:', np.sum(peak)

    plt.show()


def sync_test():
    samp_rate = 10e6  # an assumption to work with
    alpha = .5
    M = 33
    K = 32
    block_len = M * K
    L = 2
    cp_len = K
    ramp_len = cp_len / 2

    test_cfo = -.2
    snr_dB = 0.0

    print M, '*', K, '=', block_len, '(', block_len + cp_len, ')'
    print 'assumed samp_rate:', samp_rate, ' with sc_bw:', samp_rate / K
    print 'relative subcarrier CFO:', test_cfo

    data = get_random_qpsk(block_len)
    x = get_gfdm_frame(data, alpha, M, K, L, cp_len, ramp_len)

    preamble, x_preamble = generate_sync_symbol(get_random_qpsk(K), 'rrc', alpha, K, L, cp_len, ramp_len)
    frame = np.concatenate((preamble, x))
    avg_frame_ampl = np.sqrt(calculate_average_signal_energy(frame))
    # frame /= avg_frame_ampl
    print 'frame duration: ', 1e6 * len(frame) / samp_rate, 'us'

    # simulate Noise and frequency offset!
    # frame = correct_frequency_offset(frame, test_cfo / (-2. * K))
    noise_variance = calculate_awgn_noise_variance(frame, snr_dB)
    # noise_variance = 0.0
    s = get_complex_noise_vector(2 * block_len + len(frame), noise_variance)
    s[block_len:block_len + len(frame)] += frame

    print 'frame_start:', block_len, ', short preamble_start:', block_len + cp_len

    find_frame_start(s, x_preamble, K, cp_len)


def preamble_auto_corr_test():
    K = 32
    pn_seq = get_random_qpsk(K)
    pn_symbols = np.tile(pn_seq, 2)
    D = get_data_matrix(pn_symbols, K, True)
    print np.shape(D)
    print np.all(D[0] == D[1])

    pl, p = generate_sync_symbol(pn_seq, 'rrc', .5, K, 2, K, K / 2)
    print np.shape(p)
    acc = auto_correlate_halfs(p)
    print acc, np.angle(acc)

    taps = gfdm_filter_taps('rrc', .5, 2, K, 1)
    A = gfdm_modulation_matrix(taps, 2, K)
    x = A.dot(pn_symbols)
    print np.shape(x)
    acc = auto_correlate_halfs(x)
    print acc, np.angle(acc)


def main():
    np.set_printoptions(precision=4)
    sync_test()
    # preamble_auto_corr_test()
    # n_sc = 9
    # pn = get_random_qpsk(n_sc // 2)
    # print pn
    # s = np.zeros(n_sc, dtype=np.complex)
    # s[1::2] = pn
    # print s
    # x = np.fft.ifft(s, 12)
    # print x[0:6]
    # print x[6:]
    # print auto_correlate_halfs(x)



if __name__ == '__main__':
    main()
