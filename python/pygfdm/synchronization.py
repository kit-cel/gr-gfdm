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
from mapping import get_data_matrix, map_to_waveform_resources
from utils import get_random_qpsk, get_complex_noise_vector, calculate_awgn_noise_variance, calculate_average_signal_energy, calculate_signal_energy, magnitude_squared


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


def calculate_packet_average(d, avg_len):
    s = np.concatenate((d, np.zeros(avg_len - (len(d) % avg_len))))
    t = np.reshape(s, (-1, avg_len))
    return np.sum(t, axis=1)


def detect_energy_ramps(d, alpha):
    tr = alpha * np.roll(d, 1)
    thr = d > tr
    return np.arange(len(d))[thr]


def detect_frame_energy(data, alpha=50., avg_len=32):
    d = magnitude_squared(data)
    t = calculate_packet_average(d, avg_len)
    peaks = detect_energy_ramps(t, alpha)[1:]  # just dump first peak, because it's probably a np.roll artefact.
    peak_pos = peaks * avg_len
    return peak_pos


def generate_seed(my_str):
    return abs(hash(my_str))  #% (2 ** 32), seed must be a positive integer.


def mapped_preamble(seed, filtertype, alpha, active_subcarriers, fft_len, subcarrier_map, overlap, cp_len, ramp_len):
    pn_vals = get_random_qpsk(active_subcarriers, seed)
    pn_sym = map_to_waveform_resources(pn_vals, active_subcarriers, fft_len, subcarrier_map)
    return generate_sync_symbol(pn_sym, filtertype, alpha, fft_len, overlap, cp_len, ramp_len)


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
    pivot = len(s) / 2
    return np.sum(np.conjugate(s[0:pivot]) * s[pivot:])


def auto_correlate_signal(signal, K):
    plen = K * 2
    slen = len(signal)
    ac = np.zeros(slen - plen, dtype=np.complex)

    for i in range(slen - plen):
        # calc auto-correlation
        c = signal[i:i+plen]
        #normalize value
        p = calculate_signal_energy(c)
        ac[i] = 2 * auto_correlate_halfs(c) / p
    return ac


def abs_integrate(nc, cp_len):
    ic = np.zeros(len(nc), dtype=float)
    for n in range(cp_len, len(nc)):
        p_int = nc[n-cp_len:n + 1]
        ic[n] = (1. / len(p_int)) * np.sum(p_int)
    return ic


def auto_correlation_sync(s, K, cp_len):
    plen = K * 2

    ac = auto_correlate_signal(s, K)
    nc = np.abs(ac)
    ic = abs_integrate(nc, cp_len)

    nm = np.argmax(ic)
    cfo = np.angle(ac[nm]) / np.pi
    # print 'max corr val', ic[nm], 'with energy:', calculate_average_signal_energy(s[nm: nm + plen])
    return nm, cfo, ic, ac


def cross_correlate_naive(s, p):
    # naive implementation of the algorithm described in [0]
    cc = np.zeros(len(s) - len(p), dtype=np.complex)
    p = np.conjugate(p)
    p_len = len(p)
    buf_len = len(s) - p_len
    for i in range(buf_len):
        cc[i] = np.sum(s[i:i + p_len] * p)
    return cc


def cross_correlate_signal_full(s, p):
    return signal.correlate(s, p, 'full')
    # return signal.correlate(s, p, 'valid')[0:len(s)-len(p)]


def cross_correlate_signal_valid(s, p):
    # scipy.signal.correlate way of doing things!
    return signal.correlate(s, p, 'valid')


def cross_correlate_signal(s, p):
    # scipy.signal.correlate way of doing things!
    return cross_correlate_signal_valid(s, p)[0:len(s) - len(p)]  # drop last sample.


def cross_correlate_fft(s, p):
    C = np.conjugate(np.fft.fft(s)) * np.fft.fft(p, len(s))
    cf = np.fft.ifft(C)[0:len(s) - len(p)]
    return cf


def cross_correlate_fft_full(s, p):
    # function is equivalent to scipy.signal.correlate(s, p, 'full')
    pad_head = len(s) - len(p)
    s_zero_pad = len(s) - 1
    p_zero_pad = s_zero_pad + len(s) - len(p)
    s = np.append(s, np.zeros(s_zero_pad))
    p = np.append(p, np.zeros(p_zero_pad))
    S = np.fft.fft(s)
    P = np.fft.fft(p)
    P = np.conjugate(P)
    C = S * P
    cf = np.fft.ifft(C)
    cf = np.fft.fftshift(cf)
    cf = cf[pad_head:]
    if s.dtype == np.float:  # in case input was float.
        cf = np.real(cf)
    return cf


def cross_correlate_fft_valid(s, p):
    # function is equivalent to scipy.signal.correlate(s, p, 'valid')
    valid_res_len = len(s) - len(p) + 1
    S = np.fft.fft(s)
    p = np.append(p, np.zeros(len(s) - len(p)))
    P = np.fft.fft(p)
    P = np.conjugate(P)
    C = S * P
    cf = np.fft.ifft(C)
    cf = cf[0:valid_res_len]
    if s.dtype == np.float:  # in case input was float.
        cf = np.real(cf)
    return cf


def validate_full_cross_correlation(s, p, tolerance):
    cs = cross_correlate_signal_full(s, p)
    cf = cross_correlate_fft_full(s, p)
    check_results(cs, cf, tolerance, err_msg='Full cross correlation algorithms produce differing results!')


def validate_valid_cross_correlation(s, p, tolerance):
    cs = cross_correlate_signal_valid(s, p)
    cf = cross_correlate_fft_valid(s, p)
    check_results(cs, cf, tolerance, err_msg='Valid cross correlation algorithms produce differing results!')


def check_results(s, p, tolerance, err_msg):
    err = np.abs(s - p) < tolerance
    if not np.all(err):
        print s
        print p
        print err
        raise ValueError('check_results: ' + err_msg)
    else:
        # print 'passed: ' + err_msg
        pass


def validate_cross_correlation_algorithms():
    l = 3
    N = l * 3
    tolerance = 1e-12

    a = np.array([1., 2., 3., 4.])
    b = np.array([3., 2., 0., 2.])
    validate_full_cross_correlation(a, b, tolerance)
    validate_valid_cross_correlation(a, b, tolerance)

    s = np.random.randn(N) + 1j * np.random.randn(N)
    x = s[0:l]
    validate_full_cross_correlation(s, x, tolerance)
    validate_valid_cross_correlation(s, x, tolerance)


def cross_correlation_paper_def(s, preamble):
    cc = cross_correlate_signal(s, preamble)
    # cc = cross_correlate_naive(s, preamble)
    cc /= len(preamble)
    return cc


def improved_cross_correlation_peak(s, preamble, abs_auto_corr_vals):
    apcc = np.abs(cross_correlation_paper_def(s, preamble))
    # Find peak and suppress minor peaks at +-N/2
    min_array_len = min(len(apcc), len(abs_auto_corr_vals))
    napcc = apcc[0:min_array_len] * abs_auto_corr_vals[0:min_array_len]
    nc = np.argmax(napcc)
    # energy = calculate_signal_energy(s[nc:nc + len(preamble)])
    # napcc *= 1. / np.sqrt(energy)
    return nc, napcc, apcc


def correct_frequency_offset(signal, cfo):
    n = np.arange(len(signal))
    return signal * np.exp(-2j * np.pi * cfo * n)


def initialize_sync_algorithm(preamble, K):
    # initialization part
    if not len(preamble) == 2 * K:
        raise ValueError('Preamble length must be equal to 2K!')

    # normalize preamble, maybe it helps
    avg_preamble_ampl = np.sqrt(calculate_average_signal_energy(preamble))
    preamble /= avg_preamble_ampl
    if np.abs(calculate_average_signal_energy(preamble) - 1.0) > 1e-6:
        raise ValueError('preamble not properly normalized!')

    # print 'average preamble amplitude:', avg_preamble_ampl, 'normalized preamble:', calculate_average_signal_energy(preamble)
    return preamble


def calculate_threshold_factor(false_alarm_prob):
    # obviously: false_alarm_prob < 1.0
    if not false_alarm_prob < 1.0:
        raise ValueError('False alarm probability MUST be smaller 1.0!')
    return np.sqrt(-(4 / np.pi) * np.log(false_alarm_prob))



def find_frame_start(s, preamble, K, cp_len):
    # initialization part
    preamble = initialize_sync_algorithm(preamble, K)

    # THIS PARTS represents the live algorithm!
    # first part of the algorithm. rough STO sync and CFO estimation!
    nm, cfo, abs_corr_vals, corr_vals = auto_correlation_sync(s, K, cp_len)
    # print 'auto correlation nm:', nm, ', cfo:', cfo, ', val:', auto_corr_vals[nm]

    # Fix CFO!
    s = correct_frequency_offset(s, cfo / (2. * K))

    # find exact peak over fine STO estimation
    # ONLY nc required, everything else only used for evaluation!
    nc, napcc, apcc = improved_cross_correlation_peak(s, preamble, abs_corr_vals)
    # ALGORITHM FINISHED!

    return nc, cfo, abs_corr_vals, corr_vals, napcc, apcc


def generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo):
    block_len = M * K
    data = get_random_qpsk(block_len)
    x = get_gfdm_frame(data, alpha, M, K, L, cp_len, ramp_len)

    preamble, x_preamble = generate_sync_symbol(get_random_qpsk(K), 'rrc', alpha, K, L, cp_len, ramp_len)
    frame = np.concatenate((preamble, x))

    # simulate Noise and frequency offset!
    frame = correct_frequency_offset(frame, test_cfo / (-2. * K))
    noise_variance = calculate_awgn_noise_variance(frame, snr_dB)
    s = get_complex_noise_vector(2 * block_len + len(frame), noise_variance)
    s[block_len:block_len + len(frame)] += frame
    return s, x_preamble


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
    false_alarm_probability = 1e-3
    print 'Channel parameters, SNR:', snr_dB, 'dB with a relative subcarrier offset:', test_cfo


    print 'assumed samp_rate:', samp_rate, ' with sc_bw:', samp_rate / K

    data = get_random_qpsk(block_len)
    x = get_gfdm_frame(data, alpha, M, K, L, cp_len, ramp_len)

    preamble, x_preamble = generate_sync_symbol(get_random_qpsk(K), 'rrc', alpha, K, L, cp_len, ramp_len)
    frame = np.concatenate((preamble, x))

    print 'frame duration: ', 1e6 * len(frame) / samp_rate, 'us'
    print M, '*', K, '=', block_len, '(', block_len + cp_len, ')'
    # simulate Noise and frequency offset!
    frame = correct_frequency_offset(frame, test_cfo / (-2. * K))
    noise_variance = calculate_awgn_noise_variance(frame, snr_dB)
    # noise_variance = 0.0
    s = get_complex_noise_vector(2 * block_len + len(frame), noise_variance)
    s[block_len:block_len + len(frame)] += frame
    print 'frame start in test vector: ', block_len + cp_len

    s *= 100. / np.sqrt(len(x_preamble))
    nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(s, x_preamble, K, cp_len)
    print 'FOUND FRAMESTART nc:', nc, np.abs(napcc[nc]), abs_corr_vals[nc]
    # print 'signal_len:', len(s), ', auto_corr_len:', len(auto_corr_vals), ', cross_corr_len:', len(napcc), len(s) - len(napcc)
    thr = calculate_threshold_factor(false_alarm_probability) * np.sum(apcc[nc - K:nc + K]) / (2 * K)
    print 'threshold: ', thr
    plt.plot(abs_corr_vals)
    plt.plot(apcc)# * (np.abs(napcc[nc] / np.abs(apcc[nc]))))
    plt.plot(napcc)
    threshold_peak = np.zeros(len(napcc), dtype=float)
    threshold_peak[nc] = napcc[nc]
    plt.plot((threshold_peak > thr) * (napcc[nc] / thr))
    print 'threshold exceeded points: ', np.sum(threshold_peak > thr)

    for thr in (.3, .4, .5, .6):
        peak = abs_corr_vals > thr
        plt.plot(peak * thr)
        print 'n exceed thr: ', thr, ':', np.sum(peak)

    # plt.xlim((1000, 1200))
    plt.show()


def preamble_auto_corr_test():
    K = 32
    pn_seq = get_random_qpsk(K)
    pn_symbols = np.tile(pn_seq, 2)
    D = get_data_matrix(pn_symbols, K, True)
    # print np.shape(D)
    print 'subcarriers bear same symbols:', np.all(D[0] == D[1])

    pl, p = generate_sync_symbol(pn_seq, 'rrc', .5, K, 2, K, K / 2)
    # print np.shape(p)
    acc = auto_correlate_halfs(p)
    print acc, np.angle(acc)

    taps = gfdm_filter_taps('rrc', .5, 2, K, 1)
    A = gfdm_modulation_matrix(taps, 2, K)
    x = A.dot(pn_symbols)
    # print np.shape(x)
    acc = auto_correlate_halfs(x)
    print acc, np.angle(acc)


def main():
    np.set_printoptions(precision=4, suppress=True)
    preamble_auto_corr_test()
    validate_cross_correlation_algorithms()
    sync_test()

    # for i in np.arange(0.7, 8):
    #     fap = 10 ** -i
    #     print fap, calculate_threshold_factor(fap)








if __name__ == '__main__':
    main()
