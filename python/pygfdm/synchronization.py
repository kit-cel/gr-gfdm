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
import matplotlib.pyplot as plt
from modulation import gfdm_modulation_matrix
from filters import get_frequency_domain_filter, gfdm_filter_taps
from gfdm_modulation import gfdm_modulate_block, gfdm_modulate_fft
from cyclic_prefix import add_cyclic_prefix, pinch_block, get_raised_cosine_ramp, get_window_len, get_root_raised_cosine_ramp
from mapping import get_data_matrix, map_to_waveform_resources
from utils import get_random_qpsk, get_complex_noise_vector, calculate_awgn_noise_variance, calculate_average_signal_energy, calculate_signal_energy, magnitude_squared
from utils import generate_seed
from correlation import auto_correlate_halfs, cross_correlate_signal, cross_correlate_fft_cyclic
from preamble import generate_sync_symbol


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


def get_gfdm_frame(data, alpha, M, K, L, cp_len, ramp_len):
    window_len = get_window_len(cp_len, M, K)
    b = gfdm_modulate_fft(data, alpha, M, K, L)
    b = add_cyclic_prefix(b, cp_len)
    window_ramp = get_root_raised_cosine_ramp(ramp_len, window_len)
    return pinch_block(b, window_ramp)


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
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    # print 'max corr val', ic[nm], 'with energy:', calculate_average_signal_energy(s[nm: nm + plen])
    return nm, cfo, ic, ac


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


def correct_frequency_offset(signal, cfo, fft_len=1.):
    phase_inc = cfo_to_phase_increment(cfo, fft_len)
    n = np.arange(len(signal))
    return signal * np.exp(1j * phase_inc * n)


def freq_to_cfo(freq, fft_len, samp_rate):
    # freq and samp_rate should be intuitive!
    # fft_len takes into account how many samples are used for a correlation in order to calculate a CFO.
    # OR: calculate the CFO relative to subcarrier bandwidth! Which is important in this case here!
    sc_bw = samp_rate / float(fft_len)
    return freq / sc_bw


def phase_increment(freq, samp_rate):
    return 2. * np.pi * freq / samp_rate


def cfo_to_phase_increment(cfo, fft_len):
    return (2. * np.pi) * cfo / float(fft_len)


def get_complex_sine(freq, samp_rate, n_samps):
    phase_inc = phase_increment(freq, samp_rate)
    return complex_sine(phase_inc, n_samps)


def complex_sine(phase_inc, n_samps):
    if phase_inc < np.finfo(type(phase_inc)).eps:
        v = np.zeros(n_samps)
    else:
        v = np.arange(0, n_samps * phase_inc, phase_inc)
    return np.cos(v) + 1j * np.sin(v)


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

    # Fix CFO!
    s = correct_frequency_offset(s, cfo / (2. * K))

    # find exact peak over fine STO estimation
    # ONLY nc required, everything else only used for evaluation!
    nc, napcc, apcc = improved_cross_correlation_peak(s, preamble, abs_corr_vals)
    # ALGORITHM FINISHED!
    print 'find_frame_start nc: {} (nm: {}), cfo: {:.5f}, abs_corr_val: {:.5f}'.format(nc, nm, cfo, abs_corr_vals[nm])

    return nc, cfo, abs_corr_vals, corr_vals, napcc, apcc


def multiply_valid(first, second):
    if len(first) >= len(second):
        return first[:-(len(first) - len(second))] * second
    else:
        return first * second[:-(len(second) - len(first))]


def simplified_sync_algo(rx, x_preamble, subcarriers, cp_len):
    x_preamble = initialize_sync_algorithm(x_preamble, subcarriers)
    oac = auto_correlate_signal(rx, subcarriers)
    # this 2 divisor is up to debate. Seems necessary for larger cp_len relative to fft_len
    ac = np.roll(oac, cp_len // 2)

    nm = np.argmax(np.abs(ac))
    cfo = np.angle(ac[nm]) / (2. * np.pi)

    # s = correct_frequency_offset(rx, cfo / (2. * subcarriers))
    s = rx
    xc = np.correlate(s, x_preamble, 'valid')
    cc = multiply_valid(np.abs(ac), np.abs(xc))
    nc = np.argmax(np.abs(cc))

    plt.plot(np.abs(ac[nm - subcarriers:nm + subcarriers]))
    sc = s[nm - subcarriers:nm + subcarriers]
    print(len(sc))
    xcc = cross_correlate_fft_cyclic(sc, x_preamble)
    print(len(xcc))
    plt.plot(np.abs(xcc) * 0.001)
    ccc = np.abs(ac[nm - subcarriers:nm + subcarriers]) * np.abs(xcc)
    plt.plot(ccc * 0.001)
    ncc = np.argmax(np.abs(ccc))
    print 'cyclic: ', ncc, 'ncc: ', nm + ncc - subcarriers, 'cp_len: ', cp_len
    plt.show()

    print 'simplified: nc: {} (nm: {}), cfo: {:.5f}, abs_corr_val: {:.5f}'.format(nc, nm, cfo, np.abs(ac[nm]))
    return nc, cfo, cc


def generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo, init_phase=0.0):
    block_len = M * K
    data = get_random_qpsk(block_len, seed=generate_seed('awesomepayloadblabla'))
    x = get_gfdm_frame(data, alpha, M, K, L, cp_len, ramp_len)

    preamble, x_preamble = generate_sync_symbol(get_random_qpsk(K, seed=generate_seed('awesome')), 'rrc', alpha, K, L, cp_len, ramp_len)
    print 'frame energy:', calculate_average_signal_energy(x), 'preamble energy:', calculate_average_signal_energy(preamble)
    preamble *= np.sqrt(calculate_average_signal_energy(x) / calculate_average_signal_energy(preamble))

    frame = np.concatenate((preamble, x))

    # simulate Noise and frequency offset!
    wave = complex_sine(cfo_to_phase_increment(test_cfo, K), len(frame))
    phase_shift = np.repeat(np.exp(1j * init_phase), len(frame))
    wave *= phase_shift
    frame *= wave
    noise_variance = calculate_awgn_noise_variance(frame, snr_dB)
    s = get_complex_noise_vector(2 * block_len + len(frame), noise_variance)
    s[block_len:block_len + len(frame)] += frame
    return s, x_preamble


def sync_test():
    samp_rate = 10e6  # an assumption to work with
    alpha = .3
    M = 33
    K = 32
    block_len = M * K
    L = 2
    cp_len = K
    ramp_len = cp_len / 2

    test_cfo = -.2
    snr_dB = 40.0
    false_alarm_probability = 1e-3
    print 'Channel parameters, SNR:', snr_dB, 'dB with a relative subcarrier offset:', test_cfo
    print 'assumed samp_rate:', samp_rate, ' with sc_bw:', samp_rate / K

    frame, x_preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo, .5)

    print np.correlate(x_preamble, x_preamble), auto_correlate_halfs(x_preamble)
    print 'frame duration: ', 1e6 * len(frame) / samp_rate, 'us'
    print M, '*', K, '=', block_len, '(', block_len + cp_len, ')'
    print 'frame start in test vector: ', block_len + cp_len
    s = frame

    s *= 100. / np.sqrt(len(x_preamble))
    nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(s, x_preamble, K, cp_len)
    print 'FOUND FRAMESTART nc:', nc, np.abs(napcc[nc]), abs_corr_vals[nc]
    snc, scfo, scc = simplified_sync_algo(s, x_preamble, K, cp_len)
    # print 'signal_len:', len(s), ', auto_corr_len:', len(auto_corr_vals), ', cross_corr_len:', len(napcc), len(s) - len(napcc)
    thr = calculate_threshold_factor(false_alarm_probability) * np.sum(apcc[nc - K:nc + K]) / (2 * K)
    print 'threshold: ', thr
    plt.plot(abs_corr_vals, label='auto corr')
    plt.plot(apcc, label='cross corr')# * (np.abs(napcc[nc] / np.abs(apcc[nc]))))
    plt.plot(napcc, label='combined')
    threshold_peak = np.zeros(len(napcc), dtype=float)
    threshold_peak[nc] = napcc[nc]
    plt.plot((threshold_peak > thr) * (napcc[nc] / thr))
    print 'threshold exceeded points: ', np.sum(threshold_peak > thr)

    for thr in (.3, .4, .5, .6):
        peak = abs_corr_vals > thr
        plt.plot(peak * thr)
        print 'n exceed thr: ', thr, ':', np.sum(peak)

    plt.xlim((1000, 1200))
    plt.legend()
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
    # preamble_auto_corr_test()
    sync_test()
    return
    # cfo = 1.024e-5
    samp_rate = 12.5e6
    freq = 20.
    fft_len = 256
    sc_bw = samp_rate / fft_len
    cfo = freq_to_cfo(freq, fft_len, samp_rate)
    ph_i = cfo_to_phase_increment(cfo, fft_len)
    phase_inc = phase_increment(freq, samp_rate)
    print 'samp_rate: {}, frequency: {}, fft_len: {}'.format(samp_rate, freq, fft_len)
    print 'subcarrier bandwidth: {}, cfo: {}, phase increment: {}/{}'.format(sc_bw, cfo, ph_i, phase_inc)


    # wave = get_complex_sine(freq, samp_rate, 129)
    # s = np.ones(129)
    # s = correct_frequency_offset(s, cfo, fft_len)
    # # print wave
    # print np.abs(wave - s)

    preamble, x_preamble = generate_sync_symbol(get_random_qpsk(fft_len, seed=generate_seed('awesome')), 'rrc', .5, fft_len, 2, fft_len // 2, fft_len // 8)

    init_phase = 0.8
    wave = complex_sine(cfo_to_phase_increment(cfo, fft_len), len(preamble))
    # phase_shift = np.repeat(, len(preamble))
    wave *= np.exp(1j * init_phase)
    preamble *= wave

    print phase_inc

    ac = auto_correlate_halfs(preamble[64:64 + 2 * fft_len])
    fp = np.angle(ac)
    ac_phase_inc = fp / fft_len
    print 'auto corr phase: {}, AC phase {}, phase_inc: {}'.format(phase_inc, fp, ac_phase_inc)

    xc = cross_correlate_signal(preamble, x_preamble)
    ap = np.angle(xc[64])
    # print ap, (ap - init_phase) / (cfo * 2 * np.pi)

    residual = ap - init_phase
    res_phase = residual / (2 * fft_len)
    print residual, res_phase, phase_inc / res_phase
    # print fp, fp / (2 * np.pi)
    # print (ap - init_phase) - fp

    plt.plot(np.abs(xc) / 10000.)
    plt.plot(np.angle(xc))
    plt.show()


if __name__ == '__main__':
    main()
