#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Andrej Rode.
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
from modulation import gfdm_modulation_matrix
from mapping import get_data_stream
from filters import get_frequency_domain_filter, gfdm_freq_taps, gfdm_filter_taps, gfdm_freq_taps_sparse
from gfdm_modulation import get_random_GFDM_block
from utils import get_random_qpsk, calculate_average_signal_energy, map_qpsk_stream
import matplotlib.pyplot as plt

'''
[0] Low Complexity GFDM Receiver Based On Sparse Frequency Domain Processing
'''

def gfdm_transform_input_to_fd(R):
    '''
    :param R: Received symbols in time domain
    :return: Received symbols in frequency domain (Apply N-length FFT)
    '''
    return np.fft.fft(R.astype(np.complex))


def gfdm_extract_subcarriers(R, K, M, L):
    '''
    :param R: Received symbols in frequency domain (DC on 0th bin)
    :param K: number of subcarriers
    :param M: number of symbols per subcarrier
    :param L: overlapping factor
    :return: extracted subcarrier data (length M*L) row-wise in D

    '''
    D = np.empty((K, M * L), np.complex)
    for k in xrange(K):
        for l in xrange(L):
            start_pos = ((k + l + K - 1) % K) * M
            copy_pos = ((l + L / 2) % L) * M
            D[k][copy_pos:copy_pos + M] = R[start_pos:start_pos + M]
    return D


def gfdm_filter_subcarriers(R, H, K, M, L):
    '''
    :param R: data matrix with row-wise subcarrier data
    :param H: filter taps in frequency domain (length M*L) (DC on 0th bin)
    :param K: number of subcarriers
    :param M: number of symbols per subcarrier
    :param L: overlapping factor
    :return: filtered subcarriers in fd
    '''
    H = H / float(K)
    D = R.flatten()
    F = D * np.tile(H, K)
    return np.reshape(F, (-1, L * M)).T


def gfdm_superposition_subcarriers(R, K, M, L):
    '''
    :param R: filtered subcarriers in fd
    :param K: number of subcarriers
    :param M: number of symbols per subcarrier
    :param L: overlapping factor
    :return: L times superpositioned/decimated subcarriers

    '''
    S = np.zeros((K, M), np.complex)
    D = R.T
    for k in xrange(K):
        S[k] = np.sum(np.reshape(D[k], (L, -1)), axis=0)
    return S.T

def gfdm_transform_subcarriers_to_tdomain(R, K, M, L):
    S = np.zeros((K,M), np.complex)
    D = R.T
    for k in xrange(K):
        S[k] = np.fft.ifft(D[k])
    return S.T

def gfdm_get_ic_f_taps(f_taps, M):
    return np.multiply(f_taps[0:M], f_taps[-M:])


def gfdm_map_subcarriers(R, K, M, L):
    D = map_qpsk_stream(R)
    return np.reshape(D, (-1, K))

def gfdm_remove_sc_interference(R, D, K, M, L, H_sic):
    D = D.T
    R = R.T
    R_new = np.empty((K, M), dtype=np.complex)
    for k in xrange(K):
        R_new[k] = R[k] - H_sic * np.fft.fft(D[(k-1) % K] + D[(k+1) % K])
    return R_new.T


def gfdm_demodulate_block(R, H, K, M, L):
    D_0 = gfdm_transform_input_to_fd(R)
    D_1 = gfdm_extract_subcarriers(D_0, K, M, L)
    D_2 = gfdm_filter_subcarriers(D_1, H, K, M, L)
    D_3 = gfdm_superposition_subcarriers(D_2, K, M, L)
    D_4 = gfdm_transform_subcarriers_to_tdomain(D_3, K, M, L)
    print(D_4)
    print(D_4.shape)
    return get_data_stream(D_4)


def gfdm_demodulate_block_sic(R, H, K, M, L, J=0):
    H_sic = gfdm_get_ic_f_taps(H/float(K), M)
    D_0 = gfdm_transform_input_to_fd(R)
    D_1 = gfdm_extract_subcarriers(D_0, K, M, L)
    D_2 = gfdm_filter_subcarriers(D_1, H, K, M, L)
    D_3 = gfdm_superposition_subcarriers(D_2, K, M, L)
    D_4 = gfdm_transform_subcarriers_to_tdomain(D_3, K, M, L)
    for j in xrange(J):
        D_5 = gfdm_map_subcarriers(D_4, K, M, L)
        D_6 = gfdm_remove_sc_interference(D_3, D_5, K, M, L, H_sic)
        D_4 = gfdm_transform_subcarriers_to_tdomain(D_6, K, M, L)
    return get_data_stream(D_4)


def gfdm_gr_receiver(data, filtertype, alpha, M, K, overlap, compat_mode=True):
    H = get_frequency_domain_filter(filtertype, alpha, M, K, overlap) / float(K)
    return gfdm_demodulate_block(data, H.conj(), K, M, overlap)


def gfdm_demodulate_fft(data, alpha, M, K, overlap, sic_rounds=0):
    H = get_frequency_domain_filter('rrc', alpha, M, K, overlap)
    return gfdm_demodulate_block_sic(data, H.conj(), K, M, overlap)


def get_repetition_matrix(timeslots, overlap):
    '''
    [0] has a good representation of the matrices.
    '''
    im = np.identity(timeslots)
    r = np.tile(im, overlap)
    r = r.T
    return r


def get_permutation_matrix(timeslots, subcarriers, overlap):
    im = np.identity(overlap * timeslots / 2)
    zm = np.zeros((overlap * timeslots / 2, (subcarriers - 1) * overlap * timeslots / 2))
    p0 = np.hstack((im, zm))
    p1 = np.hstack((zm, im))
    p = np.vstack((p0, p1))
    return p.T


def get_nth_sc_permutation_matrix(sc, timeslots, subcarriers, overlap):
    p0 = get_permutation_matrix(timeslots, subcarriers, overlap)
    pn = np.roll(p0, sc * timeslots, axis=0)
    return pn


def gfdm_rx_to_fd(rx):
    return np.fft.fft(rx)


def gfdm_rx_sc_bins(rx, sc_num, timeslots, overlap):
    sc = np.roll(rx, int(overlap * timeslots // 2 - timeslots * sc_num))
    return np.fft.fftshift(sc[:overlap * timeslots])


def gfdm_downsample_fft(rx, overlap):
    # this would be equivalent! np.fft.ifft(rx)[::overlap]
    return np.fft.ifft(np.sum(np.reshape(rx, (overlap, -1)), axis=0)) / overlap


def gfdm_demodulate_fft_loop(rx, timeslots, subcarriers, overlap, sparse_freq_taps):
    frx = gfdm_rx_to_fd(rx)
    rxdata = np.zeros((subcarriers, timeslots), dtype=rx.dtype)
    for sc in range(subcarriers):
        sc_syms = gfdm_rx_sc_bins(frx, sc, timeslots, overlap)
        filtered = sc_syms * sparse_freq_taps / (1. * subcarriers)
        downsampled = gfdm_downsample_fft(filtered, overlap)
        downsampled *= overlap  # this is weird, but GR compat please!
        rxdata[sc, :] = downsampled
    return rxdata.flatten()


def main():
    '''
    This is a comparison for 3 different demodulation approaches.
    matched filter matrix being the 'benchmark'
    The other two should converge towards the matrix approach for overlap -> subcarriers
    Actually, there's a bug in the 'GR' approach, thus it only works for overlap==2
    '''
    timeslots = 25
    subcarriers = 16
    overlap = 2
    time_taps = gfdm_filter_taps('rrc', .5, timeslots, subcarriers, 1)
    freq_taps = gfdm_freq_taps(time_taps)
    sparse_freq_taps = gfdm_freq_taps_sparse(freq_taps, timeslots, overlap)
    A = gfdm_modulation_matrix(time_taps, timeslots, subcarriers, 1, True)
    Ainv = np.linalg.inv(A)
    Amf = np.conjugate(A).T

    tx_syms = get_random_qpsk(timeslots * subcarriers)
    rx_syms = A.dot(tx_syms)

    mf_matrix_rx = Amf.dot(rx_syms)
    inv_matrix_rx = Ainv.dot(rx_syms)
    gr_res = gfdm_demodulate_block(rx_syms, sparse_freq_taps, subcarriers, timeslots, overlap)
    fft_res = gfdm_demodulate_fft_loop(rx_syms, timeslots, subcarriers, overlap, sparse_freq_taps)

    mf_matrix_rx *= np.sqrt(calculate_average_signal_energy(fft_res) / calculate_average_signal_energy(mf_matrix_rx))
    inv_matrix_rx *= np.sqrt(calculate_average_signal_energy(fft_res) / calculate_average_signal_energy(inv_matrix_rx))
    gr_res *= np.sqrt(calculate_average_signal_energy(fft_res) / calculate_average_signal_energy(gr_res))

    print 'compare demodulation accuracy for different approaches'
    for e in range(11):
        em = 10 ** (-1. * e)
        matrixvsloop = np.all(np.abs(fft_res - mf_matrix_rx) < em)
        grvsmatrix = np.all(np.abs(gr_res - mf_matrix_rx) < em)
        grvsloop = np.all(np.abs(gr_res - fft_res) < em)
        print 'error margin {:.1e}\tMFmatriXvsGR: {}\tMFmatriXvsLoop: {}\tGRvsLoop: {}'.format(em, grvsmatrix, matrixvsloop, grvsloop)


if __name__ == '__main__':
    main()
