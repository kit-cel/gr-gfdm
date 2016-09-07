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
from filters import get_frequency_domain_filter
from gfdm_modulation import get_random_GFDM_block
from utils import map_qpsk_stream
import matplotlib.pyplot as plt


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


def gfdm_demodulate_block_sic(R, H, K, M, L, J=1):
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


def gfdm_demodulate_fft(data, alpha, M, K, overlap, sic_rounds=1):
    H = get_frequency_domain_filter('rrc', alpha, M, K, overlap)
    return gfdm_demodulate_block_sic(data, H.conj(), K, M, overlap,sic_rounds)


def main():
    K = 32
    M = 15
    overlap = 2
    alpha = .5
    (data, block) = get_random_GFDM_block(M, K, overlap, alpha)
    rx_data = gfdm_demodulate_fft(block*(1./K), alpha, M, K , overlap, 10)
    plt.plot(rx_data)
    plt.plot(data)
    axes = plt.gca()
    axes.set_xlim([0,20])
    plt.show()

if __name__ == "__main__":
    main()
