#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Johannes Demel.
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
A few hints on used papers, consider them to be a good read.
[0] Generalized frequency division multiplexing: Analysis of an alternative multi-carrier technique for next generation cellular systems
[1] Generalized Frequency Division Multiplexing for 5th Generation Cellular Networks
[2] "Bit Error Rate Performance of Generalized Frequency Division Multiplexing"

'''

import numpy as np
import commpy as cp
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from gfdm_plot_utils import plot_gfdm_matrix
from modulation import gfdm_modulation_matrix, gfdm_tx_fft2, reshape_input, get_data_matrix
from filters import gfdm_filter_taps, gfdm_freq_taps, gfdm_freq_taps_sparse, get_frequency_domain_filter


def gfdm_transform_subcarriers_to_fd(D, M):
    '''
    compare [0] Section IIIA: M-point FFT on each vector dk, aka each column of D, use periodicity of circular convolution.
    :param D: data matrix, each column of D represents the symbols located on one subcarrier.
    :param M: number of symbols on a subcarrier and FFT size.
    :return: data on subcarriers in frequency domain representation. DC on 0th bin.
    '''
    F = np.fft.fft(D.astype(np.complex), M, axis=0)
    return F


def gfdm_upsample_subcarriers_in_fd(D, K, L):
    '''
    compare [0] Section IIIA: upsampling with with R(L). aka repeat each dk L-times. L is called overlap
    :param D: matrix with FD representation of data.
    :param K: number of subcarriers
    :param L: overlapping factor. aka upsampling integer value.
    :return: upsampled data in frequency domain.
    '''
    F = np.reshape(np.tile(D.flatten(), L), (-1, K))
    # F = np.fft.fftshift(F, axes=0)
    return F


def gfdm_filter_subcarriers_in_fd(D, H, M, K, L):
    '''
    compare [0] Section IIIA: multiply with filter in freq domain. Denoted as capital gamma in paper.
    :param D: matrix with upsampled data in FD
    :param H: transfer function of filter. DC on 0-th bin.
    :param M: number of time slots
    :param K: number of subcarriers
    :param L: overlap factor
    :return: filtered subcarriers in FD
    '''
    # H = np.fft.fftshift(H)
    F = D.T.flatten() * np.tile(H, K)
    return np.reshape(F, (-1, L * M)).T


def gfdm_subcarrier_modulator_in_fd(D, H, M, K, L):
    # function sums up portion of FFT-based modulation which is exactly as in 'gfdm_tx_fft2'.
    W = gfdm_transform_subcarriers_to_fd(D, M)
    R = gfdm_upsample_subcarriers_in_fd(W, K, L)
    F = gfdm_filter_subcarriers_in_fd(R, H, M, K, L)
    return F


def gfdm_combine_subcarriers_in_fd(F, M, K, L, compat_mode=True):
    '''
    compare [0] Section IIIA: shift samples to match desired frequency
    :param F: matrix with k-th subcarrier in k-th column in FD.
    :param M: number of time slots
    :param K: number of subcarriers
    :param L: overlap factor
    :return: GFDM vector block in frequency domain.
    '''
    # FFT-shift necessary here.
    F = np.fft.fftshift(F, axes=0)
    tail_length = (L - 1) * M
    X = np.zeros(M * K + tail_length, dtype=np.complex)
    for k in range(K):
        X[k * M: k * M + L * M] += F[:, k]
    X[0: tail_length] += X[-tail_length:]
    X = X[0:-tail_length]

    # This last step confuses things! put one on DC!
    if compat_mode:
        X = np.roll(X, -M / 2)
    else:
        X = np.roll(X, -M * L / 2)
    return X


def gfdm_modulate_block(D, H, M, K, L, compat_mode=True):
    '''
    this function aims to reproduce [0] Section IIIA
    It is a combination of all steps necessary for FFT based modulation.
    :param D: data matrix. symbols for subcarrier on columns. different grouping in compat_mode
    :param H: filter transfer function
    :param M: number of time slots
    :param K: number of subcarriers
    :param L: overlap factor
    :param compat_mode: define if function modulates in accordance with 'gfdm_tx_fft2' for tests.
    :return:
    '''
    F = gfdm_subcarrier_modulator_in_fd(D, H, M, K, L)
    x_out = gfdm_combine_subcarriers_in_fd(F, M, K, L, compat_mode=compat_mode)
    x_t = x_out

    if compat_mode:
        x_t = np.fft.ifftshift(x_out)

    x_t = np.fft.ifft(x_t)

    if compat_mode:
        x_t *= 1.0 / K
    return x_t


def gfdm_gr_modulator(x, filtertype, alpha, M, K, L, compat_mode=True):
    # this function aims to reproduce [0] Section IIIA
    H = get_frequency_domain_filter(filtertype, alpha, M, K, L)

    # compare [0]: D holds symbols grouped by subcarrier in each column.
    D = get_data_matrix(x, K, group_by_subcarrier=compat_mode)

    return gfdm_modulate_block(D, H, M, K, L, compat_mode=compat_mode)


def gfdm_modulate_fft(data, alpha, M, K, overlap):
    # this function aims to reproduce [0] Section IIIA

    H = get_frequency_domain_filter('rrc', alpha, M, K, overlap)
    D = get_data_matrix(data, K, group_by_subcarrier=False)
    return gfdm_modulate_block(D, H, M, K, overlap, False)
    return x


def get_random_samples(nsamples):
    d = np.random.standard_normal(2 * nsamples)
    d = np.reshape(d, (2, -1))
    return d[0] + 1j * d[1]

def gr_conformity_validation():
    M = 32
    K = 8
    alpha = .5
    oversampling_factor = 1
    overlap = 2

    tests = 100
    for t in range(tests):
        # d = np.random.standard_normal(2 * M * K)
        # d = np.reshape(d, (2, -1))
        # d = d[0] + 1j * d[1]
        d = get_random_samples(M * K)

        xo = gfdm_tx_fft2(d, 'rrc', alpha, M, K, overlap, oversampling_factor)
        xn = gfdm_gr_modulator(d, 'rrc', alpha, M, K, overlap)

        if not np.all(xo == xn):
            raise RuntimeError('Function results deviate')


def get_zero_f_data(k, K, M):
    data = np.zeros(K)
    data[k] = 1.
    # data = np.tile(data, M)
    data = np.repeat(data, M)
    return data


def validate_subcarrier_location(alpha, M, K, overlap, oversampling_factor):
    goofy_ordering = False
    taps = gfdm_filter_taps('rrc', alpha, M, K, oversampling_factor)
    A0 = gfdm_modulation_matrix(taps, M, K, oversampling_factor, group_by_subcarrier=goofy_ordering)

    n = np.arange(M * K * oversampling_factor, dtype=np.complex)
    for k in range(K):
        f = np.exp(1j * 2 * np.pi * (float(k) / (K * oversampling_factor)) * n)
        F = abs(np.fft.fft(f))
        fm = 1. * np.argmax(F) / M

        data = get_zero_f_data(k, K, M)

        x0 = gfdm_gr_modulator(data, 'rrc', alpha, M, K, overlap, compat_mode=goofy_ordering) * (2. / K)
        f0 = 1. * np.argmax(abs(np.fft.fft(x0))) / M

        xA = A0.dot(get_data_matrix(data, K, group_by_subcarrier=goofy_ordering).flatten()) * (1. / K)
        fA = 1. * np.argmax(abs(np.fft.fft(xA))) / M
        if not fm == fA == f0:
            raise RuntimeError('ERROR: subcarriers are not located at the same bins!')


def compare_subcarrier_location(alpha, M, K, overlap, oversampling_factor):
    goofy_ordering = False
    taps = gfdm_filter_taps('rrc', alpha, M, K, oversampling_factor)
    A0 = gfdm_modulation_matrix(taps, M, K, oversampling_factor, group_by_subcarrier=goofy_ordering)
    n = np.arange(M * K * oversampling_factor, dtype=np.complex)
    colors = iter(cm.rainbow(np.linspace(0, 1, K)))

    for k in range(K):
        color = next(colors)
        f = np.exp(1j * 2 * np.pi * (float(k) / (K * oversampling_factor)) * n)
        F = abs(np.fft.fft(f))
        fm = np.argmax(F) / M
        plt.plot(F, '-.', label=k, color=color)

        data = get_zero_f_data(k, K, M)

        x0 = gfdm_gr_modulator(data, 'rrc', alpha, M, K, overlap, compat_mode=goofy_ordering) * (2. / K)
        f0 = 1. * np.argmax(abs(np.fft.fft(x0))) / M
        plt.plot(abs(np.fft.fft(x0)), label='FFT' + str(k), color=color)

        xA = A0.dot(get_data_matrix(data, K, group_by_subcarrier=goofy_ordering).flatten()) * (1. / K)
        fA = np.argmax(abs(np.fft.fft(xA))) / M
        plt.plot(abs(np.fft.fft(xA)), '-', label='matrix' + str(k), color=color)
        print fm, fA, f0


def compare_subcarrier_combination():
    M = 8
    K = 4
    overlap = 2
    data = get_zero_f_data(0, K, M)
    data += 2 * get_zero_f_data(3, K, M)
    print data
    F = get_data_matrix(data, K, False)
    print F
    F = np.reshape(np.tile(F.flatten(), overlap), (-1, K))

    X0 = np.real(gfdm_combine_subcarriers_in_fd(F, M, K, overlap))
    print X0
    X = np.zeros(M * K, dtype=np.complex)
    for k in range(K):
        s = np.zeros(M * K, dtype=np.complex)
        s[0:M * overlap] = F[:, k]
        s = np.roll(s, k * M - M / 2)
        X += s

    print np.real(np.roll(X, -M / 2))
    print np.all(np.roll(X, -M / 2) == X0)
    print np.all(X == X0)
    print np.real(X)


def main():
    gr_conformity_validation()
    compare_subcarrier_combination()
    M = 8
    K = 16
    alpha = .5
    oversampling_factor = 1
    overlap = 2
    validate_subcarrier_location(alpha, M, K, overlap, oversampling_factor)
    compare_subcarrier_location(alpha, M, K, overlap, oversampling_factor)

    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()