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

import numpy as np
import commpy as cp
import scipy.signal as signal
import matplotlib.pyplot as plt
from gfdm_plot_utils import plot_gfdm_matrix
from modulation import gfdm_filter_taps, gfdm_modulation_matrix, gfdm_tx_fft2, reshape_input, gfdm_freq_taps, gfdm_freq_taps_sparse


def get_data_matrix(data, K):
    return np.reshape(data, (-1, K))


def gfdm_gr_modulator(x, filtertype, alpha, M, K, L):
    h = gfdm_filter_taps(filtertype, alpha, M, K, 1)
    H = gfdm_freq_taps(h)
    H_sparse = gfdm_freq_taps_sparse(H, M, L)

    # Sort Input subcarrierwise
    D = np.reshape(x, (-1, K))

    W = np.fft.fft(D.astype(np.complex), M, axis=0)

    R = np.reshape(np.tile(W.flatten(), L), (-1, K))

    F = R.T.flatten() * np.tile(H_sparse, K)
    Fd = np.reshape(F, (-1, L * M)).T

    #this seems important!
    Fs = np.fft.fftshift(Fd, axes=0)

    x_out = np.zeros((M * K) + (L - 1) * M, dtype='complex')
    # print np.shape(x_out)

    # FIXME: Correct overlap handling is key!
    # Add data-vector to correct position -max neg frequency : 0 :
    # max_pos_frequency
    for k in range(K):
        x_out[k * M:(k + L) * M] += Fs[:, k]
    # Add 'oversampled' parts of first subcarrier to end and 'oversampled' parts
    # of last subcarrier to start
    tail_length = (L - 1) * M / 2
    x_first = x_out[0:tail_length]
    x_last = x_out[-tail_length:]

    x_out = x_out[tail_length:-tail_length]
    # print np.shape(x_out)
    x_out[0:tail_length] += x_last
    x_out[-tail_length:] += x_first

    x_s = np.fft.ifftshift(x_out)

    x_t = np.fft.ifft(x_s)
    x_t *= 1.0 / K
    return x_t


# [0] Generalized frequency division multiplexing: Analysis of an alternative multi-carrier technique for next generation cellular systems
def gfdm_modulate_fft(data, alpha, M, K, overlap, oversampling):
    # this function aims to reproduce [0] Section IIIA
    print M, K, overlap, oversampling
    print data

    # compare [0]: D holds the data symbols with dk being the kth column holding M data symbols.
    D = get_data_matrix(data, K)
    print D
    # plt.plot(D[:, 0])
    # compare [0]: M-point FFT on each vector dk, aka each column of D, use periodicity of circular convolution.
    W = np.fft.fft(D.astype(np.complex), M, axis=0) #* (1 / np.sqrt(M))
    # plt.plot(abs(W[:, 0]))

    # compare [0]: upsampling with with R(L). aka repeat each dk L-times. L is called overlap
    R = np.reshape(np.tile(W.flatten(), overlap), (-1, K))
    # plt.plot(abs(R[:, 0]))


    # compare [0]: multiply with filter in freq domain. Denoted as capital gamma in paper.
    taps = gfdm_filter_taps('rrc', alpha, M, overlap, oversampling)
    Gamma = np.fft.fft(taps, overlap * M) #* (1 / np.sqrt(overlap * M))
    # plt.plot(abs(Gamma))
    plt.plot(taps)
    plt.show()

    F = R.T.flatten() * np.tile(Gamma, K)
    Fd = np.reshape(F, (-1, overlap * M)).T
    # Fd = np.fft.fftshift(Fd, axes=0)

    # compare [0]: shift samples to match desired frequency
    X = np.zeros(M * K, dtype=np.complex)
    for k in range(K):
        s = np.zeros(M * K, dtype=np.complex)
        p = Fd[:, k]
        s[0:len(p)] = p
        s = np.roll(s, k * M)
        X += s
        # plt.plot(s)

    # plt.plot(abs(Fd[:, 0]))
    #
    # plt.show()
    # X = np.fft.fftshift(X)
    x = np.fft.ifft(X)
    return x, Gamma


def my_compare_function():
    M = 32
    K = 8
    alpha = .5
    oversampling_factor = 1
    overlap = 2

    tests = 100
    for t in range(tests):
        d = np.random.standard_normal(2 * M * K)
        d = np.reshape(d, (2, -1))
        d = d[0] + 1j * d[1]

        xo = gfdm_tx_fft2(d, 'rrc', alpha, M, K, overlap, oversampling_factor)
        xn = gfdm_gr_modulator(d, 'rrc', alpha, M, K, overlap)

        if not np.all(xo == xn):
            raise RuntimeError('Function results deviate')



def main():
    my_compare_function()
    M = 8
    K = 4
    alpha = .5
    oversampling_factor = 1
    overlap = 4


    # taps = gfdm_filter_taps('rrc', alpha, M, K, oversampling_factor)
    # A0 = gfdm_modulation_matrix(taps, M, K, oversampling_factor, rearrange_indices=False)
    # print 'GFDM shape: ', np.shape(A0), 'tap length:', len(taps)
    # # plot_gfdm_matrix(A0)
    #
    #
    # fig = plt.figure()
    # # data = signal.gaussian(K, 1.0)
    # data = np.arange(0, K / 4, 1. / K) + 1
    # data = np.tile(data, M)
    # print data
    # xA = A0.dot(data)
    #
    # xA *= (1. / K)
    # x1 = gfdm_tx_fft2(data, 'rrc', alpha, M, K, overlap, oversampling_factor)
    #
    # x0, Gamma = gfdm_modulate_fft(data, alpha, M, K, overlap, oversampling_factor)
    # x0 *= (1. / K)
    #
    #
    # plt.plot(np.real(x0), 'b')
    # plt.plot(np.real(x1), 'g')
    # plt.plot(np.real(xA), 'r')
    # plt.plot(np.imag(x0), 'b--')
    # plt.plot(np.imag(x1), 'g--')
    # plt.plot(np.imag(xA), 'r--')
    # # plt.show()
    # #
    # # H = np.arange(M * K)
    # # print H
    # # L = overlap
    # # print (H[0:(M * L) / 2], H[-(M * L) / 2:])
    # # plt.plot(Gamma)
    # # plt.plot(H)
    # # plt.plot(H_sparse)
    # plt.show()


if __name__ == '__main__':
    main()