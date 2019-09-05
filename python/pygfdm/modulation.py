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
from .filters import gfdm_filter_taps, gfdm_freq_taps, gfdm_freq_taps_sparse
from .mapping import reshape_input


def gfdm_modulation_matrix(filter_taps, M, K, oversampling_factor=1, group_by_subcarrier=False):
    '''
    This function returns a GFDM modulation matrix
    :param filter_taps: M*K*N length filter tap array
    :param M: number of time slots
    :param K: number of subcarriers
    :param oversampling_factor: factor for oversampling
    :param group_by_subcarrier: if True group by time symbol not subcarrier
    In order to comply with vector-element-sorting used for the efficient impl, this parameter should be 'True'!
    :return: modulation matrix
    [0] Generalized Frequency Division Multiplexing for 5th Generation Cellular Networks
    [1] Generalized frequency division multiplexing: Analysis of an alternative multi-carrier technique for next generation cellular systems

    CAREFUL: It is tricky for a number of reasons.
    First off, only in [1] oversampling is considered. Mostly it's undesirable.
    Secondly, the definitions differ slightly. In [0] frequency modulation is defined by -1j * 2 * np.pi,
    whereas in [1] it is 1j * 2 * np.pi. Notice the sign change!
    '''
    N = M * K

    filter_taps = np.roll(filter_taps, (N * oversampling_factor) // 2)
    A = np.zeros((N * oversampling_factor, N), dtype=np.complex)

    n = np.arange(N * oversampling_factor, dtype=np.complex)
    for m in range(M):
        for k in range(K):
            f_mod = np.exp(1j * 2 * np.pi * (float(k) / (K * oversampling_factor)) * n)
            g = filter_taps * f_mod
            g = np.roll(g, m * K * oversampling_factor)
            A[:, m * K + k] = g

    if group_by_subcarrier:
        indices = np.arange(M * K)
        indices = np.reshape(indices, (-1, K)).T.flatten()
        A = A[:, indices]
    return A


def transmitMatrix(filtertype, alpha, M, K, oversampling_factor=1):
    # replaces old definition.
    taps = gfdm_filter_taps(filtertype, alpha, M, K, oversampling_factor)
    return gfdm_modulation_matrix(taps, M, K, oversampling_factor, False)


def gfdm_tx(x, filtertype, alpha, M, K, L, N):
    '''
    x: Input-Symbols (length M*K)
    filtertype: ['rrc','rc']
    alpha: rolloff-factor
    M: number of timeslots
    K: number of subcarrier
    oversampling_factor: sometimes referred to as N
    '''
    A = transmitMatrix(filtertype, alpha, M, K, N)
    A *= M
    tx = A.dot(x)
    return tx


def gfdm_tx_fft2(x, filtertype, alpha, M, K, L, N):
    '''
    x: Input-Array (length: M*K symbols)
    filtertype: ('rrc'|'rc')
    alpha: (0,1) float
    M: number of slots
    K: number of subcarriers
    L: freq-domain length of filter

    Low-complexity transmitter implementation as proposed by G. Fettweis
    '''
    h = gfdm_filter_taps(filtertype, alpha, M, K, 1)
    H = gfdm_freq_taps(h)
    H_sparse = gfdm_freq_taps_sparse(H, M, L)

    # Sort Input subcarrierwise
    x = reshape_input(x, M, K)
    x_out = np.zeros((M * K) + (L - 1) * M, dtype='complex')
    for k in xrange(K):
        # M rows and L columns with respective FFT output
        # pick symbols per subcarrier
        x_k = x[k * M:((k + 1) * M)]
        # perform fft and switch to frequency domain
        x_f = np.fft.fft(x_k)
        # copy values of M-point DFT to obtain MK-point DFT
        x_f_L = np.tile(x_f, L)
        # make data-vector 'sparse'
        # x_f_L = np.concatenate((x_f_K[0:(M*L)/2], x_f_K[-(M*L)/2:]))
        # filter with sparse filter taps in frequency domain
        x_fil = np.multiply(x_f_L, H_sparse)
        # Add data-vector to correct position -max neg frequency : 0 :
        # max_pos_frequency
        x_out[k * M:(k + L) * M] = x_out[k * M:(L + k) * M] + np.fft.fftshift(x_fil)
    # Add 'oversampled' parts of first subcarrier to end and 'oversampled' parts
    # of last subcarrier to start
    x_first = x_out[0:(L - 1) * M / 2]
    x_last = x_out[-(L - 1) * M / 2:]
    x_out = x_out[(L - 1) * M / 2:-(L - 1) * M / 2]
    x_out[0:(L - 1) * M / 2] = x_out[0:(L - 1) * M / 2] + x_last
    x_out[-(L - 1) * M / 2:] = x_out[-(L - 1) * M / 2:] + x_first

    x_t = np.fft.ifft(np.fft.ifftshift(x_out))
    x_t *= 1.0 / K
    return x_t


def check_matrix_invertibility(A, Ainv, err_margin=1e-10):
    '''
    Use this function to get an idea if the chosen parameters (timeslots, subcarriers) result in an invertible matrix.
    This is a 'hands-on' approach at estimating if a given matrix is invertible.
    The idea is that given an invertible matrix, if multiplied with its inverse, the result should be an almost perfect identity matrix.
    Imperfections are caused be numerical errors.
    Maximum deviation from the identity matrix is used as a measure.
    '''
    Icheck = A.dot(Ainv)
    msize = np.shape(Icheck)[0]
    Iexact = np.identity(msize, dtype=A.dtype)
    Idev = np.abs(Icheck - Iexact)
    if not np.all(Idev < err_margin):
        m = np.max(np.max(Idev))
        cond_number = np.linalg.cond(A)
        print('Inversion failed, maximum deviation', m, ' exceeds error margin', err_margin)
        print('Matrix condition number {:.2f} results in approx. {:.2f}bit precision loss'.format(cond_number, np.log10(cond_number)))
        print(np.abs(Icheck[0:10, 0:10]))
        return False
    return True


def main():
    A = transmitMatrix('rrc', .5, 9, 16)
    Ainv = np.linalg.inv(A)
    check_matrix_invertibility(A, Ainv)


if __name__ == '__main__':
    main()
