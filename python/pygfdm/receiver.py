#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Andrej, Rode.
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
from filters import gfdm_filter_taps
from mapping import reshape_input
from modulation import transmitMatrix


def gfdm_rx(y, filtertype, alpha, M, K, L, N, QAM, J):
    '''
    y: Transmit-Symbols (length M*K*N)
    filtertype: ['rrc','rc']
    alpha: rolloff-factor
    rx_strat: ['zf','mf']
    M: number of timeslots
    K: numbor of subcarrier
    N: oversampling-factor
    '''
    A = transmitMatrix(filtertype, alpha, M, K, N)
    A_rx = A.conj().transpose()
    rx = A_rx.dot(y)
    return rx


def gfdm_rx_fft2(y, filtertype, alpha, M, K, L, N, QAM, J):
    '''
    y: transmitted gfdm-block (length: M*K samples)
    filtertype: ('rrc'|'rc')
    alpha: (0,1) float
    M: number of slots
    K: number of subcarriers
    L: freq-domain length of filter
    Low-complexity receiver implementation as proposed by G.Fettweis
    (based on sparse frequency Domain Processing)

    output: demodulated samples in original order (first K samples in timeslot 1, second K ...)
    '''
    h = gfdm_filter_taps(filtertype, alpha, M, K, 1)
    h = np.roll(h, h.shape[-1] / 2)
    H_rx = np.fft.fft(h)
    H_sparse = np.concatenate((H_rx[0:M * L / 2], H_rx[-M * L / 2:]))
    y_ifft = np.array([])
    y = (1.0 / K) * y
    # Transfer input to frequency domain and center around 0th frequency bin
    y_f = np.fft.fftshift(np.fft.fft(y))
    # Filter and superposition in frequency domain
    Y_fs = gfdm_rx_filsup(y_f, H_sparse, M, K, L)
    # Demodulate per subcarrier
    y_ifft = gfdm_rx_demod(Y_fs, K)
    if J > 0:
        y_ifft = gfdm_rx_sic(K, M, J, H_sparse, y_ifft, Y_fs, QAM)
        y_ifft = np.reshape(y_ifft, (K * M))
    # Sort output in timeslot,subcarrier order
    y_ifft = reshape_input(y_ifft, K, M)
    return y_ifft


def gfdm_rx_demod(Y_fs, K):
    '''
    Y_fs: received samples filtered and superpositioned in frequency domain (not centered) KxM-array
    K: Number of subcarriers

    output: demodulated samples in subcarrier order (first M samples are on subcarrier 1, second M....)
    '''
    y_ifft = np.array([])
    for k in xrange(K):
        y_ifft = np.concatenate((y_ifft, np.fft.ifft(Y_fs[k])))
    return y_ifft


def gfdm_rx_filsup(y_f, H_sparse, M, K, L):
    '''
    y_f: input samples centered in frequency domain 1xK*M-array
    H_sparse: Rx-filter per subcarrier - length (M*L)
    M: number of time slots
    K: number of subcarrier
    L: width of sparse Rx-filter in number of subcarrier

    output: (K,M) - array
    '''
    y_out = np.empty((K, M), dtype='complex')
    y_f = np.concatenate((y_f[-(L - 1) * M / 2:], y_f, y_f[0:(L - 1) * M / 2]))
    for k in xrange(K):
        # select kth subcarrier
        y_down = y_f[k * M:(k + L) * M]
        # 'uncenter' in Frequency domain
        y_down = np.fft.ifftshift(y_down)
        # apply filter in frequency domain (not centered)
        y_filter = np.multiply(y_down, H_sparse)
        # Superposition L samples in frequency domain
        y_out[k] = np.sum(y_filter.reshape(L, M), axis=0)
    return y_out


def gfdm_rx_sic(K, M, J, H_sparse, d_rx, Y_fs, QAM):
    '''
    K: Number of subcarriers
    M: Number of slots
    J: Number of Iterations for Interference Cancellation
    H_sparse: Rx-filter of length M*L in frequency domain
    d_rx: mapped symbols before Interference cancellation (sorted by subcarrier)
    Y_fs: filtered, superpositioned input samples in frequency domain (not centered) KxM-array
    QAM: QAM order
    '''
    import commpy as cp
    # Receive all subcarriers in F-Domain
    # map each symbol to closest QAM - Point
    # d_rx s
    qam_mod = cp.QAMModem(QAM)
    # Calculate rising/falling flank interference coefficients
    H_rf = np.multiply((H_sparse[0:M] / K), (H_sparse[M:] / K))
    # Reshape mapped symbols into per-subcarrier array
    d_p = np.empty((K, M), dtype='complex')
    # d_p (K,M)
    for k in xrange(K):
        d_p[k] = qam_mod.mapping(d_rx[k * M:(k + 1) * M], 'hard')
    for j in xrange(J):
        y = np.empty((K, M), dtype='complex')
        for k in xrange(K):
            y[k] = Y_fs[k] - H_rf * np.fft.fft(d_p[(k - 1) % K] + d_p[(k + 1) % K])
            # Recalculate d_rx
        d_rx = gfdm_rx_demod(y, K)
        for k in xrange(K):
            d_p[k] = d_rx[k * M:(k + 1) * M]
            d_p[k] = qam_mod.mapping(d_p[k], 'hard')
    return d_rx