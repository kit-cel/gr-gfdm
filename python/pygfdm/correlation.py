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
Relevant paper section
[0] A synchronization technique for generalized frequency division multiplexing
[1] Improved Preamble-Aided Timing Estimation for OFDM Systems

COMMENT
[1] describes the algorithm while [0] additionally explains GFDM preamble generation.
'''

import numpy as np
import matplotlib.pyplot as plt


def auto_correlate_halfs(s):
    pivot = len(s) / 2
    return np.sum(np.conjugate(s[0:pivot]) * s[pivot:])


def cross_correlate_naive(s, p):
    # naive implementation of the algorithm described in [0]
    p = np.conjugate(p)
    p_len = len(p)
    buf_len = len(s) - p_len + 1
    cc = np.zeros(buf_len, dtype=np.complex)
    for i in range(buf_len):
        cc[i] = np.sum(s[i:i + p_len] * p)
    return cc


def cross_correlate_signal_full(s, p):
    return np.correlate(s, p, 'full')
    # return signal.correlate(s, p, 'valid')[0:len(s)-len(p)]


def cross_correlate_signal_valid(s, p):
    # scipy.signal.correlate way of doing things!
    return np.correlate(s, p, 'valid')


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


def validate_naive_cross_correlation(s, p, tolerance):
    cs = cross_correlate_signal_valid(s, p)
    cf = cross_correlate_naive(s, p)
    check_results(cs, cf, tolerance, err_msg='Naive cross correlation algorithms produce differing results!')


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
    l = 30
    N = l * 3
    tolerance = 1e-12

    a = np.array([1., 2., 3., 4.])
    b = np.array([3., 2., 0., 2.])

    validate_naive_cross_correlation(a, b, tolerance)
    validate_full_cross_correlation(a, b, tolerance)
    validate_valid_cross_correlation(a, b, tolerance)

    s = np.random.randn(N) + 1j * np.random.randn(N)
    x = s[0:l]
    validate_full_cross_correlation(s, x, tolerance)
    validate_valid_cross_correlation(s, x, tolerance)
    validate_naive_cross_correlation(s, x, tolerance)


def cross_correlate_fft_cyclic(s, p):

    # valid_res_len = len(s) - len(p) + 1
    S = np.fft.fft(s[0: len(p)])
    # p = np.append(p, np.zeros(len(s) - len(p)))
    P = np.fft.fft(p)
    P = np.conjugate(P)
    C = S * P
    cf = np.fft.ifft(C)
    # cf = cf[0:valid_res_len]
    if s.dtype == np.float:  # in case input was float.
        cf = np.real(cf)
    return cf


def main():
    np.set_printoptions(precision=4, suppress=True)
    validate_cross_correlation_algorithms()
    l = 30
    N = l * 3
    offset = 5
    s = np.random.randn(N) + 1j * np.random.randn(N)
    x = s[offset:offset + l]

    cf = cross_correlate_fft_valid(s, x)
    plt.plot(np.abs(cf))
    ccc = cross_correlate_fft_cyclic(s, x)
    plt.plot(np.abs(ccc))
    plt.show()


if __name__ == '__main__':
    main()