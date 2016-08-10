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


import numpy as np


def get_random_qpsk(nsamples, seed=None, dtype=np.complex):
    if seed:
        np.random.seed(seed)
    d = np.random.randint(0, 2, 2 * nsamples) * -2. + 1.
    d = np.reshape(d, (2, -1))
    d = d[0] + 1j * d[1]
    return d.astype(dtype=dtype)


def get_random_samples(nsamples, seed=None, dtype=np.complex):
    if seed:
        np.random.seed(seed)
    d = np.random.standard_normal(2 * nsamples)
    d = np.reshape(d, (2, -1))
    d = d[0] + 1j * d[1]
    return d.astype(dtype=dtype)


def randomQAMSymbols(length, M):
    '''
     length: number of symbols to generate
     M: M-QAM - Order (4,16,64,...)
    '''
    n = np.sqrt(M / 4)
    if np.around(n) - n > 1e-10:
        raise Exception('M must be power of 4')
    n = int(n)
    n_M_pos = np.array([1 + 2 * i for i in xrange(n)])
    n_M_neg = np.array([-1 - 2 * i for i in xrange(n)])
    choices = np.concatenate((n_M_pos, n_M_neg))
    return np.array(
        [np.random.choice(choices) + 1j * np.random.choice
        (choices) for i in xrange(length)])


def get_zero_f_data(k, K, M):
    data = np.zeros(K)
    data[k] = 1.
    # data = np.tile(data, M)
    data = np.repeat(data, M)
    return data


def magnitude_squared(input_signal):
    return input_signal.real ** 2 + input_signal.imag ** 2


def calculate_signal_energy(input_signal):
    return np.sum(magnitude_squared(input_signal))


def calculate_average_signal_energy(input_signal):
    return calculate_signal_energy(input_signal) / len(input_signal)


# function adopted from scikit-commpy. Separated noise variance calculation and noise vector generation.
def calculate_awgn_noise_variance(input_signal, snr_dB, rate=1.0):
    avg_energy = calculate_average_signal_energy(input_signal)
    snr_linear = 10. ** (snr_dB / 10.0)
    noise_variance = avg_energy/(2*rate*snr_linear)
    return noise_variance


# corresponds to 'calculate_awgn_noise_variance.
def get_complex_noise_vector(nsamples, noise_variance):
    if noise_variance == 0.0:
        return np.zeros(nsamples, dtype=np.complex)
    return (np.sqrt(noise_variance) * np.random.randn(nsamples)) + (np.sqrt(noise_variance) * np.random.randn(nsamples) * 1j)

