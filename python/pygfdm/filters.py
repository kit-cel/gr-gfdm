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
import matplotlib.pyplot as plt


def gfdm_filter_taps(filtertype, alpha, M, K, oversampling_factor=1.):
    N = oversampling_factor
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M * K * N, alpha, 1. * K * N, 1.)
    elif filtertype == "rc":
        time_h, h = cp.rcosfilter(M * K * N, alpha, 1. * K * N , 1. )
    return h


def gfdm_freq_taps(h):
    h = np.roll(h, h.shape[-1] / 2)
    H = np.fft.fft(h)
    return H


def gfdm_freq_taps_sparse(H, M, L):
    H_sparse = np.concatenate((H[0:(M * L) / 2], H[-(M * L) / 2:]))
    return H_sparse


def get_frequency_domain_filter(filtertype, alpha, M, K, L):
    h = gfdm_filter_taps(filtertype, alpha, M, K, 1)
    H = gfdm_freq_taps(h)
    return gfdm_freq_taps_sparse(H, M, L)


def tapered_cosine(t, alpha):
    d = (1 - 4 * (alpha ** 2) * (t ** 2))
    z = np.where(d == 0.)
    d[z] = 1.
    f = np.cos(np.pi * alpha * t) / d
    f[z] = 0.
    return f


def sinc(t):
    z = np.where(t == 0.)
    it = np.ones(len(t))
    indices = np.delete(np.arange(len(t)), z)
    it[indices] = t[indices]
    f = np.sin(np.pi * it) / (np.pi * it)
    f[z] = 1.
    return f


def freq_tapered_raised_cosine(t, alpha):
    f = sinc(t) * tapered_cosine(t, alpha)
    return f


def main():
    sc = 16
    ts = 8
    overlap = 2
    alpha = .5
    time_taps = gfdm_filter_taps('rc', alpha, ts, sc)
    freq_taps = gfdm_freq_taps(time_taps)
    freq_taps_sparse = gfdm_freq_taps_sparse(freq_taps, ts, overlap)

    fig = plt.figure()
    tp = fig.add_subplot('211')
    t = np.arange(0, ts, 1. / sc)
    plt.plot(t, np.abs(time_taps))
    plt.plot(t, np.abs(freq_tapered_raised_cosine(t - 1. * ts / 2., alpha)))
    plt.xlabel('time samples')

    fp = fig.add_subplot('212')
    f = np.arange(0, sc, 1. / ts)
    plt.plot(f, np.abs(freq_taps))
    plt.plot(f, np.abs(np.fft.fft(freq_tapered_raised_cosine(t - 1. * ts / 2., alpha))))
    plt.plot(f, np.abs(np.concatenate((freq_taps_sparse[0:len(freq_taps_sparse)/2], np.zeros(sc * ts - len(freq_taps_sparse)), freq_taps_sparse[len(freq_taps_sparse)/2:]))), linestyle='--')
    plt.xlabel('frequency samples')
    plt.show()


if __name__ == '__main__':
    main()
