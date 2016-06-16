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


def gfdm_filter_taps(filtertype, alpha, M, K, oversampling_factor):
    N = oversampling_factor
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M * K * N, alpha, N * K, 1)
    elif filtertype == "rc":
        time_h, h = cp.rcosfilter(M * K * N, alpha, N * K, 1)
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