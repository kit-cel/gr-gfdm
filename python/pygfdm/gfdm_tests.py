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
import matplotlib.pyplot as plt
from gfdm_modulation import gfdm_gr_modulator
from receiver import gfdm_rx_fft2
from utils import get_random_samples, get_random_qpsk, get_zero_f_data

def test_transceiver_00():
    K = 32
    M = 8
    overlap = 2
    alpha = 0.5
    oversampling_factor = 1

    tests = 100
    for t in xrange(tests):
        d = get_random_qpsk(M*K)
        tx = gfdm_gr_modulator(d, 'rrc', alpha, M, K, overlap)
        rx = gfdm_rx_fft2(tx,'rrc',alpha,M,K,overlap,oversampling_factor,4,16)
        print(np.max(np.abs(d-rx)))
        if not (np.max(np.abs(d-rx)) < 1e-2):
            raise RuntimeError('Input and Output deviate')

def test_transceiver_01():
    K = 32
    M = 8
    overlap = 2
    alpha = 0.5
    oversampling_factor = 1

    tests = 100
    for t in xrange(tests):
        d = get_random_qpsk(M*K)
        tx = gfdm_gr_modulator(d, 'rrc', alpha, M, K, overlap,compat_mode=False)
        rx = gfdm_rx_fft2(tx,'rrc',alpha,M,K,overlap,oversampling_factor,4,16)
        print(np.max(np.abs(d-rx)))
        if not (np.max(np.abs(d-rx)) < 1e-2):
            raise RuntimeError('Input and Output deviate')



