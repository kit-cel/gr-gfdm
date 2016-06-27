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

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import gfdm_swig as gfdm
import numpy as np
from pygfdm.utils import calculate_signal_energy
from pygfdm.synchronization import auto_correlate_halfs, generate_test_sync_samples, find_frame_start
import matplotlib.pyplot as plt

class qa_sync_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        alpha = .5
        M = 33
        K = 32
        block_len = M * K
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.1
        snr_dB = 15.0
        print "SNR:", snr_dB, "[dB], CFO:", test_cfo

        pre_gen = gfdm.preamble_generator(K, alpha, K * 2)
        preamble = np.array(pre_gen.get_preamble())

        signal = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        nc, cfo, auto_corr_vals, napcc, apcc = find_frame_start(signal, preamble, K, cp_len)

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)
        iqres = kernel.find_preamble(signal)
        res = np.zeros(len(iqres), dtype=float)
        # iqres = np.array(iqres[0:len(iqres) // 2])
        # res[0::2] = iqres.real
        # res[1::2] = iqres.imag
        # print np.shape(res)
        p = iqres[1088:1088+2*K]
        print p
        val = auto_correlate_halfs(p)
        print 'CFO residual:', np.angle(val) / np.pi, ', PYTHON residual: ', test_cfo - cfo

if __name__ == '__main__':
    # gr_unittest.run(qa_sync_cc, "qa_sync_cc.xml")
    gr_unittest.run(qa_sync_cc)
