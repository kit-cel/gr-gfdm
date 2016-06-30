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
import scipy as sp
from pygfdm.utils import calculate_signal_energy
from pygfdm.synchronization import auto_correlate_halfs, generate_test_sync_samples, find_frame_start
from pygfdm.synchronization import auto_correlate_signal, abs_integrate, correct_frequency_offset
from pygfdm.synchronization import cross_correlate_naive, cross_correlate_signal
import matplotlib.pyplot as plt

class qa_sync_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_norm_preamble(self):
        alpha = .5
        K = 32
        cp_len = K

        pre_gen = gfdm.preamble_generator(K, alpha, K * 2)
        preamble = np.array(pre_gen.get_preamble())

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)
        kernel_preamble = np.array(kernel.preamble())

        # normalize preamble
        preamble /= np.sqrt(np.abs(calculate_signal_energy(preamble)))
        print 'norm_preamble energy:', calculate_signal_energy(preamble)

        self.assertComplexTuplesAlmostEqual(kernel_preamble, preamble)

    def test_002_auto_correlate(self):
        print 'auto correlation test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.1
        snr_dB = 15.0

        pre_gen = gfdm.preamble_generator(K, alpha, K * 2)
        preamble = np.array(pre_gen.get_preamble())
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        ac = auto_correlate_signal(signal, K)
        kac = np.array(kernel.auto_correlate_preamble(signal))
        self.assertComplexTuplesAlmostEqual(ac, kac, 6)

    def test_003_abs_integrate(self):
        print 'abs_integrate test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.1
        snr_dB = 15.0

        pre_gen = gfdm.preamble_generator(K, alpha, K * 2)
        preamble = np.array(pre_gen.get_preamble())
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        ac = auto_correlate_signal(signal, K)
        ic = abs_integrate(np.abs(ac), cp_len)
        kic = np.array(kernel.abs_integrate_preamble(ac))
        self.assertFloatTuplesAlmostEqual(ic[cp_len:], kic[cp_len:], 6)

    def test_004_find_auto_correlation_peak(self):
        print 'find auto correlation peak test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.1
        snr_dB = 15.0

        pre_gen = gfdm.preamble_generator(K, alpha, K * 2)
        preamble = np.array(pre_gen.get_preamble())
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        ac = auto_correlate_signal(signal, K)
        ic = abs_integrate(np.abs(ac), cp_len)

        nm = np.argmax(ic)
        knm = kernel.find_peak_preamble(ic)
        self.assertEqual(nm, knm)

        cfo = np.angle(ac[nm]) / np.pi
        kcfo = kernel.calculate_normalized_cfo_preamble(ac[knm])
        self.assertAlmostEqual(cfo, kcfo)

    def test_005_remove_cfo(self):
        print 'remove CFO test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.1
        snr_dB = 15.0

        pre_gen = gfdm.preamble_generator(K, alpha, K * 2)
        preamble = np.array(pre_gen.get_preamble())
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        ac = auto_correlate_signal(signal, K)
        ic = abs_integrate(np.abs(ac), cp_len)
        nm = np.argmax(ic)
        cfo = np.angle(ac[nm]) / np.pi

        s = correct_frequency_offset(signal, cfo / (2. * K))
        ks = kernel.remove_cfo_preamble(signal, cfo)
        self.assertComplexTuplesAlmostEqual(s, ks, 3)  # VOLK rotator is inaccurate. Thus, accuracy == 3

    def test_006_cross_correlation(self):
        print 'cross correlation test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = 0.0
        snr_dB = 15.0
        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)
        preamble = np.array(kernel.preamble())

        kcc = np.array(kernel.cross_correlate_preamble(signal))
        cc = cross_correlate_signal(signal, preamble)
        self.assertComplexTuplesAlmostEqual(kcc, cc, 4)

    def test_007_frame(self):
        print 'find frame test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.2
        snr_dB = 0.0

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)
        preamble = np.array(kernel.preamble())
        nc, cfo, auto_corr_vals, corr_vals, napcc, apcc = find_frame_start(signal, preamble, K, cp_len)

        knc = np.array(kernel.find_preamble(signal))
        print knc, nc
        self.assertEqual(nc, knc)


if __name__ == '__main__':
    # gr_unittest.run(qa_sync_cc, "qa_sync_cc.xml")
    gr_unittest.run(qa_sync_cc)
