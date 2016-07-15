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
from pygfdm.utils import calculate_signal_energy, get_complex_noise_vector
from pygfdm.synchronization import auto_correlate_halfs, generate_test_sync_samples, find_frame_start
from pygfdm.synchronization import auto_correlate_signal, abs_integrate, correct_frequency_offset
from pygfdm.synchronization import cross_correlate_naive, cross_correlate_signal


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

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)
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
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)

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
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        ac = auto_correlate_signal(signal, K)
        ic = abs_integrate(np.abs(ac), cp_len)
        window_size = 1024
        kic = np.array([], dtype=np.complex)
        for i in range(0, len(signal), window_size):
            w = ac[i:i + window_size]
            kic = np.concatenate((kic, kernel.abs_integrate_preamble(w)))
        # kic = np.array(kernel.abs_integrate_preamble(ac))
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
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)

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
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)

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

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)
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
        snr_dB = 20.0

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        # print 'preamble size:', len(preamble), len(preamble) == 2 * K
        # print 'init kernel'
        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)
        preamble = np.array(kernel.preamble())
        nc, cfo, auto_corr_vals, corr_vals, napcc, apcc = find_frame_start(signal, preamble, K, cp_len)

        # print 'kernel.find_preamble'
        knc = np.array(kernel.find_preamble(signal))
        print knc, nc
        self.assertEqual(nc, knc)

    def test_008_auto_correlation_stepped(self):
        print 'test window buffering for auto correlation'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.2
        snr_dB = 10.0

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        signal *= .001

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)
        preamble = np.array(kernel.preamble())
        nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(signal, preamble, K, cp_len)

        step_size = 1088
        for i in range(0, len(signal), step_size):
            w = signal[i:i + step_size + 3 * K]
            a = abs_corr_vals[i:i + step_size]
            ref_buf = abs_corr_vals[i + step_size - 2 * K:i + step_size]
            cref_buf = corr_vals[i + step_size - 2 * K:i + step_size]
            inref_buf = signal[i + step_size - 2 * K:i + step_size]
            if len(w) > 4 * K and len(ref_buf) == 2 * K:
                bi_buf = np.array(kernel.integration_buffer())
                kbuf = np.concatenate((bi_buf, a))
                print 'KBEF argmax', np.argmax(kbuf)
                snc = int(kernel.find_preamble(w))  # necessary to update kernel state!
                ai_buf = np.array(kernel.integration_buffer())
                self.assertFloatTuplesAlmostEqual(ref_buf, ai_buf, 5)
                c_buf = np.array(kernel.auto_corr_buffer())
                self.assertComplexTuplesAlmostEqual(cref_buf, c_buf, 5)
                in_buf = np.array(kernel.input_buffer())
                self.assertComplexTuplesAlmostEqual(inref_buf, in_buf)

    def test_009_step_window(self):
        print 'long windowed test'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2

        test_cfo = -.2
        snr_dB = 10.0

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        signal *= .001

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)
        preamble = np.array(kernel.preamble())
        nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(signal, preamble, K, cp_len)
        slen = len(signal)
        n_rep = 10
        signal = np.tile(signal, n_rep)

        nc_vec = (np.arange(0, n_rep) * slen) + nc

        step_size = nc + 4
        window_nc = np.array([], dtype=int)
        dumped = np.array([], dtype=int)
        for i in range(0, len(signal), step_size):
            w = signal[i:i + step_size + 3 * K]
            if len(w) > 4 * K:
                # it = i // step_size
                # print 'PROCESS WINDOW ', it
                snc = int(kernel.find_preamble(w))
                if not snc == -2 * len(w):
                    abs_nc = i + snc
                    n_frame = len(window_nc)
                    prev_pos = nc_vec[n_frame - 1]
                    if abs_nc == prev_pos:
                        print 'found:', abs_nc, 'previous:', prev_pos
                    elif not abs_nc == nc_vec[n_frame]:
                        print 'FAIL STATUS: reps', n_rep, 'detected frames:', n_frame, 'avail_frames:', len(nc_vec)
                        print '#frame', n_frame, 'expected:', nc_vec[n_frame], 'found:', abs_nc, 'diff:', abs_nc - nc_vec[n_frame]
                        print 'step_size:', step_size, 'step_start:', i, 'snc:', snc
                        print nc_vec[n_frame-3:n_frame+3]
                        print window_nc[-3:]
                        self.assertEqual(abs_nc, nc_vec[n_frame])
                    else:
                        # print 'SUCCESSFULLY detected frame @', abs_nc
                        window_nc = np.append(window_nc, abs_nc)
                else:
                    dumped = np.append(dumped, snc)

        self.assertTupleEqual(tuple(nc_vec), tuple(window_nc))

    def test_010_sync_block(self):
        print 'GR block sync!'
        alpha = .5
        M = 33
        K = 32
        L = 2
        cp_len = K
        ramp_len = cp_len / 2
        frame_len = 2 * K + cp_len + M * K

        test_cfo = -.2
        snr_dB = 10.0

        signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
        signal *= .001

        kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble, 4000)
        nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(signal, np.array(kernel.preamble()), K, cp_len)
        print 'Python frame start: ', nc

        src = blocks.vector_source_c(signal)
        sync = gfdm.sync_cc(K, cp_len, frame_len, preamble, 'gfdm_block')
        snk = blocks.vector_sink_c()
        self.tb.connect(src, sync, snk)
        self.tb.run()

        ref = signal[nc:nc + frame_len]
        res = np.array(snk.data())
        print 'res length:', len(res)
        print 'frame size:', frame_len
        self.assertComplexTuplesAlmostEqual(ref[32:], res[32:])

    # def test_011_block_noise_input(self):
    #     print '\n\n\ntest_011_Noise only!'
    #     alpha = .5
    #     M = 33
    #     K = 32
    #     L = 2
    #     cp_len = K
    #     ramp_len = cp_len / 2
    #     frame_len = 2 * K + cp_len + M * K
    #
    #     test_cfo = -.2
    #     snr_dB = 10.0
    #
    #     signal, preamble = generate_test_sync_samples(M, K, L, alpha, cp_len, ramp_len, snr_dB, test_cfo)
    #     # signal *= .001
    #     #
    #     # kernel = gfdm.improved_sync_algorithm_kernel_cc(K, cp_len, preamble)
    #     # nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(signal, np.array(kernel.preamble()), K, cp_len)
    #
    #     noise_variance = .5
    #     signal = get_complex_noise_vector(4 * M * K, noise_variance)
    #
    #     src = blocks.vector_source_c(signal)
    #     sync = gfdm.sync_cc(K, cp_len, frame_len, preamble, 'gfdm_block')
    #     snk = blocks.vector_sink_c()
    #     self.tb.connect(src, sync, snk)
    #     self.tb.run()
    #
    #     # ref = signal[nc:nc + frame_len]
    #     res = np.array(snk.data())
    #     print 'res length:', len(res)
    #     print 'frame size:', frame_len
    #     # self.assertComplexTuplesAlmostEqual(ref, res)


if __name__ == '__main__':
    # gr_unittest.run(qa_sync_cc, "qa_sync_cc.xml")
    gr_unittest.run(qa_sync_cc)
