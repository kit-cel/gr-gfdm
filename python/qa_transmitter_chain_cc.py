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
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.gfdm_modulation import gfdm_modulate_block
from pygfdm.mapping import get_data_matrix
from pygfdm.utils import get_random_qpsk
from pygfdm.cyclic_prefix import get_window_len, get_raised_cosine_ramp, add_cyclic_prefix, pinch_block
from pygfdm.synchronization import get_sync_symbol
import numpy as np


class qa_transmitter_chain_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        np.set_printoptions(precision=2)
        n_frames = 3
        alpha = .5
        M = 8
        K = 4
        L = 2
        cp_len = 8
        ramp_len = 4
        block_len = M * K
        window_len = get_window_len(cp_len, M, K)
        taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
        window_taps = get_raised_cosine_ramp(ramp_len, window_len)
        pn_symbols = get_random_qpsk(K)
        H_preamble = get_frequency_domain_filter('rrc', alpha, 2, K, L)
        preamble = get_sync_symbol(pn_symbols, H_preamble, K, L, cp_len, ramp_len)[0]

        ref = np.array([], dtype=np.complex)
        data = np.array([], dtype=np.complex)
        for i in range(n_frames):
            d = get_random_qpsk(block_len)
            D = get_data_matrix(d, K, group_by_subcarrier=False)
            b = gfdm_modulate_block(D, taps, M, K, L, False)
            b = add_cyclic_prefix(b, cp_len)
            b = pinch_block(b, window_taps)
            ref = np.concatenate((ref, preamble, b))
            data = np.concatenate((data, d))

        src = blocks.vector_source_c(data)
        mod = gfdm.simple_modulator_cc(M, K, L, taps)
        prefixer = gfdm.cyclic_prefixer_cc(cp_len, ramp_len, block_len, window_taps)
        preambler = blocks.vector_insert_c(preamble, window_len + len(preamble), 0)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, prefixer, preambler, dst)
        self.tb.run()
        res = np.array(dst.data())[0:len(ref)]

        self.assertComplexTuplesAlmostEqual(ref, res, 5)

    # def test_002_big_data(self):
    #     print "big data test"
    #     reps = 5
    #     alpha = .5
    #     M = 127
    #     K = 16
    #     L = 4
    #     taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
    #     data = np.array([], dtype=np.complex)
    #     ref = np.array([], dtype=np.complex)
    #     for i in range(reps):
    #         d = get_random_qpsk(M * K)
    #         D = get_data_matrix(d, K, group_by_subcarrier=False)
    #         ref = np.append(ref, gfdm_modulate_block(D, taps, M, K, L, False))
    #         data = np.append(data, d)
    #     # print data
    #     # print ref
    #     # print "MAXIMUM ref value: ", np.max(abs(ref))
    #
    #     src = blocks.vector_source_c(data)
    #     mod = gfdm.simple_modulator_cc(M, K, L, taps)
    #     dst = blocks.vector_sink_c()
    #
    #     self.tb.connect(src, mod, dst)
    #     # set up fg
    #     self.tb.run()
    #     # check data
    #     res = np.array(dst.data())
    #     # res /= M * K
    #     # print "MAXIMUM result value: ", np.max(abs(res))
    #
    #     self.assertComplexTuplesAlmostEqual(ref, res, 2)


if __name__ == '__main__':
    # gr_unittest.run(qa_simple_modulator_cc, "qa_simple_modulator_cc.xml")
    gr_unittest.run(qa_transmitter_chain_cc)
