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
from pygfdm.mapping import get_data_matrix, map_to_waveform_resources, get_subcarrier_map
from pygfdm.utils import get_random_qpsk, calculate_signal_energy
from pygfdm.cyclic_prefix import get_window_len, get_raised_cosine_ramp, add_cyclic_prefix, pinch_block, add_cyclic_starfix
from pygfdm.preamble import get_sync_symbol
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
        active = 8
        M = 8
        K = 16
        L = 2
        cp_len = 8
        cs_len = 4
        ramp_len = 4
        block_len = M * K
        window_len = get_window_len(cp_len, M, K, cs_len)
        taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
        taps /= np.sqrt(calculate_signal_energy(taps) / M)
        window_taps = get_raised_cosine_ramp(ramp_len, window_len)
        pn_symbols = get_random_qpsk(K)
        H_preamble = get_frequency_domain_filter('rrc', alpha, 2, K, L)
        preamble = get_sync_symbol(pn_symbols, H_preamble, K, L, cp_len, ramp_len)[0]
        # smap = np.arange(active) + (K - active) // 2
        smap = get_subcarrier_map(K, active, dc_free=True)

        ref = np.array([], dtype=np.complex)
        data = np.array([], dtype=np.complex)
        frame_len = window_len + len(preamble)
        frame_gap = np.zeros(frame_len)
        for i in range(n_frames):
            d = get_random_qpsk(active * M)
            dd = map_to_waveform_resources(d, active, K, smap)
            D = get_data_matrix(dd, K, group_by_subcarrier=False)
            b = gfdm_modulate_block(D, taps, M, K, L, False)
            b = add_cyclic_starfix(b, cp_len, cs_len)
            b = pinch_block(b, window_taps)
            ref = np.concatenate((ref, frame_gap, preamble, b))
            data = np.concatenate((data, d))

        src = blocks.vector_source_c(data)
        mapper = gfdm.resource_mapper_cc(active, K, M, smap, True)
        mod = gfdm.simple_modulator_cc(M, K, L, taps)
        prefixer = gfdm.cyclic_prefixer_cc(block_len, cp_len, cs_len, ramp_len, window_taps)
        preambler = blocks.vector_insert_c(preamble, window_len + len(preamble), 0)
        gapper = blocks.vector_insert_c(frame_gap, frame_len + len(frame_gap), 0)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mapper, mod, prefixer, preambler, gapper, dst)
        # self.tb.connect(src, mapper, dst)
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
