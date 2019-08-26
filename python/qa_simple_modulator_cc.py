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
from pygfdm.utils import get_random_qpsk, calculate_signal_energy
import numpy as np


class qa_simple_modulator_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        alpha = .5
        M = 8
        K = 4
        L = 2
        taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
        taps /= np.sqrt(calculate_signal_energy(taps) / M)
        # data = np.repeat(np.arange(1, K + 1), M)
        data = get_random_qpsk(M * K)
        D = get_data_matrix(data, K, group_by_subcarrier=False)

        src = blocks.vector_source_c(data)
        mod = gfdm.simple_modulator_cc(M, K, L, taps)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, dst)
        # set up fg
        self.tb.run()
        # check data
        res = np.array(dst.data())
        # res /= M * K

        # F = gfdm_transform_subcarriers_to_fd(D, M)
        # F = gfdm_upsample_subcarriers_in_fd(F, K, L)
        # F = gfdm_filter_subcarriers_in_fd(F, taps, M, K, L)

        ref = gfdm_modulate_block(D, taps, M, K, L, False)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)

    def test_002_big_data(self):
        print("big data test")
        reps = 5
        alpha = .5
        M = 127
        K = 16
        L = 4
        taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
        taps /= np.sqrt(calculate_signal_energy(taps) / M)
        data = np.array([], dtype=np.complex)
        ref = np.array([], dtype=np.complex)
        for i in range(reps):
            d = get_random_qpsk(M * K)
            D = get_data_matrix(d, K, group_by_subcarrier=False)
            ref = np.append(ref, gfdm_modulate_block(D, taps, M, K, L, False))
            data = np.append(data, d)

        src = blocks.vector_source_c(data)
        mod = gfdm.simple_modulator_cc(M, K, L, taps)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, dst)
        # set up fg
        self.tb.run()
        # check data
        res = np.array(dst.data())
        # res /= M * K

        self.assertComplexTuplesAlmostEqual(ref, res, 2)


if __name__ == '__main__':
    gr_unittest.run(qa_simple_modulator_cc)
