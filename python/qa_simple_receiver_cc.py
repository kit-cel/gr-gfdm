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

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import gfdm_swig as gfdm
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.gfdm_receiver import gfdm_demodulate_block
from pygfdm.utils import get_random_qpsk, calculate_signal_energy
import numpy as np


class qa_simple_receiver_cc(gr_unittest.TestCase):
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
        data = get_random_qpsk(M * K)
        src = blocks.vector_source_c(data)
        mod = gfdm.simple_receiver_cc(M, K, L, taps)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, dst)
        self.tb.run()
        res = np.array(dst.data())

        ref = gfdm_demodulate_block(data, taps, K, M, L)
        # print calculate_signal_energy(ref), calculate_signal_energy(res)
        res *= np.sqrt(calculate_signal_energy(ref) / calculate_signal_energy(res))
        self.assertComplexTuplesAlmostEqual(ref, res, 5)

    def test_002_big_data(self):
        print("big data test")
        reps = 5
        alpha = .5
        M = 127
        K = 16
        L = 2
        taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
        data = np.array([], dtype=np.complex)
        ref = np.array([], dtype=np.complex)
        for i in xrange(reps):
            d = get_random_qpsk(M * K)
            ref = np.append(ref, gfdm_demodulate_block(d, taps, K, M, L))
            data = np.append(data, d)
        # print data
        # print ref
        # print "MAXIMUM ref value: ", np.max(abs(ref))

        src = blocks.vector_source_c(data)
        mod = gfdm.simple_receiver_cc(M, K, L, taps)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, dst)
        # set up fg
        self.tb.run()
        # check data
        res = np.array(dst.data())

        self.assertComplexTuplesAlmostEqual(ref, res, 4)


if __name__ == '__main__':
    # gr_unittest.run(qa_simple_receiver_cc, "qa_simple_receiver_cc.xml")
    gr_unittest.run(qa_simple_receiver_cc)
