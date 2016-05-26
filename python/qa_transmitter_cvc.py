#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2015 Andrej Rode.
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
import numpy as np
import gfdm_swig as gfdm
from pygfdm.modulation import gfdm_tx_fft2


class qa_transmitter_cvc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        # set up fg
        nsubcarrier = 8
        ntimeslots = 16
        filter_width = 2
        filter_alpha = 0.35
        src_data = np.array([np.complex(np.random.choice([-1, 1]), np.random.choice([-1, 1])) for i in
                             xrange(nsubcarrier * ntimeslots)])
        expected_result = gfdm_tx_fft2(src_data, 'rrc', filter_alpha, ntimeslots, nsubcarrier, filter_width, 1)
        src = blocks.vector_source_c(src_data)
        tm = gfdm.transmitter_cvc(nsubcarrier, ntimeslots, filter_width, filter_alpha)
        dst = blocks.vector_sink_c(vlen=nsubcarrier * ntimeslots)
        self.tb.connect(src, tm)
        self.tb.connect(tm, dst)
        self.tb.run()
        # check data
        result_data = dst.data()
        self.assertComplexTuplesAlmostEqual(expected_result, result_data, 2)

    def test_002_t(self):
        # set up fg
        nsubcarrier = 64
        ntimeslots = 32
        filter_width = 2
        filter_alpha = 0.35
        src_data = np.array([np.complex(np.random.choice([-1, 1]), np.random.choice([-1, 1])) for i in
                             xrange(nsubcarrier * ntimeslots)])
        expected_result = gfdm_tx_fft2(src_data, 'rrc', filter_alpha, ntimeslots, nsubcarrier, filter_width, 1)
        src = blocks.vector_source_c(src_data)
        tm = gfdm.transmitter_cvc(nsubcarrier, ntimeslots, filter_width, filter_alpha)
        dst = blocks.vector_sink_c(vlen=nsubcarrier * ntimeslots)
        self.tb.connect(src, tm)
        self.tb.connect(tm, dst)
        self.tb.run()
        # check data
        result_data = dst.data()
        self.assertComplexTuplesAlmostEqual(expected_result, result_data, 3)


if __name__ == '__main__':
    gr_unittest.run(qa_transmitter_cvc, "qa_transmitter_cvc.xml")
