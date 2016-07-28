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
from gnuradio import gr, gr_unittest
import gfdm_swig as gfdm
from gfdm.pygfdm import gfdm_modulation, gfdm_receiver, utils, receiver

class qa_pygfdm_module (gr_unittest.TestCase):

    def setUp (self):
        return 0

    def tearDown (self):
        return 0

    def test_001_transceiver_legacy_00(self):
        return 0
    def test_transceiver_legacy_00():
        K = 32
        M = 8
        overlap = 2
        alpha = 0.5
        oversampling_factor = 1

        tests = 100
        for t in xrange(tests):
            d = utils.get_random_qpsk(M*K)
            tx = gfdm_modulation.gfdm_gr_modulator(d, 'rrc', alpha, M, K, overlap)
            rx = receiver.gfdm_rx_fft2(
                tx,
                'rrc',
                alpha,
                M,
                K,
                overlap,
                oversampling_factor,
                4,
                16)
            if not (np.max(np.abs(d-rx)) < 1e-2):
                raise RuntimeError('Input and Output deviate')



if __name__ == '__mail__':
    gr_unittest.run(qa_pygfdm_module, "qa_pygfdm_module.xml")
