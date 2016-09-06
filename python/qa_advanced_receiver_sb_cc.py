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
from gnuradio import blocks, digital
import gfdm_swig as gfdm
from pygfdm import filters, utils
from pygfdm.gfdm_receiver import gfdm_demodulate_block
import numpy as np


def struct(data): return type('Struct', (object,), data)()


class qa_advanced_receiver_sb_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_simple_receiver(self):
        # make sure advanced receiver works like simple receiver in case no IC iterations are applied!
        reps = 5
        alpha = .5
        M = 127
        K = 16
        L = 2
        taps = filters.get_frequency_domain_filter('rrc', alpha, M, K, L)
        data = np.array([], dtype=np.complex)
        ref = np.array([], dtype=np.complex)
        for i in xrange(reps):
            d = utils.get_random_qpsk(M * K)
            ref = np.append(ref, gfdm_demodulate_block(d, taps, K, M, L))
            data = np.append(data, d)
        # print data
        # print ref
        # print "MAXIMUM ref value: ", np.max(abs(ref))

        src = blocks.vector_source_c(data)
        gfdm_constellation = digital.constellation_qpsk().base()
        mod = gfdm.advanced_receiver_sb_cc(M, K, L, 0,
                                           taps, gfdm_constellation)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, dst)
        # set up fg
        self.tb.run()
        # check data
        res = np.array(dst.data())

        self.assertComplexTuplesAlmostEqual(ref, res, 4)

    def test_002_t(self):
        n_frames = 1
        self.gfdm_var = gfdm_var = struct({'subcarriers': 64, 'timeslots': 9, 'alpha': 0.5, 'overlap': 2,})
        self.gfdm_constellation = gfdm_constellation = digital.constellation_qpsk().base()
        self.f_taps = f_taps = filters.get_frequency_domain_filter('rrc', 1.0, gfdm_var.timeslots, gfdm_var.subcarriers,
                                                                   gfdm_var.overlap)
        self.random_bits = blocks.vector_source_b(map(int, np.random.randint(0, len(gfdm_constellation.points()),
                                                                                n_frames * gfdm_var.timeslots * gfdm_var.subcarriers)),
                                                  False)
        self.bits_to_symbols = digital.chunks_to_symbols_bc((gfdm_constellation.points()), 1)
        self.mod = gfdm.simple_modulator_cc(gfdm_var.timeslots, gfdm_var.subcarriers, gfdm_var.overlap, f_taps)
        self.demod = gfdm.advanced_receiver_sb_cc(gfdm_var.timeslots, gfdm_var.subcarriers, gfdm_var.overlap, 64,
                                                  f_taps, gfdm_constellation)
        self.tx_symbols = blocks.vector_sink_c()
        self.rx_symbols = blocks.vector_sink_c()
        self.tb.connect((self.random_bits, 0), (self.bits_to_symbols, 0))
        self.tb.connect((self.bits_to_symbols, 0), (self.tx_symbols, 0))
        self.tb.connect((self.bits_to_symbols, 0), (self.mod, 0))
        self.tb.connect((self.mod, 0), (self.demod, 0))
        self.tb.connect((self.demod, 0), (self.rx_symbols, 0))
        self.tb.run()
        ref = np.array(self.tx_symbols.data())
        res = np.array(self.rx_symbols.data())
        # more or less make sure all symbols have their correct sign.
        self.assertComplexTuplesAlmostEqual(ref, res, 2)


if __name__ == '__main__':
    # gr_unittest.run(qa_advanced_receiver_sb_cc, "qa_advanced_receiver_sb_cc.xml")
    gr_unittest.run(qa_advanced_receiver_sb_cc)
