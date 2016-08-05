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

from gfdm.pygfdm import filters
from gnuradio import gr, gr_unittest
from gnuradio import blocks,digital
import gfdm_swig as gfdm
import numpy
def struct(data): return type('Struct', (object,), data)()

class qa_advanced_receiver_sb_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        self.gfdm_var = gfdm_var = struct({'subcarriers': 64, 'timeslots': 8, 'alpha': 0.5, 'overlap': 2, })
        self.gfdm_constellation = gfdm_constellation = digital.constellation_qpsk().base()
        self.f_taps = f_taps = filters.get_frequency_domain_filter('rrc', 1.0, gfdm_var.timeslots, gfdm_var.subcarriers, gfdm_var.overlap)
        self.random_bits = blocks.vector_source_b(map(int, numpy.random.randint(0, len(gfdm_constellation.points()), 100*gfdm_var.timeslots*gfdm_var.subcarriers)), False)
        self.bits_to_symbols = digital.chunks_to_symbols_bc((gfdm_constellation.points()), 1)
        self.mod = gfdm.simple_modulator_cc(gfdm_var.timeslots, gfdm_var.subcarriers, gfdm_var.overlap, f_taps)
        self.scale = blocks.multiply_const_vcc((1./gfdm_var.subcarriers, ))
        self.demod = gfdm.advanced_receiver_sb_cc(gfdm_var.timeslots, gfdm_var.subcarriers, gfdm_var.overlap, 64, f_taps, gfdm_constellation)
        self.tx_symbols = blocks.vector_sink_c()
        self.rx_symbols = blocks.vector_sink_c()
        self.tb.connect((self.random_bits,0),(self.bits_to_symbols,0))
        self.tb.connect((self.bits_to_symbols,0),(self.tx_symbols,0))
        self.tb.connect((self.bits_to_symbols,0),(self.mod,0))
        self.tb.connect((self.mod,0),(self.scale,0))
        self.tb.connect((self.scale,0),(self.demod,0))
        self.tb.connect((self.demod,0),(self.rx_symbols,0))
        self.tb.run ()
        ref = numpy.array(self.tx_symbols.data())
        res = numpy.array(self.rx_symbols.data())
        self.assertComplexTuplesAlmostEqual(ref, res, 2)


if __name__ == '__main__':
    gr_unittest.run(qa_advanced_receiver_sb_cc, "qa_advanced_receiver_sb_cc.xml")
