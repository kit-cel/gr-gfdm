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
import gfdm_swig as gfdm

class qa_transmitter_cvc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        # set up fg
        nsubcarrier = 4
        ntimeslots = 4
        filter_width = 2
        filter_alpha = 0.35
        src_data = (1+1j,1-1j,-1-1j,-1+1j,1+1j,1-1j,1+1j,1+1j,-1-1j,1-1j,-1+1j,-1-1j,1+1j,1-1j,1+1j,1-1j)
        expected_result = (-1.92239433e-01-0.06435655j,  -7.79145448e-02-1.42576238j,
        -1.03862842e-01+0.06672548j,  -2.01631722e-01+0.47144911j,
        -1.18147020e+00-0.59032j   ,   3.42442003e-01+0.06587994j,
         4.74627623e-01+0.3040393j ,   7.49309212e-01+0.60528329j,
        -7.18202886e-01-0.59032j   ,   6.36604711e-02+0.08897928j,
        -8.45392403e-01-0.3040393j ,  -1.41490093e-01-0.15499258j,
        -1.18147020e+00-0.06435655j,  -8.49657987e-04+0.28888835j,
         4.74627623e-01-0.06672548j,  -7.88491259e-02+0.06027499j)
        src = blocks.vector_source_c(src_data)
        tm = gfdm.transmitter_cvc(nsubcarrier,ntimeslots,filter_width,filter_alpha)
        dst = blocks.vector_sink_f()
        self.tb.connect(src,tm)
        self.tb.connect(tm,dst)
        self.tb.run ()
        # check data
        result_data = dst.data()
        self.assertComplexTuplesAlmostEqual(expected_result, result_data, 6)


if __name__ == '__main__':
    gr_unittest.run(qa_transmitter_cvc, "qa_transmitter_cvc.xml")
