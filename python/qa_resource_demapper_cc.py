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
import gfdm_python as gfdm
from pygfdm.mapping import get_subcarrier_map
from pygfdm.utils import get_random_qpsk
import numpy as np


class qa_resource_demapper_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_timeslot_first(self):
        timeslots = 9
        subcarriers = 32
        active_subcarriers = 20
        subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers)

        data = get_random_qpsk(10 * timeslots * active_subcarriers)
        src = blocks.vector_source_c(data)
        mapper = gfdm.resource_mapper_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, True)
        demapper = gfdm.resource_demapper_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, True)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, mapper, demapper, snk)
        self.tb.run()
        # check data
        res = np.array(snk.data())

        self.assertComplexTuplesAlmostEqual(data, res)

    def test_002_subcarrier_first(self):
        timeslots = 9
        subcarriers = 32
        active_subcarriers = 20
        subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers)

        data = get_random_qpsk(10 * timeslots * active_subcarriers)
        src = blocks.vector_source_c(data)
        mapper = gfdm.resource_mapper_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, False)
        demapper = gfdm.resource_demapper_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, False)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, mapper, demapper, snk)
        self.tb.run()
        # check data
        res = np.array(snk.data())

        self.assertComplexTuplesAlmostEqual(data, res)


if __name__ == '__main__':
    gr_unittest.run(qa_resource_demapper_cc, "qa_resource_demapper_cc.xml")
