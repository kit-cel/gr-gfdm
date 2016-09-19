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
import numpy as np
from pygfdm.mapping import map_to_waveform_resources
from pygfdm.utils import get_random_qpsk

class qa_resource_mapper_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        # set up fg
        n_frames = 3
        active_subcarriers = 110
        subcarriers = 128
        timeslots = 205
        smap = np.arange(active_subcarriers) + (subcarriers - active_subcarriers) // 2

        data = get_random_qpsk(active_subcarriers * timeslots)
        ref = map_to_waveform_resources(data, active_subcarriers, subcarriers, smap, True)
        for i in range(n_frames - 1):
            d = get_random_qpsk(active_subcarriers * timeslots)
            data = np.concatenate((data, d))
            ref = np.concatenate((ref, map_to_waveform_resources(d, active_subcarriers, subcarriers, smap, True)))

        src = blocks.vector_source_c(data)
        mapper = gfdm.resource_mapper_cc(timeslots, subcarriers, active_subcarriers, smap, True)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, mapper, snk)
        self.tb.run()
        # check data
        res = snk.data()
        self.assertComplexTuplesAlmostEqual(ref, res)


if __name__ == '__main__':
    # gr_unittest.run(qa_resource_mapper_cc, "qa_resource_mapper_cc.xml")
    gr_unittest.run(qa_resource_mapper_cc)
