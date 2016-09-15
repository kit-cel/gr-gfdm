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
from pygfdm.preamble import mapped_preamble
from pygfdm.mapping import get_subcarrier_map
import pmt


class qa_simple_preamble_sync_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        subcarriers = 128
        active_subcarriers = 110
        cp_len = subcarriers // 2
        smap = get_subcarrier_map(subcarriers, active_subcarriers)
        preamble, x_preamble = mapped_preamble(4711, 'rrc', .5, active_subcarriers, subcarriers, smap, 2, cp_len, cp_len // 2)
        data = np.zeros(4000, dtype=np.complex)
        frame_start = 800
        backoff = 80
        data[frame_start - cp_len:frame_start - cp_len + len(preamble)] = preamble

        tag = gr.tag_t()
        tag.key = pmt.string_to_symbol('energy')
        tag.offset = frame_start - cp_len - backoff
        tag.srcid = pmt.string_to_symbol('qa')
        tag.value = pmt.from_long(len(preamble) + 2 * backoff)

        # FIXME: the following 2 lines are subject to a bugreport at the moment. Will see how to fix the QA test later.
        # src = blocks.vector_source_c(data, (tag, ))
        src = blocks.vector_source_c(data)

        detector = gfdm.simple_preamble_sync_cc(len(preamble), subcarriers, cp_len, x_preamble, 'energy', 'frame')
        snk = blocks.vector_sink_c()

        self.tb.connect(src, detector, snk)
        self.tb.run()
        # check data
        # res = snk.data()
        # rtags = snk.tags()
        # for t in rtags:
        #     print gr.tag_to_python(t)


if __name__ == '__main__':
    # gr_unittest.run(qa_simple_preamble_sync_cc, "qa_simple_preamble_sync_cc.xml")
    gr_unittest.run(qa_simple_preamble_sync_cc)
