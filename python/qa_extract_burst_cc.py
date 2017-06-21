#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2017 Johannes Demel.
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
import pmt


class qa_extract_burst_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        n_frames = 500
        burst_len = 383
        gap_len = 53
        tag_key = 'energy_start'

        data = np.arange(burst_len)
        ref = np.array([], dtype=np.complex)
        tags = []
        for i in range(n_frames):
            frame = np.ones(burst_len) * (i + 1)
            ref = np.concatenate((ref, frame))
            tag = gr.tag_t()
            tag.key = pmt.string_to_symbol(tag_key)
            tag.offset = burst_len + i * (burst_len + gap_len)
            tag.srcid = pmt.string_to_symbol('qa')
            tag.value = pmt.PMT_T
            tags.append(tag)
            data = np.concatenate((data, frame, np.zeros(gap_len)))
        # print(np.reshape(data, (-1, burst_len)))
        # print('data len', len(data), 'ref len', len(ref))

        src = blocks.vector_source_c(data, False, 1, tags)
        burster = gfdm.extract_burst_cc(burst_len, tag_key)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, burster, snk)
        self.tb.run()

        res = np.array(snk.data())
        rx_tags = snk.tags()
        for i, t in enumerate(rx_tags):
            assert pmt.symbol_to_string(t.key) == tag_key
            assert pmt.to_long(t.value) == burst_len
            assert pmt.symbol_to_string(t.srcid) == burster.name()
            assert t.offset == i * burst_len
            # print t.offset, t.value

        # check data
        self.assertComplexTuplesAlmostEqual(ref, res)


if __name__ == '__main__':
    gr_unittest.run(qa_extract_burst_cc)
