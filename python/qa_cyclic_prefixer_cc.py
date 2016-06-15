#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2016 <+YOU OR YOUR COMPANY+>.
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
from pygfdm.cyclic_prefix import add_cyclic_prefix, pinch_block, get_raised_cosine_ramp, get_window_len


class qa_cyclic_prefixer_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_init(self):
        # check if prefixer is properly ctor'ed / dtor'ed
        prefixer = gfdm.cyclic_prefixer_cc(4, 'gfdm_frame')
        prefixer = gfdm.cyclic_prefixer_cc(4, 4, 16 * 8, np.arange(4 * 2), 'gfdm_frame')
        prefixer = gfdm.cyclic_prefixer_cc(4, 4, 16 * 8, np.arange(16 * 8 + 4), 'gfdm_frame')
        try:
            prefixer = gfdm.cyclic_prefixer_cc(4, 4, 16 * 8, np.arange(16 * 8), 'gfdm_frame')
            raise ValueError('invalid parameter set, but passed anyway!')
        except:
            # expected behavior!
            pass

    def test_002_simple_cp(self):
        tag_key = 'gfdm_block'
        block_len = 48
        cp_len = 8
        data = np.arange(block_len, dtype=np.complex) + 1
        ref = add_cyclic_prefix(data, cp_len)

        prefixer = gfdm.cyclic_prefixer_cc(cp_len, tag_key)
        tagger = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 1, block_len, tag_key)
        src = blocks.vector_source_c(data)
        dst = blocks.vector_sink_c()
        self.tb.connect(src, tagger, prefixer, dst)
        self.tb.run()
        res = np.array(dst.data())
        self.assertComplexTuplesAlmostEqual(res, ref)

    def test_003_block_pinching(self):
        tag_key = 'gfdm_block'
        n_subcarriers = 8
        n_timeslots = 8
        block_len = n_subcarriers * n_timeslots
        cp_len = 8
        ramp_len = 4
        window_len = get_window_len(cp_len, n_timeslots, n_subcarriers)
        window_taps = get_raised_cosine_ramp(ramp_len, window_len)
        data = np.arange(block_len, dtype=np.complex) + 1
        ref = add_cyclic_prefix(data, cp_len)
        ref = pinch_block(ref, window_taps)

        prefixer = gfdm.cyclic_prefixer_cc(cp_len, ramp_len, block_len, window_taps, tag_key)
        tagger = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 1, block_len, tag_key)
        src = blocks.vector_source_c(data)
        dst = blocks.vector_sink_c()
        self.tb.connect(src, tagger, prefixer, dst)
        self.tb.run()

        res = np.array(dst.data())

        self.assertComplexTuplesAlmostEqual(res, ref, 4)


if __name__ == '__main__':
    # gr_unittest.run(qa_cyclic_prefixer_cc, "qa_cyclic_prefixer_cc.xml")
    gr_unittest.run(qa_cyclic_prefixer_cc)
