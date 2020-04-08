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
import pmt
import gfdm_swig as gfdm
import numpy as np
from pygfdm.preamble import mapped_preamble
from pygfdm.mapping import get_subcarrier_map
from pygfdm.gfdm_modulation import modulate_mapped_gfdm_block
from pygfdm.utils import get_random_qpsk
from pygfdm.cyclic_prefix import pinch_cp_add_block


class qa_simple_preamble_sync_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        n_frames = 20
        timeslots = 9
        subcarriers = 128
        active_subcarriers = 110
        cp_len = subcarriers // 2
        smap = get_subcarrier_map(subcarriers, active_subcarriers)
        seed = 4711
        ftype = 'rrc'
        falpha = .5
        tag_key = 'frame_start'

        preamble, x_preamble = mapped_preamble(seed, ftype, falpha, active_subcarriers, subcarriers, smap, 2, cp_len, cp_len // 2)
        block_len = timeslots * subcarriers
        offset = len(preamble) + cp_len
        frame_len = len(preamble) + timeslots * subcarriers + cp_len

        data = np.array([], dtype=np.complex)
        ref = np.array([], dtype=np.complex)
        tags = []
        print('frame_len: ', frame_len)
        for i in range(n_frames):
            d_block = modulate_mapped_gfdm_block(get_random_qpsk(timeslots * active_subcarriers), timeslots, subcarriers, active_subcarriers, 2, falpha)
            frame = pinch_cp_add_block(d_block, timeslots, subcarriers, cp_len, cp_len // 2)
            frame = np.concatenate((preamble, frame))
            r = frame[offset:offset + block_len]
            ref = np.concatenate((ref, r))
            tag = gr.tag_t()
            tag.key = pmt.string_to_symbol(tag_key)
            tag.offset = len(data)
            tag.srcid = pmt.string_to_symbol('qa')
            tag.value = pmt.from_long(block_len)
            tags.append(tag)
            data = np.concatenate((data, frame))

        src = blocks.vector_source_c(data, False, 1, tags)
        cp_rm = gfdm.remove_prefix_cc(frame_len, block_len, offset, tag_key)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, cp_rm, snk)
        self.tb.run()

        # # check data
        res = np.array(snk.data())
        tags = snk.tags()
        self.assertComplexTuplesAlmostEqual(res, ref, 5)


if __name__ == '__main__':
    gr_unittest.run(qa_simple_preamble_sync_cc)
