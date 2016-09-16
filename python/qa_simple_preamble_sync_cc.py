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
from pygfdm.gfdm_modulation import modulate_mapped_gfdm_block
from pygfdm.utils import get_random_qpsk
from pygfdm.cyclic_prefix import pinch_cp_add_block
import pmt


class qa_simple_preamble_sync_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        n_frames = 100
        timeslots = 9
        subcarriers = 128
        active_subcarriers = 110
        cp_len = subcarriers // 2
        smap = get_subcarrier_map(subcarriers, active_subcarriers)
        seed = 4711
        ftype = 'rrc'
        falpha = .5

        preamble, x_preamble = mapped_preamble(seed, ftype, falpha, active_subcarriers, subcarriers, smap, 2, cp_len, cp_len // 2)
        frame_len = len(preamble) + timeslots * subcarriers + cp_len
        frame_gap = np.zeros(frame_len, dtype=np.complex)
        data = frame_gap
        ref = np.array([], dtype=np.complex)
        for i in range(n_frames):
            d_block = modulate_mapped_gfdm_block(get_random_qpsk(timeslots * active_subcarriers), timeslots, subcarriers, active_subcarriers, 2, falpha)
            frame = pinch_cp_add_block(d_block, timeslots, subcarriers, cp_len, cp_len // 2)
            frame = np.concatenate((preamble, frame))
            ref = np.concatenate((ref, frame))
            data = np.concatenate((data, frame, frame_gap))
        # data = np.zeros(4000, dtype=np.complex)
        # frame_start = 800
        backoff = 80
        # data[frame_start - cp_len:frame_start - cp_len + len(preamble)] = preamble
        print 'qa', len(data)

        # tag = gr.tag_t()
        # tag.key = pmt.string_to_symbol('energy')
        # tag.offset = frame_start - cp_len - backoff
        # tag.srcid = pmt.string_to_symbol('qa')
        # tag.value = pmt.from_long(len(preamble) + 2 * backoff)

        # FIXME: the following 2 lines are subject to a bugreport at the moment. Will see how to fix the QA test later.
        # src = blocks.vector_source_c(data, (tag, ))
        src = blocks.vector_source_c(data)
        e_detector = gfdm.frame_energy_detector_cc(20., 32, frame_len, backoff, 'energy')
        detector = gfdm.simple_preamble_sync_cc(frame_len, subcarriers, cp_len, x_preamble, 'energy', 'frame')
        snk = blocks.vector_sink_c()

        self.tb.connect(src, e_detector, detector, snk)
        self.tb.run()
        # check data
        res = np.array(snk.data())
        tags = snk.tags()
        for t in tags:
            print 'srcid {}, key {}, offset {}, value {}'.format(t.srcid, t.key, t.offset, t.value)

        self.assertComplexTuplesAlmostEqual(res, ref, 5)
        # rtags = snk.tags()
        # for t in rtags:
        #     print gr.tag_to_python(t)
        # import matplotlib.pyplot as plt
        # plt.plot(np.abs(data))
        # plt.plot(np.abs(res))
        # plt.show()


if __name__ == '__main__':
    # gr_unittest.run(qa_simple_preamble_sync_cc, "qa_simple_preamble_sync_cc.xml")
    gr_unittest.run(qa_simple_preamble_sync_cc)
