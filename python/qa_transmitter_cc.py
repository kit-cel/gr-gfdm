#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 Johannes Demel.
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
import gfdm_python as gfdm
from pygfdm.cyclic_prefix import get_window_len, get_raised_cosine_ramp
from pygfdm.cyclic_prefix import add_cyclic_starfix, pinch_block
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.utils import get_random_qpsk
from pygfdm.preamble import mapped_preamble
from pygfdm.mapping import get_subcarrier_map, map_to_waveform_resources
from pygfdm.mapping import get_data_matrix
from pygfdm.gfdm_modulation import gfdm_modulate_block
import numpy as np

# import os
# print('Blocked waiting for GDB attach (pid = %d)' % (os.getpid(),))
# input('Press Enter to continue: ')


def generate_reference_frame(symbols, timeslots, subcarriers, active_subcarriers, subcarrier_map,
                             taps, overlap, cp_len, cs_len, window_taps, cyclic_shifts, preambles):
    dd = map_to_waveform_resources(symbols, active_subcarriers,
                                   subcarriers, subcarrier_map)
    D = get_data_matrix(dd, subcarriers, group_by_subcarrier=False)
    b = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                            overlap, False)
    frame = []
    for cs, p in zip(cyclic_shifts, preambles):
        data = np.roll(b, cs)
        data = add_cyclic_starfix(data, cp_len, cs_len)
        data = pinch_block(data, window_taps)
        frame.append(np.concatenate((p, data)))
    return frame


def get_full_preambles(filtertype, alpha, active_subcarriers, subcarriers,
                       subcarrier_map, overlap, cp_len, ramp_len,
                       cyclic_shifts):
    seed = 4711
    preambles = []
    for cs in cyclic_shifts:
        p = mapped_preamble(seed, filtertype, alpha, active_subcarriers, subcarriers,
                            subcarrier_map, overlap, cp_len, ramp_len,
                            use_zadoff_chu=True, cyclic_shift=cs)
        preambles.append(p[0])
    return preambles


class qa_transmitter_cc (gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()
        self.filter_type = 'rrc'
        self.alpha = .5
        self.overlap = 2

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        np.set_printoptions(precision=2)
        n_frames = 7
        active_subcarriers = 52
        timeslots = 9
        subcarriers = 64
        cp_len = subcarriers // 4
        cs_len = cp_len // 2

        taps = get_frequency_domain_filter(self.filter_type, self.alpha, timeslots,
                                           subcarriers, self.overlap)

        window_taps = get_raised_cosine_ramp(
            cs_len, get_window_len(cp_len, timeslots, subcarriers, cs_len))

        smap = get_subcarrier_map(subcarriers, active_subcarriers, True)
        preambles = get_full_preambles(self.filter_type, self.alpha, active_subcarriers,
                                       subcarriers, smap, self.overlap, cp_len, cs_len, [0, ])

        frame_size = preambles[0].size + cp_len + timeslots * subcarriers + cs_len

        ref = np.array([], dtype=np.complex)
        data = np.array([], dtype=np.complex)

        for i in range(n_frames):
            d = get_random_qpsk(active_subcarriers * timeslots)
            frame = generate_reference_frame(d, timeslots, subcarriers, active_subcarriers, smap,
                                             taps, self.overlap, cp_len, cs_len, window_taps, [0, ], preambles)

            ref = np.concatenate((ref, frame[0]))
            data = np.concatenate((data, d))

        src = blocks.vector_source_c(data)
        dut = gfdm.transmitter_cc(timeslots, subcarriers, active_subcarriers,
                                  cp_len, cs_len, cs_len, smap, True,
                                  self.overlap, taps, window_taps, [0, ], preambles, "packet_len")
        dst = blocks.vector_sink_c()

        self.tb.connect(src, dut, dst)
        self.tb.run()
        res = np.array(dst.data())[0:len(ref)]
        self.assertComplexTuplesAlmostEqual(ref, res, 5)

        tags = dst.tags()
        for i, t in enumerate(tags):
            print(f't={i}, offset={t.offset}, value={pmt.to_python(t.value)}')
            self.assertEqual(t.offset, i * frame_size)
            self.assertEqual(pmt.to_python(t.value), frame_size)

    def test_002_cyclic_delay_diversity(self):
        np.set_printoptions(precision=2)
        n_frames = 7
        active_subcarriers = 52
        timeslots = 9
        subcarriers = 64
        cp_len = subcarriers // 4
        cs_len = cp_len // 2
        cyclic_shifts = [0, 3, 7, 8]

        taps = get_frequency_domain_filter(self.filter_type, self.alpha, timeslots,
                                           subcarriers, self.overlap)

        window_taps = get_raised_cosine_ramp(
            cs_len, get_window_len(cp_len, timeslots, subcarriers, cs_len))

        smap = get_subcarrier_map(subcarriers, active_subcarriers, True)
        preambles = get_full_preambles(self.filter_type, self.alpha, active_subcarriers,
                                       subcarriers, smap, self.overlap, cp_len, cs_len, cyclic_shifts)

        frame_size = preambles[0].size + cp_len + timeslots * subcarriers + cs_len

        ref = [np.array([], dtype=np.complex) for _ in cyclic_shifts]
        data = np.array([], dtype=np.complex)

        for i in range(n_frames):
            d = get_random_qpsk(active_subcarriers * timeslots)
            frame = generate_reference_frame(d, timeslots, subcarriers, active_subcarriers, smap,
                                             taps, self.overlap, cp_len, cs_len, window_taps,
                                             cyclic_shifts, preambles)

            ref = [np.concatenate((r, f)) for r, f in zip(ref, frame)]
            data = np.concatenate((data, d))

        src = blocks.vector_source_c(data)
        dut = gfdm.transmitter_cc(timeslots, subcarriers, active_subcarriers,
                                  cp_len, cs_len, cs_len, smap, True,
                                  self.overlap, taps, window_taps, cyclic_shifts, preambles, "packet_len")
        snks = [blocks.vector_sink_c() for _ in cyclic_shifts]

        self.tb.connect(src, dut)
        for i, s in enumerate(snks):
            self.tb.connect((dut, i), s)
        self.tb.run()

        for snk, refport in zip(snks, ref):
            res = np.array(snk.data())[0:refport.size]
            self.assertComplexTuplesAlmostEqual(refport, res, 5)

        for j, snk in enumerate(snks):
            tags = snk.tags()
            for i, t in enumerate(tags):
                print(f'p={j}, t={i}, offset={t.offset}, value={pmt.to_python(t.value)}')
                self.assertEqual(t.offset, i * frame_size)
                self.assertEqual(pmt.to_python(t.value), frame_size)


if __name__ == '__main__':
    gr_unittest.run(qa_transmitter_cc)
