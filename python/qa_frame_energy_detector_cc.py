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
import pmt

class qa_frame_energy_detector_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_setup(self):
        alpha = 3.
        average_len = 8
        frame_len = 25
        backoff_len = average_len
        tag_key = 'enertest'
        detector = gfdm.frame_energy_detector_cc(alpha, average_len, frame_len, backoff_len, tag_key)

        self.assertEqual(backoff_len, detector.backoff_len())
        self.assertAlmostEqual(alpha, detector.alpha())

        alpha = 5.0
        detector.set_alpha(alpha)
        backoff_len = 16
        detector.set_backoff_len(backoff_len)
        self.assertEqual(backoff_len, detector.backoff_len())
        self.assertAlmostEqual(alpha, detector.alpha())

    def test_002_simple_passthru(self):
        alpha = 3.
        average_len = 8
        frame_len = 25
        backoff_len = 0
        tag_key = 'enertest'
        # set up fg
        detector = gfdm.frame_energy_detector_cc(alpha, average_len, frame_len, backoff_len, tag_key)

        data = np.zeros(100)
        data[50:] = 100
        src = blocks.vector_source_c(data)
        snk = blocks.vector_sink_c()

        self.tb.connect(src, detector, snk)
        self.tb.run()

        # check data
        res = snk.data()
        # make sure output is a multiple of average_len!
        self.assertTrue(len(res) == frame_len)
        tags = snk.tags()

        for t in tags:
            print 'srcid {}, key {}, offset {}, value {}'.format(t.srcid, t.key, t.offset, t.value)

        p = tags[0].offset
        self.assertEqual(p, 0)

        self.assertComplexTuplesAlmostEqual(res, data[48:48 + len(res)])

    def test_003_frame_detect(self):
        alpha = 3.
        average_len = 8
        frame_len = 25
        backoff_len = 4
        tag_key = 'enertest'
        # set up fg
        detector = gfdm.frame_energy_detector_cc(alpha, average_len, frame_len, backoff_len, tag_key)

        data = np.zeros(100)
        data[50:] = 100
        src = blocks.vector_source_c(data)
        snk = blocks.vector_sink_c()

        self.tb.connect(src, detector, snk)
        self.tb.run()

        # check data
        res = np.array(snk.data())
        tags = snk.tags()

        for t in tags:
            print 'srcid {}, key {}, offset {}, value {}'.format(t.srcid, t.key, t.offset, t.value)

        p = tags[0].offset
        self.assertEqual(p, 0)
        v = pmt.to_long(tags[0].value)
        fl = frame_len + 2 * backoff_len
        self.assertEqual(int(v), int(fl))

        self.assertComplexTuplesAlmostEqual(res, data[44:44 + len(res)])

    def test_004_frames(self):
        alpha = 3.
        average_len = 8
        frame_len = 25
        backoff_len = 2 * average_len
        tag_key = 'enertest'
        # set up fg
        detector = gfdm.frame_energy_detector_cc(alpha, average_len, frame_len, backoff_len, tag_key)

        test_frame = np.zeros(frame_len + 2 * backoff_len)
        test_frame[backoff_len:-backoff_len] = 100.
        data = np.zeros(1000)
        p = 6 * average_len
        data[p:p + len(test_frame)] = test_frame
        p = 60 * average_len
        data[p:p + len(test_frame)] = test_frame
        ref = np.tile(test_frame, 2)

        src = blocks.vector_source_c(data)
        snk = blocks.vector_sink_c()

        self.tb.connect(src, detector, snk)
        self.tb.run()

        # check data
        res = np.array(snk.data())
        # make sure output is a multiple of average_len!
        tags = snk.tags()

        p = tags[0].offset
        self.assertEqual(p, 0)
        v = pmt.to_long(tags[0].value)
        fl = frame_len + 2 * backoff_len
        self.assertEqual(int(v), int(fl))
        self.assertComplexTuplesAlmostEqual(res[0:len(test_frame)], test_frame)

    def test_005_long_frame(self):
        alpha = 3.
        average_len = 32
        frame_len = 5000
        backoff_len = 2 * average_len
        tag_key = 'enertest'
        # set up fg
        detector = gfdm.frame_energy_detector_cc(alpha, average_len, frame_len, backoff_len, tag_key)

        test_frame = np.zeros(frame_len + 2 * backoff_len)
        test_frame[backoff_len:-backoff_len] = 100.
        data = np.zeros(10 * frame_len)
        p = 6 * average_len
        data[p:p + len(test_frame)] = test_frame

        src = blocks.vector_source_c(data)
        snk = blocks.vector_sink_c()

        self.tb.connect(src, detector, snk)
        self.tb.run()

        # check data
        res = np.array(snk.data())
        # make sure output is a multiple of average_len!
        tags = snk.tags()
        for t in tags:
            print 'srcid {}, key {}, offset {}, value {}'.format(t.srcid, t.key, t.offset, t.value)

        p = tags[0].offset
        self.assertEqual(p, 0)
        v = pmt.to_long(tags[0].value)
        fl = frame_len + 2 * backoff_len
        self.assertEqual(int(v), int(fl))
        self.assertComplexTuplesAlmostEqual(res[0:len(test_frame)], test_frame)


if __name__ == '__main__':
    # gr_unittest.run(qa_frame_energy_detector_cc, "qa_frame_energy_detector_cc.xml")
    gr_unittest.run(qa_frame_energy_detector_cc)
