#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Andrej Rode.
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
import numpy as np
import gfdm.modulation as mod
import gfdm_swig as gfdms
from pprint import pprint

class qa_framer_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        nsubcarrier = 64
        ntimeslots = 16
        src_data = np.array([np.complex(np.random.choice([-1,1]),np.random.choice([-1,1])) for i in xrange(nsubcarrier*ntimeslots)])
        src = blocks.vector_source_c(src_data,vlen=1)
        stts = blocks.stream_to_tagged_stream(
            gr.sizeof_gr_complex,vlen=1,
            packet_len=nsubcarrier*ntimeslots,len_tag_key="frame_len")
        fr = gfdms.framer_cc(nsubcarrier,ntimeslots,False,[],"frame_len")
        dst = blocks.vector_sink_c(vlen=1)
        expected_result = mod.reshape_input(src_data, ntimeslots, nsubcarrier)
        self.tb.connect(src,stts)
        self.tb.connect(stts,fr)
        self.tb.connect(fr,dst)
        self.tb.run()
        result_data = dst.data()
        tags = dst.tags()
        for tag in tags:
            ptag = gr.tag_to_python(tag)
            pprint(vars(ptag))
        self.assertComplexTuplesAlmostEqual(expected_result,result_data, 6)

    def test_002_t (self):
        nsubcarrier = 16
        ntimeslots = 64
        sync_data = np.array([np.complex(np.random.choice([-1,1]),np.random.choice([-1,1])) for i in xrange(nsubcarrier)])
        src_data = np.array([np.complex(np.random.choice([-1,1]),np.random.choice([-1,1])) for i in xrange(nsubcarrier*ntimeslots)])
        src = blocks.vector_source_c(src_data,vlen=1)
        stts = blocks.stream_to_tagged_stream(
            gr.sizeof_gr_complex,vlen=1,
            packet_len=nsubcarrier*ntimeslots,len_tag_key="frame_len")
        fr = gfdms.framer_cc(nsubcarrier,ntimeslots,True,sync_data,"frame_len")
        dst = blocks.vector_sink_c(vlen=1)
        expected_result = np.concatenate(
            (mod.reshape_input(np.tile(sync_data,2),2,nsubcarrier),
             mod.reshape_input(src_data, ntimeslots, nsubcarrier)))
        self.tb.connect(src,stts)
        self.tb.connect(stts,fr)
        self.tb.connect(fr,dst)
        self.tb.run()
        result_data = dst.data()
        tags = dst.tags()
        for tag in tags:
            ptag = gr.tag_to_python(tag)
            pprint(vars(ptag))
        self.assertComplexTuplesAlmostEqual(expected_result,result_data[0:nsubcarrier*ntimeslots+2*nsubcarrier],6)


if __name__ == '__main__':
    gr_unittest.run(qa_framer_cc, "qa_framer_cc.xml")
