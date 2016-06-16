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
import gfdm_swig as gfdms
from pygfdm.gfdm_modulation import gfdm_modulate_block
from pygfdm.mapping import get_data_matrix
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.utils import get_random_qpsk
import numpy as np


class qa_modulator_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        nsubcarrier = 4
        ntimeslots = 16
        filter_alpha = 0.35
        tag_key = "frame_len"
        fft_length = nsubcarrier * ntimeslots
        taps = get_frequency_domain_filter('rrc', filter_alpha, ntimeslots, nsubcarrier, 2)

        data = get_random_qpsk(nsubcarrier * ntimeslots)
        # data = np.zeros(nsubcarrier)
        # data[0] = 1
        # data = np.repeat(data, ntimeslots)
        D = get_data_matrix(data, nsubcarrier, group_by_subcarrier=False)
        print D

        md = gfdms.modulator_cc(nsubcarrier, ntimeslots, filter_alpha, fft_length, 1, tag_key)
        tagger = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 1, fft_length, tag_key)
        src = blocks.vector_source_c(data)
        dst = blocks.vector_sink_c()
        self.tb.connect(src, tagger, md, dst)
        self.tb.run()

        res = np.array(dst.data())
        print res
        ref = gfdm_modulate_block(D, taps, ntimeslots, nsubcarrier, 2, True)
        print ref

        self.assertComplexTuplesAlmostEqual(ref, res, 2)


if __name__ == '__main__':
    # gr_unittest.run(qa_modulator_cc, "qa_modulator_cc.xml")
    gr_unittest.run(qa_modulator_cc)
