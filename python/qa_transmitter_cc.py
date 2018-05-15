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
import gfdm_swig as gfdm
from pygfdm.cyclic_prefix import get_window_len, get_raised_cosine_ramp
from pygfdm.cyclic_prefix import add_cyclic_starfix, pinch_block
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.utils import get_random_qpsk
from pygfdm.preamble import get_sync_symbol
from pygfdm.mapping import get_subcarrier_map, map_to_waveform_resources
from pygfdm.mapping import get_data_matrix
from pygfdm.gfdm_modulation import gfdm_modulate_block
import numpy as np


class qa_transmitter_cc (gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        np.set_printoptions(precision=2)
        n_frames = 7
        alpha = .5
        active_subcarriers = 52
        timeslots = 9
        subcarriers = 64
        overlap = 2
        cp_len = 8
        cs_len = 4
        ramp_len = 4

        window_len = get_window_len(cp_len, timeslots, subcarriers, cs_len)
        taps = get_frequency_domain_filter('rrc', alpha, timeslots,
                                           subcarriers, overlap)
        # taps /= np.sqrt(calculate_signal_energy(taps) / time)
        window_taps = get_raised_cosine_ramp(ramp_len, window_len)
        pn_symbols = get_random_qpsk(subcarriers)
        H_preamble = get_frequency_domain_filter('rrc', alpha, 2,
                                                 subcarriers, overlap)
        preamble = get_sync_symbol(pn_symbols, H_preamble, subcarriers,
                                   overlap, cp_len, ramp_len)[0]
        smap = get_subcarrier_map(subcarriers, active_subcarriers, True)

        ref = np.array([], dtype=np.complex)
        data = np.array([], dtype=np.complex)

        for i in range(n_frames):
            d = get_random_qpsk(active_subcarriers * timeslots)
            dd = map_to_waveform_resources(d, active_subcarriers,
                                           subcarriers, smap)
            D = get_data_matrix(dd, subcarriers, group_by_subcarrier=False)
            b = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                    overlap, False)
            b = add_cyclic_starfix(b, cp_len, cs_len)
            b = pinch_block(b, window_taps)
            ref = np.concatenate((ref, preamble, b))
            data = np.concatenate((data, d))

        src = blocks.vector_source_c(data)
        dut = gfdm.transmitter_cc(timeslots, subcarriers, active_subcarriers,
                                  cp_len, cs_len, ramp_len, smap, True,
                                  overlap, taps, window_taps, preamble)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, dut, dst)
        self.tb.run()
        res = np.array(dst.data())[0:len(ref)]
        self.assertComplexTuplesAlmostEqual(ref, res, 5)


if __name__ == '__main__':
    gr_unittest.run(qa_transmitter_cc)
