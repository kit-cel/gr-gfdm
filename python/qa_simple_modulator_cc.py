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
from gi.overrides.Gtk import Label
from gnuradio import gr, gr_unittest
from gnuradio import blocks
import gfdm_swig as gfdm
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.gfdm_modulation import gfdm_modulate_block, gfdm_transform_subcarriers_to_fd,\
    gfdm_upsample_subcarriers_in_fd, gfdm_filter_subcarriers_in_fd, gfdm_combine_subcarriers_in_fd
from pygfdm.modulation import get_data_matrix
import numpy as np


class qa_simple_modulator_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_t(self):
        alpha = .5
        M = 8
        K = 4
        L = 2
        taps = get_frequency_domain_filter('rrc', alpha, M, K, L)
        data = np.repeat(np.arange(1, K + 1), M)
        D = get_data_matrix(data, K, group_by_subcarrier=False)
        print data
        print D
        print np.reshape(data, (-1, M)).T
        src = blocks.vector_source_c(data)
        mod = gfdm.simple_modulator_cc(M, K, L, taps)
        dst = blocks.vector_sink_c()

        self.tb.connect(src, mod, dst)
        # set up fg
        self.tb.run()
        # check data
        print 'NOW: Check results!'
        res = np.array(dst.data())

        F = gfdm_transform_subcarriers_to_fd(D, M)
        F = gfdm_upsample_subcarriers_in_fd(F, K, L)
        F = gfdm_filter_subcarriers_in_fd(F, taps, M, K, L)
        print F[:, 0:2]
        # print np.reshape(res, (-1, M * L)).T

        X = gfdm_combine_subcarriers_in_fd(F, M, K, L)
        print X[0:M]
        print res[0:M]
        print
        print X
        print res
        # print np.reshape(res, (-1, M)).T

        # ref = gfdm_modulate_block(D, taps, M, K, L)
        # print np.reshape(ref, (-1, M)).T
        # self.assertComplexTuplesAlmostEqual(F[:, 0:2].T.flatten(), res, 5)

if __name__ == '__main__':
    # gr_unittest.run(qa_simple_modulator_cc, "qa_simple_modulator_cc.xml")
    gr_unittest.run(qa_simple_modulator_cc)
