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

from __future__ import print_function, division
import numpy as np
import unittest
from pygfdm.mapping import map_to_waveform_resource_grid, get_subcarrier_map

import gfdmlib as vc


class MapperTests(unittest.TestCase):
    def setUp(self):
        self.params = vc.defaultGFDM.get_defaultGFDM('BER')
        self.params.Non = self.params.Kon * self.params.Mon
        print(self.params.__dict__)
        self.subcarrier_map = np.arange(self.params.Kon)

    def tearDown(self):
        self.params = None
        self.subcarrier_map = None

    def test_001_map_full(self):
        # mod = vc.Modulator.DefaultModulator(params)
        mapper = vc.mapping.Mapper(self.params)

        d = np.arange(self.params.Non, dtype=np.complex64) + 1
        f = mapper.doMap(d)
        ref = map_to_waveform_resource_grid(d, self.params.Kon, self.params.K,
                                            self.subcarrier_map, True)
        self.assertTrue(np.all(f == ref))

    def test_002_map_active_subcarriers(self):
        self.params.Kon = 112
        self.params.Non = self.params.Kon * self.params.Mon

        self.subcarrier_map = np.arange(4, self.params.Kon + 4)
        self.params.Kset = self.subcarrier_map

        mapper = vc.mapping.Mapper(self.params)

        d = np.arange(self.params.Non, dtype=np.complex64) + 1
        f = mapper.doMap(d)
        ref = map_to_waveform_resource_grid(d, self.params.Kon, self.params.K,
                                            self.subcarrier_map, True)
        self.assertTrue(np.all(f == ref))

    def test_003_map_complex_config(self):
        self.params.K = 64
        self.params.Kon = 52
        self.params.Non = self.params.Kon * self.params.Mon
        self.params.N = self.params.K * self.params.M

        print(self.params.__dict__)

        self.subcarrier_map = get_subcarrier_map(self.params.K,
                                                 self.params.Kon, dc_free=True)
        self.params.Kset = self.subcarrier_map

        mapper = vc.mapping.Mapper(self.params)

        d = np.arange(self.params.Non, dtype=np.complex64) + 1
        f = mapper.doMap(d)
        ref = map_to_waveform_resource_grid(d, self.params.Kon, self.params.K,
                                            self.subcarrier_map, True)

        self.assertTrue(np.all(f == ref))


if __name__ == '__main__':
    unittest.main(failfast=True)
