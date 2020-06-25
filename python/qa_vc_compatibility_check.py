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

import numpy as np
import unittest
import itertools
from pygfdm.mapping import map_to_waveform_resource_grid, get_subcarrier_map
from pygfdm.filters import gfdm_filter_taps, get_frequency_domain_filter
from pygfdm.gfdm_modulation import gfdm_modulate_block

import gfdmlib as vc


class MapperTests(unittest.TestCase):
    def setUp(self):
        self.params = vc.get_defaultGFDM('BER')
        self.params.Non = self.params.Kon * self.params.Mon
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

        self.subcarrier_map = get_subcarrier_map(self.params.K,
                                                 self.params.Kon, dc_free=True)
        self.params.Kset = self.subcarrier_map

        mapper = vc.mapping.Mapper(self.params)

        d = np.arange(self.params.Non, dtype=np.complex64) + 1
        f = mapper.doMap(d)
        ref = map_to_waveform_resource_grid(d, self.params.Kon, self.params.K,
                                            self.subcarrier_map, True)

        self.assertTrue(np.all(f == ref))


class ModulatorTests(unittest.TestCase):
    def setUp(self):
        self.params = vc.get_defaultGFDM('BER')
        self.params.Non = self.params.Kon * self.params.Mon
        self.subcarrier_map = np.arange(self.params.Kon)

    def tearDown(self):
        self.params = None
        self.subcarrier_map = None

    def test_001_td_filter(self):
        self.params.K = 64
        self.params.Kon = 52
        self.params.Non = self.params.Kon * self.params.Mon
        self.params.N = self.params.K * self.params.M
        taps = vc.gfdmutil.get_transmitter_pulse(self.params)
        my_taps = gfdm_filter_taps(self.params.pulse, self.params.a,
                                   self.params.M, self.params.K, 1)
        my_taps = np.fft.fftshift(my_taps)
        my_taps *= taps[0] / my_taps[0]
        self.assertTrue(np.all(np.abs(my_taps - taps) < 1e-12))

    def test_002_fd_filter(self):
        self.params.K = 64
        self.params.Kon = 52
        self.params.Non = self.params.Kon * self.params.Mon
        self.params.N = self.params.K * self.params.M
        self.params.pulse = 'rrc'
        taps = vc.gfdmutil.get_transmitter_pulse(self.params)

        g2 = taps[::self.params.K // self.params.L]
        G2 = np.fft.fft(g2)

        freq_taps = get_frequency_domain_filter('rrc',
                                                self.params.a, self.params.M,
                                                self.params.K, self.params.L)
        freq_taps *= np.abs(G2[0]) / np.abs(freq_taps[0])
        self.assertTrue(np.all(np.abs(freq_taps - G2) < 1e-3))

    def test_003_rc_fd_filter(self):
        self.params.K = 64
        self.params.Kon = 52
        self.params.Non = self.params.Kon * self.params.Mon
        self.params.N = self.params.K * self.params.M
        self.params.pulse = 'rc'
        taps = vc.gfdmutil.get_transmitter_pulse(self.params)

        g2 = taps[::self.params.K // self.params.L]
        G2 = np.fft.fft(g2)

        freq_taps = get_frequency_domain_filter('rc',
                                                self.params.a, self.params.M,
                                                self.params.K, self.params.L)
        freq_taps *= np.abs(G2[0]) / np.abs(freq_taps[0])

        # plot_taps(G2.real, freq_taps.real)
        # print(np.max(np.abs(G2 - freq_taps)))
        self.assertTrue(np.all(np.abs(freq_taps - G2) < 1e-3))

    def test_004_std_config(self):
        self.validate_parameter_set(64, 52, 5, 5e-14)

    def test_005_vc_video(self):
        self.validate_parameter_set(64, 52, 9, 5e-14)

    def test_006_vc_ll(self):
        self.validate_parameter_set(32, 28, 3, 5e-14)

    def test_007_parameter_list(self):
        Ms = np.array([3, 5, 6, 9, 12, 15])
        Ks = np.array([8, 16, 32, 64, 128, 256, 512])
        myset = itertools.product(Ks, Ms)
        for K, M in myset:
            Kon = int(0.8125 * K)
            Kon = 2 * (Kon // 2)
            print('K={}\tKon={}\tM={}'.format(K, Kon, M))
            self.validate_parameter_set(K, Kon, M, 5e-3)

    def validate_parameter_set(self, K, Kon, M, tolerance=5e-4):
        self.params.K = K
        self.params.Kon = Kon
        self.params.M = self.params.Mon = M
        self.params.Non = self.params.Kon * self.params.Mon
        self.params.N = self.params.K * self.params.M
        self.params.pulse = 'rc_fd'

        self.subcarrier_map = get_subcarrier_map(self.params.K,
                                                 self.params.Kon, dc_free=True)
        self.params.Kset = self.subcarrier_map
        # vtaps = vc.gfdmutil.get_transmitter_pulse(self.params)
        # vtaps = vtaps[::self.params.K // self.params.L]
        # taps = np.fft.fft(vtaps)
        # taps = get_frequency_domain_filter(self.params.pulse, self.params.a,
        #                                    self.params.M, self.params.K,
        #                                    self.params.L)
        taps = vc.gfdmutil.get_transmitter_pulse(self.params)
        g2 = taps[::self.params.K // self.params.L]
        taps = np.fft.fft(g2)

        mod = vc.DefaultModulator(self.params)

        d = np.random.randn(self.params.Non) + 1.j * np.random.randn(self.params.Non)
        dframe = map_to_waveform_resource_grid(d, self.params.Kon,
                                               self.params.K,
                                               self.subcarrier_map, True)

        vdata = mod.modulate(dframe)
        vdata *= 1. / np.linalg.norm(vdata)
        gdata = gfdm_modulate_block(dframe.T, taps, self.params.M,
                                    self.params.K, self.params.L, False)
        gdata *= 1. / np.linalg.norm(gdata)
        print(np.max(np.abs(vdata - gdata)))
        self.assertTrue(np.all(np.abs(vdata - gdata) < tolerance))


def plot_taps(ataps, btaps):
    import matplotlib.pyplot as plt
    plt.plot(ataps)
    plt.plot(btaps)
    plt.show()


def main():
    import matplotlib.pyplot as plt
    params = vc.defaultGFDM.get_defaultGFDM('BER')
    params.Non = params.Kon * params.Mon
    subcarrier_map = np.arange(params.Kon)
    params.K = 64
    params.Kon = 52
    params.Non = params.Kon * params.Mon
    params.N = params.K * params.M
    print(params.__dict__)
    taps = vc.gfdmutil.get_transmitter_pulse(params)
    my_taps = gfdm_filter_taps(params.pulse, params.a, params.M, params.K, 1)
    my_taps = np.fft.fftshift(my_taps)
    my_taps *= taps[0] / my_taps[0]

    print(taps.dtype, my_taps.dtype)
    print(np.max(np.abs(taps - my_taps)))

    print(np.all(np.abs(my_taps - taps) < 1e-12))

    plt.plot(taps)
    plt.plot(my_taps)
    plt.show()


if __name__ == '__main__':
    # main()
    unittest.main(failfast=True)
