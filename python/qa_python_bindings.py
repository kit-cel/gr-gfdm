#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
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

from gnuradio import gr_unittest
import numpy as np
from pygfdm.mapping import map_to_waveform_resources, get_subcarrier_map
from pygfdm.mapping import get_data_matrix
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.gfdm_modulation import gfdm_modulate_block
from pygfdm.gfdm_receiver import gfdm_demodulate_block
from pygfdm.utils import get_random_qpsk
from pygfdm.cyclic_prefix import add_cyclic_prefix, pinch_block, get_raised_cosine_ramp, get_window_len, add_cyclic_starfix

from gfdm_python import Resource_mapper, Modulator, Demodulator, Cyclic_prefixer


class BindingTests(gr_unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_001_map_timeslots(self):
        timeslots = 15
        subcarriers = 32
        active_subcarriers = 24
        subcarrier_map = np.arange(4, 28, dtype=np.int)
        mapper = Resource_mapper(timeslots, subcarriers,
                                 active_subcarriers,
                                 subcarrier_map,
                                 True)

        self.assertEqual(mapper.block_size(),
                         timeslots * active_subcarriers)
        self.assertEqual(mapper.frame_size(),
                         timeslots * subcarriers)

        d = np.arange(timeslots * active_subcarriers,
                      dtype=np.complex64) + 1
        f = mapper.map_to_resources(d)
        ref = map_to_waveform_resources(d, active_subcarriers,
                                        subcarriers,
                                        subcarrier_map,
                                        True)
        self.assertComplexTuplesAlmostEqual(f, ref)


    def test_002_map_subcarriers(self):
        timeslots = 15
        subcarriers = 32
        active_subcarriers = 24
        subcarrier_map = np.arange(4, 28, dtype=np.int)
        mapper = Resource_mapper(timeslots, subcarriers,
                                 active_subcarriers,
                                 subcarrier_map,
                                 False)

        self.assertEqual(mapper.block_size(),
                         timeslots * active_subcarriers)
        self.assertEqual(mapper.frame_size(),
                         timeslots * subcarriers)

        d = np.arange(timeslots * active_subcarriers,
                      dtype=np.complex64) + 1
        f = mapper.map_to_resources(d)
        ref = map_to_waveform_resources(d, active_subcarriers,
                                        subcarriers,
                                        subcarrier_map,
                                        False)
        self.assertComplexTuplesAlmostEqual(f, ref)

    def test_003_demap_timeslots(self):
        timeslots = 15
        subcarriers = 32
        active_subcarriers = 24
        subcarrier_map = get_subcarrier_map(subcarriers,
                                            active_subcarriers,
                                            True)
        # subcarrier_map = np.arange(4, 28, dtype=np.int)
        mapper = Resource_mapper(timeslots, subcarriers,
                                 active_subcarriers,
                                 subcarrier_map,
                                 True)

        self.assertEqual(mapper.block_size(),
                         timeslots * active_subcarriers)
        self.assertEqual(mapper.frame_size(),
                         timeslots * subcarriers)

        d = np.arange(timeslots * active_subcarriers,
                      dtype=np.complex64) + 1

        ref = map_to_waveform_resources(d, active_subcarriers,
                                        subcarriers,
                                        subcarrier_map,
                                        True)
        f = mapper.demap_from_resources(ref)
        self.assertComplexTuplesAlmostEqual(f, d)

    def test_004_demap_subcarriers(self):
        timeslots = 15
        subcarriers = 32
        active_subcarriers = 24
        subcarrier_map = np.arange(4, 28, dtype=np.int)
        mapper = Resource_mapper(timeslots, subcarriers,
                                 active_subcarriers,
                                 subcarrier_map,
                                 False)

        self.assertEqual(mapper.block_size(),
                         timeslots * active_subcarriers)
        self.assertEqual(mapper.frame_size(),
                         timeslots * subcarriers)

        d = np.arange(timeslots * active_subcarriers,
                      dtype=np.complex64) + 1

        ref = map_to_waveform_resources(d, active_subcarriers,
                                        subcarriers,
                                        subcarrier_map,
                                        False)
        f = mapper.demap_from_resources(ref)
        self.assertComplexTuplesAlmostEqual(f, d)


class PrefixerTests(gr_unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_001_prefix(self):
        timeslots = 19
        subcarriers = 32

        block_len = timeslots * subcarriers
        cp_len = 16
        ramp_len = 4
        cs_len = ramp_len * 2
        window_len = get_window_len(cp_len, timeslots, subcarriers,
                                    cs_len)
        window_taps = get_raised_cosine_ramp(ramp_len, window_len)
        data = np.arange(block_len, dtype=np.complex) + 1
        ref = add_cyclic_starfix(data, cp_len, cs_len)
        ref = pinch_block(ref, window_taps)
        ref = ref.astype(np.complex64)

        prefixer = Cyclic_prefixer(block_len, cp_len, cs_len, ramp_len,
                                   window_taps)
        res = prefixer.add_cyclic_prefix(data)

        self.assertComplexTuplesAlmostEqual(res, ref, 6)

    def test_002_prefix_shifted(self):
        print('shifted')
        timeslots = 3
        subcarriers = 32
        cyclic_shift = 4

        block_len = timeslots * subcarriers
        cp_len = 16
        cs_len = cp_len // 2
        ramp_len = 4
        window_len = get_window_len(cp_len, timeslots, subcarriers,
                                    cs_len)
        window_taps = get_raised_cosine_ramp(ramp_len, window_len)
        data = np.arange(block_len, dtype=np.complex) + 1
        ref = add_cyclic_starfix(data, cp_len, cs_len)
        ref = np.concatenate((data[-(cp_len - cyclic_shift):], data, data[0:cs_len + cyclic_shift]))
        ref = pinch_block(ref, window_taps)
        ref = ref.astype(np.complex64)

        prefixer = Cyclic_prefixer(block_len, cp_len, cs_len, ramp_len,
                                   window_taps, cyclic_shift)
        res = prefixer.add_cyclic_prefix(data)
        self.assertEqual(res.size, cp_len + block_len + cs_len)
        self.assertComplexTuplesAlmostEqual(res, ref, 5)


class ModulatorTests(gr_unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_001_small(self):
        timeslots = 16
        subcarriers = 4
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)

        ref = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                  overlap, False)


        mod = Modulator(timeslots, subcarriers, overlap, taps)
        res = mod.modulate(data)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)

    def test_002_big(self):
        timeslots = 21
        subcarriers = 128
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)

        ref = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                  overlap, False)


        mod = Modulator(timeslots, subcarriers, overlap, taps)
        res = mod.modulate(data)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)


class DemodulatorTests(gr_unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_001_init(self):
        timeslots = 25
        subcarriers = 96
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        self.assertEqual(timeslots, demod.timeslots())
        self.assertEqual(subcarriers, demod.subcarriers())
        self.assertEqual(overlap, demod.overlap())
        self.assertEqual(timeslots * subcarriers, demod.block_size())
        self.assertComplexTuplesAlmostEqual(taps, demod.filter_taps())

    def test_002_small(self):
        timeslots = 16
        subcarriers = 4
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)
        frame = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                    overlap, False)
        ref = gfdm_demodulate_block(frame, taps, subcarriers, timeslots, overlap)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        res = demod.demodulate(frame)

        self.assertComplexTuplesAlmostEqual(ref, res, 6)

    def test_003_big(self):
        timeslots = 21
        subcarriers = 128
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)
        frame = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                    overlap, False)
        ref = gfdm_demodulate_block(frame, taps, subcarriers, timeslots, overlap)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        res = demod.demodulate(frame)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)

    def test_004_big_equalize(self):
        timeslots = 21
        subcarriers = 128
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)
        frame = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                    overlap, False)
        ref = gfdm_demodulate_block(frame, taps, subcarriers, timeslots, overlap)
        eq_vals = np.ones(ref.size, ref.dtype) * np.exp(1.j)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        res = demod.demodulate_equalize(frame * np.exp(1.j), eq_vals)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)

    def test_005_steps(self):
        timeslots = 5
        subcarriers = 32
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)
        frame = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                    overlap, False)
        ref = gfdm_demodulate_block(frame, taps, subcarriers, timeslots, overlap)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        fd_res = demod.fft_filter_downsample(frame)
        res = demod.transform_subcarriers_to_td(fd_res)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)

        for _ in range(2):
            ic_res = demod.cancel_sc_interference(data, fd_res)
            res = demod.transform_subcarriers_to_td(ic_res)

        self.assertComplexTuplesAlmostEqual(data, res, 1)

    def test_006_steps_equalize(self):
        timeslots = 5
        subcarriers = 32
        overlap = 2
        filter_alpha = 0.35

        taps = get_frequency_domain_filter('rrc', filter_alpha,
                                           timeslots, subcarriers,
                                           overlap)

        data = get_random_qpsk(timeslots * subcarriers)
        D = get_data_matrix(data, subcarriers, group_by_subcarrier=False)
        frame = gfdm_modulate_block(D, taps, timeslots, subcarriers,
                                    overlap, False)
        ref = gfdm_demodulate_block(frame, taps, subcarriers, timeslots, overlap)
        eq_vals = np.ones(ref.size, ref.dtype) * np.exp(1.j)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        fd_res = demod.fft_equalize_filter_downsample(frame * np.exp(1.j), eq_vals)
        res = demod.transform_subcarriers_to_td(fd_res)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)




if __name__ == '__main__':
    gr_unittest.run(BindingTests)
    gr_unittest.run(ModulatorTests)
    gr_unittest.run(PrefixerTests)
    gr_unittest.run(DemodulatorTests)
