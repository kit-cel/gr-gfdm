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
from pygfdm.preamble import mapped_preamble
from pygfdm.symbolmapping import bits2symbols, symbols2bits, generate_constellation

from gfdm_python import Resource_mapper, Modulator, Demodulator, Cyclic_prefixer, Preamble_channel_estimator


def calculate_energy(vec):
    return np.sum(calculate_element_energy(vec))


def calculate_element_energy(vec):
    return vec.real ** 2 + vec.imag ** 2


def get_noise_vector(size, scale):
    noise = np.random.randn(size) + 1.j * np.random.randn(size)
    noise /= np.abs(noise)
    return noise * scale


def calculate_noise_scale(snr_lin, signalenergy,
                          activecarrier_ratio, noise_vector_length):
    nscale = 1. / np.sqrt(snr_lin)
    nscale *= np.sqrt(activecarrier_ratio * 2. *
                      signalenergy / noise_vector_length)
    return nscale


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
        ref = np.concatenate(
            (data[-(cp_len + cyclic_shift):], data, data[0:cs_len - cyclic_shift]))
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
        ref = gfdm_demodulate_block(
            frame, taps, subcarriers, timeslots, overlap)

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
        ref = gfdm_demodulate_block(
            frame, taps, subcarriers, timeslots, overlap)

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
        ref = gfdm_demodulate_block(
            frame, taps, subcarriers, timeslots, overlap)
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
        ref = gfdm_demodulate_block(
            frame, taps, subcarriers, timeslots, overlap)

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
        ref = gfdm_demodulate_block(
            frame, taps, subcarriers, timeslots, overlap)
        eq_vals = np.ones(ref.size, ref.dtype) * np.exp(1.j)

        demod = Demodulator(timeslots, subcarriers, overlap, taps)
        fd_res = demod.fft_equalize_filter_downsample(
            frame * np.exp(1.j), eq_vals)
        res = demod.transform_subcarriers_to_td(fd_res)

        self.assertComplexTuplesAlmostEqual(ref, res, 5)


class EstimatorTests(gr_unittest.TestCase):
    def setUp(self):
        self.filtertype = 'rrc'
        self.filteralpha = .5
        self.seed = int(3660365253)

    def tearDown(self):
        pass

    def test_001_selective(self):
        timeslots = 5
        subcarriers = 64
        active_subcarriers = 52
        overlap = 2
        cp_len = subcarriers // 2
        ramp_len = cp_len // 2
        active_symbols = timeslots * active_subcarriers

        subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers,
                                            dc_free=True)
        preambles = mapped_preamble(self.seed, self.filtertype, self.filteralpha,
                                    active_subcarriers, subcarriers,
                                    subcarrier_map, overlap, cp_len, ramp_len)
        full_preamble = preambles[0]
        core_preamble = preambles[1]
        h = np.array([1., .5, .1j, .1+.05j], dtype=np.complex)
        data = np.convolve(full_preamble, h, 'full')[0:full_preamble.size]
        data = data[cp_len:-ramp_len]
        self.assertEqual(data.size, core_preamble.size)

        estimator = Preamble_channel_estimator(
            timeslots, subcarriers, active_subcarriers, True, 1, core_preamble)
        self.assertEqual(estimator.timeslots(), timeslots)
        self.assertEqual(estimator.subcarriers(), subcarriers)
        self.assertEqual(estimator.active_subcarriers(), active_subcarriers)
        self.assertEqual(estimator.frame_len(), timeslots * subcarriers)
        self.assertEqual(estimator.is_dc_free(), True)

        res = estimator.estimate_frame(data)
        lowres = res[0:active_symbols // 2]
        hires = res[-active_symbols // 2:]

        fh = np.fft.fft(h, timeslots * subcarriers)
        lowfh = fh[0:active_symbols // 2]
        hifh = fh[-active_symbols // 2:]

        self.assertComplexTuplesAlmostEqual(lowres, lowfh, 1)
        self.assertComplexTuplesAlmostEqual(hires, hifh, 1)

    def test_002_snr(self):
        timeslots = 5
        subcarriers = 1024
        active_subcarriers = 936
        overlap = 2
        cp_len = subcarriers // 2
        ramp_len = cp_len // 2
        active_ratio = subcarriers / active_subcarriers

        subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers,
                                            dc_free=True)
        preambles = mapped_preamble(self.seed, self.filtertype, self.filteralpha,
                                    active_subcarriers, subcarriers,
                                    subcarrier_map, overlap, cp_len, ramp_len)

        core_preamble = preambles[1]

        sigenergy = calculate_energy(core_preamble)

        data = np.copy(core_preamble)
        snr = 4.0
        snr_lin = 10. ** (snr / 10.)

        nscale = calculate_noise_scale(
            snr_lin, sigenergy, active_ratio, core_preamble.size)
        noise = get_noise_vector(core_preamble.size, nscale)

        data = core_preamble + noise

        estimator = Preamble_channel_estimator(
            timeslots, subcarriers, active_subcarriers, True, 1, core_preamble)

        res = estimator.estimate_snr(data)
        res_db = 10. * np.log10(res)

        print(res, snr_lin)
        print(res_db, snr)
        self.assertTrue(np.abs(res_db - snr) < 1.)


class CyclicDelayDiversityTests(gr_unittest.TestCase):
    def setUp(self):
        self.filtertype = 'rrc'
        self.filteralpha = .5
        self.seed = int(3660365253)
        self.overlap = 2
        self.ic_iterations = 2
        self.constellation_order = 2

        self.channels = [
            np.array([1., .3+.1j, .4j]),
            np.array([1.j, .4+.1j, .2]),
            np.array([.7j, .8+.1j, .1]),
            np.array([.7+.7j, .1+.3j, .6+.2j]),
        ]
        for i, chan in enumerate(self.channels):
            ec = np.mean(np.abs(chan) ** 2)
            chan /= np.sqrt(ec)
            self.channels[i] = chan

    def tearDown(self):
        pass

    def simulate_channel(self, parallel_frames):
        rx_frames = [np.convolve(frame, self.channels[i], 'full')[0:frame.size]
                     for i, frame in enumerate(parallel_frames)]
        return np.sum(rx_frames, axis=0)

    def get_effective_channel_taps(self, cyclic_shift):
        nparallel = len(cyclic_shift)
        max_shift = np.max(cyclic_shift)
        max_channel_delay = np.max(
            [self.channels[i].size for i in range(nparallel)])
        grid = np.zeros((nparallel, max_shift + max_channel_delay),
                        dtype=self.channels[0].dtype)
        for i, cs in enumerate(cyclic_shift):
            chan = self.channels[i]
            print(chan)
            grid[i, cs:cs + chan.size] += chan
        return np.sum(grid, axis=0)

    def test_001_selective(self):
        timeslots = 5
        subcarriers = 64
        active_subcarriers = 52
        cp_len = subcarriers // 2
        ramp_len = cp_len // 2
        # cyclic_shift = [0, 2, 4, 6, ]
        cyclic_shift = [0, 4, 8, 12, ]
        active_symbols = timeslots * active_subcarriers
        preambles = [self.generate_preamble(
            subcarriers, active_subcarriers, cp_len, ramp_len, cs) for cs in cyclic_shift]

        full_preambles = [p[0] for p in preambles]
        core_preambles = [p[1] for p in preambles]
        core_preamble = core_preambles[0]

        rx_preamble = self.simulate_channel(full_preambles)
        # print(rx_preamble)
        print(rx_preamble.size, full_preambles[0].size)
        rx_core_preamble = rx_preamble[cp_len:-ramp_len]

        estimator = Preamble_channel_estimator(
            timeslots, subcarriers, active_subcarriers, True, 1,
            core_preamble)

        res = estimator.estimate_frame(rx_core_preamble)
        lowres = res[0:active_symbols // 2]
        hires = res[-active_symbols // 2:]

        effective_channel = self.get_effective_channel_taps(cyclic_shift)
        print(effective_channel)
        fh = np.fft.fft(effective_channel, timeslots * subcarriers)
        lowfh = fh[0:active_symbols // 2]
        hifh = fh[-active_symbols // 2:]
        import matplotlib.pyplot as plt
        plt.plot(res.real)
        plt.plot(res.imag)
        plt.plot(fh.real, ls='dashed')
        plt.plot(fh.imag, ls='dashed')
        # plt.plot(np.abs(lowres))
        # plt.plot(np.angle(lowres))
        # plt.plot(np.abs(lowfh), ls='dashed')
        # plt.plot(np.angle(lowfh), ls='dashed')
        plt.show()

        self.assertFloatTuplesAlmostEqual(np.unwrap(np.angle(lowres)),
                                          np.unwrap(np.angle(lowfh)), 0)
        self.assertFloatTuplesAlmostEqual(np.unwrap(np.angle(hires)),
                                          np.unwrap(np.angle(hifh)), 0)
        # self.assertFloatTuplesAlmostEqual(np.abs(hires),
        #                                   np.abs(hifh), 0)
        # self.assertComplexTuplesAlmostEqual(lowres, lowfh, 1)
        self.assertComplexTuplesAlmostEqual(hires, hifh, 0)

    def generate_preamble(self, subcarriers, active_subcarriers, cp_len,
                          ramp_len, cyclic_shift):
        subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers,
                                            dc_free=True)
        return mapped_preamble(self.seed, self.filtertype, self.filteralpha,
                               active_subcarriers, subcarriers,
                               subcarrier_map, self.overlap, cp_len, ramp_len,
                               use_zadoff_chu=True, cyclic_shift=cyclic_shift)


if __name__ == '__main__':
    gr_unittest.run(BindingTests)
    gr_unittest.run(ModulatorTests)
    gr_unittest.run(PrefixerTests)
    gr_unittest.run(DemodulatorTests)
    gr_unittest.run(EstimatorTests)
    gr_unittest.run(CyclicDelayDiversityTests)
