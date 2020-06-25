#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2017, 2020 Johannes Demel.
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
import numpy as np
from pygfdm.mapping import get_subcarrier_map
from pygfdm.preamble import mapped_preamble


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


class qa_channel_estimator_cc(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()
        self.filtertype = 'rrc'
        self.filteralpha = .5
        self.seed = int(3660365253)

    def tearDown(self):
        self.tb = None

    def test_001_simple(self):
        timeslots = 3
        subcarriers = 32
        active_subcarriers = 24
        overlap = 2
        cp_len = subcarriers // 2
        ramp_len = cp_len // 2

        subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers,
                                            dc_free=True)
        preambles = mapped_preamble(self.seed, self.filtertype, self.filteralpha, active_subcarriers,
                                    subcarriers, subcarrier_map, overlap, cp_len, ramp_len)
        core_preamble = preambles[1]

        dut = gfdm.channel_estimator_cc(
            timeslots, subcarriers, active_subcarriers, True, 1, core_preamble)
        src = blocks.vector_source_c(core_preamble)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, dut, snk)
        self.tb.run()

        res = np.array(snk.data())
        self.assertComplexTuplesAlmostEqual(
            res, np.ones(res.size, dtype=res.dtype), 6)

    def test_002_selective(self):
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

        dut = gfdm.channel_estimator_cc(
            timeslots, subcarriers, active_subcarriers, True, 1, core_preamble)
        src = blocks.vector_source_c(data)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, dut, snk)
        self.tb.run()

        res = np.array(snk.data())
        lowres = res[0:active_symbols // 2]
        hires = res[-active_symbols // 2:]

        fh = np.fft.fft(h, timeslots * subcarriers)
        lowfh = fh[0:active_symbols // 2]
        hifh = fh[-active_symbols // 2:]

        self.assertComplexTuplesAlmostEqual(lowres, lowfh, 1)
        self.assertComplexTuplesAlmostEqual(hires, hifh, 1)

    def test_003_snr(self):
        nframes = 30
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
        snrs = np.arange(3, 3 * nframes, 3, dtype=np.float)
        snrs_lin = 10. ** (snrs / 10.)
        expected_snrs_lin = np.concatenate(((np.inf,), snrs_lin))

        for i, snr_lin in enumerate(snrs_lin):
            nscale = calculate_noise_scale(
                snr_lin, sigenergy, active_ratio, core_preamble.size)
            noise = get_noise_vector(core_preamble.size, nscale)

            d = core_preamble + noise
            data = np.concatenate((data, d))

        dut = gfdm.channel_estimator_cc(
            timeslots, subcarriers, active_subcarriers, True, 1, core_preamble)
        src = blocks.vector_source_c(data)
        snk = blocks.vector_sink_c()
        self.tb.connect(src, dut, snk)
        self.tb.run()

        res = np.array(snk.data())
        self.assertEqual(res.size, nframes * timeslots * subcarriers)

        tags = snk.tags()
        snr_tags = [t for t in tags if pmt.eq(t.key, pmt.mp("snr_lin"))]

        for i, t in enumerate(snr_tags):
            self.assertEqual(t.offset, i * timeslots * subcarriers)
            res_lin = pmt.to_float(t.value)
            res_db = 10. * np.log10(res_lin)
            ref_db = 10. * np.log10(expected_snrs_lin[i])
            # print(f"Reference: {ref_db:6.3f}dB\t{res_db:6.3f}dB")

            if np.isfinite(ref_db):
                self.assertTrue(np.abs(res_db - ref_db) < 1.)


if __name__ == '__main__':
    gr_unittest.run(qa_channel_estimator_cc)
