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


import numpy as np
from gfdm.pygfdm.mapping import get_subcarrier_map
from gfdm.pygfdm.preamble import mapped_preamble


def lin2db(v):
    return 10. * np.log10(v)


def db2lin(v):
    return 10. ** (v / 10.)


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


def estimate_snr0(rxvec, subcarrier_map, subcarriers):
    fd = np.fft.fft(rxvec, 2 * subcarriers)
    fds = fd[0::2]
    fdn = fd[1::2]

    # se = calculate_energy(fds)

    se = calculate_energy(fds[subcarrier_map])
    ne = calculate_energy(fdn[subcarrier_map])
    return (se - ne) / ne


def estimate_snr1(rxvec, subcarrier_map, subcarriers, preamble):
    fd = np.fft.fft(rxvec, 2 * subcarriers)
    fds = fd[0::2][subcarrier_map]
    # fdn = fd[1::2]
    fp = np.fft.fft(preamble, 2 * subcarriers)
    fps = fp[0::2][subcarrier_map]
    se = calculate_energy(fds)
    ne = calculate_energy(fds - fps)
    return (se - ne) / ne


def foo(nframes=15,
        timeslots=5,
        subcarriers=1024,
        active_subcarriers=936,
        overlap=2,
        filtertype='rrc',
        filteralpha=.5,
        seed=int(3660365253)):

    cp_len = subcarriers // 2
    ramp_len = cp_len // 2
    active_ratio = subcarriers / active_subcarriers

    subcarrier_map = get_subcarrier_map(subcarriers, active_subcarriers,
                                        dc_free=True)
    preambles = mapped_preamble(seed, filtertype, filteralpha,
                                active_subcarriers, subcarriers,
                                subcarrier_map, overlap, cp_len, ramp_len)
    core_preamble = preambles[1]

    sigenergy = calculate_energy(core_preamble)

    snrs = np.arange(3, 3 * nframes, 3, dtype=np.float)
    snrs_lin = db2lin(snrs)

    for i, snr in enumerate(snrs):
        snr_lin = snrs_lin[i]

        nscale = calculate_noise_scale(
            snr_lin, sigenergy, active_ratio, core_preamble.size)
        iterations = 1000
        snr_lin0 = np.zeros(iterations, dtype=np.float)
        snr_lin1 = np.zeros(iterations, dtype=np.float)
        for i in range(iterations):
            noise = get_noise_vector(core_preamble.size, nscale)

            d = core_preamble + noise

            snr_lin0[i] = estimate_snr0(d, subcarrier_map, subcarriers)
            snr_lin1[i] = estimate_snr1(
                d, subcarrier_map, subcarriers, core_preamble)

        snr_lin0 = np.mean(snr_lin0)
        snr_lin1 = np.mean(snr_lin1)
        snr_db0 = lin2db(snr_lin0)
        snr_db1 = lin2db(snr_lin1)
        print(f"noise {snr:4.1f}dB\t{snr_db0:6.3f}dB\t{snr_db1:6.3f}dB")


def main():
    foo()


if __name__ == '__main__':
    main()
