#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from fractions import gcd

import synchronization as sync
import preamble


try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


def load_frame(filename):
    f = np.load(filename)
    f = f.item()
    frame = f['frame']
    return frame


def get_iq_max(signal):
    sr_max = np.max(signal.real)
    si_max = np.max(signal.imag)
    s_max = max(sr_max, si_max)
    return s_max


def convert_to_sc16(signal):
    # assume 12bit ADC craze
    ii16 = np.iinfo(np.int16)
    max_i16 = ii16.max
    max_i16 = min(max_i16, 2 ** 11)
    sig_max = max_i16 * .9
    s_max = get_iq_max(signal)
    signal *= sig_max / s_max
    i16 = signal.real.astype(np.int16)
    q16 = signal.imag.astype(np.int16)
    sc16 = np.dstack((i16, q16)).flatten()
    return sc16


def convert_from_sc16(signal):
    i = signal[0::2].astype(np.float)
    q = signal[1::2].astype(np.float)
    return i + 1j * q


def convert_to_cf64(signal):
    return signal.astype(np.complex64)


def load_gr_iq_file(filename):
    return np.fromfile(filename, dtype=np.complex64)


def main():
    np.set_printoptions(precision=2, linewidth=150)
    filename = '/home/demel/iq_samples/gfdm_samples_part0.dat'
    subcarriers = 64
    frame = np.fromfile(filename, dtype=np.complex64)[25000:225000]
    print('M * K', 9 * 64, ' + 2 * K', 2 * 64, ' + 2CP + 2CS', 2 * 64 + 2 * 16)
    print(len(frame))
    ac = sync.auto_correlate_signal(frame, subcarriers)

    plt.plot(np.abs(ac))
    plt.plot(np.abs(frame) * 10)

    # plt.plot(np.abs(frame))
    #
    plt.show()

if __name__ == '__main__':
    main()
