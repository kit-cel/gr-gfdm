#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from fractions import gcd

import synchronization as sync
import preamble
import converter

import matplotlib.pyplot as plt


def main():
    np.set_printoptions(precision=2, linewidth=150)
    seed = abs(hash('awesome')) % (2 ** 32)
    subcarrier_map = 
    frame_preamble, x_preamble = preamble.mapped_preamble(seed, 'rrc', .5, 52, 64, subcarrier_map, 2, 64, 16)
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
