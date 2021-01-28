#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2021 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

'''
This file exists for debug/testing/QA!

If you want something flexible and fast:
- gr-symbolmapping

'''
import numpy as np


# dict key is constellation order
constellations = {1: np.array([1.+0.j, -1.+0.j]),
                  2: np.array([1+1j, 1-1j, -1+1j, -1-1j]) / np.sqrt(2)}


def generate_constellation(constellation_order):
    return constellations[constellation_order]


def pack_bits(bits, n_packed):
    b = np.reshape(bits, (-1, n_packed))
    bytes = np.packbits(b, axis=1)
    return np.right_shift(bytes, 8 - n_packed)
    # vals = 2 ** np.arange(n_packed - 1, -1, -1)
    # return np.sum(b * vals, axis=1)


def bits2symbols(bits, constellation):
    constellation_order = int(np.log2(len(constellation)))
    points = pack_bits(bits, constellation_order)
    symbols = np.array([constellation[i] for i in points])
    return symbols


def symbols2bits(symbols, constellation):
    constellation_order = int(np.log2(len(constellation)))
    points = np.array([np.argmin(np.abs(s - constellation) ** 2) for s in symbols], dtype=np.uint8)
    bits = np.array([np.unpackbits(p)[-constellation_order:] for p in points]).flatten()
    return bits


def main():
    print(constellations[1])
    print(constellations[2])

    bits = np.random.randint(0, 2, 24 * 4)
    symbols = bits2symbols(bits, constellations[2])
    rx_bits = symbols2bits(symbols, constellations[2])
    print(bits)
    print(rx_bits)

    assert(np.all(bits == rx_bits))


if __name__ == '__main__':
    main()
