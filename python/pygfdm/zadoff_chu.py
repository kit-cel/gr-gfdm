#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
import numpy as np


def generate_zadoff_chu_sequence(seq_length, uvalue, shift_value=0):
    '''
    Read up on Wikipedia https://en.wikipedia.org/wiki/Zadoff%E2%80%93Chu_sequence
    '''
    if not np.gcd(seq_length, uvalue) == 1:
        raise RuntimeError(f'GCD(N_ZC={seq_length}, u={uvalue}) != 1 !')

    if not 0 < uvalue < seq_length:
        raise RuntimeError(f'Does not satisfy: 0 < u={uvalue} < N_ZC!')
    c_f = seq_length % 2
    n = np.arange(seq_length)
    vec = n * (n + c_f + 2 * shift_value)
    seq = np.exp(-1.j * np.pi * vec / seq_length)
    return seq


if __name__ == "__main__":
    seq = generate_zadoff_chu_sequence(64, 19)
    print(seq)
