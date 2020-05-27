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
    import matplotlib.pyplot as plt
    from utils import get_random_qpsk
    nsc = 64
    seq = generate_zadoff_chu_sequence(nsc, 19)
    print(seq)
    nasc = 52
    sc_mask = np.concatenate((np.arange(nsc - nasc // 2, nsc), np.arange(1, 1 + nasc // 2)))
    sc_mask = np.concatenate((np.arange(1, 1 + nasc // 2), np.arange(nsc - nasc // 2, nsc)))
    print(sc_mask)
    print(sc_mask.shape)

    seq0 = np.zeros_like(seq)
    seq0[sc_mask] = seq[sc_mask]

    td_seq = np.fft.ifft(seq0)

    td_seq = np.tile(td_seq, 2)

    papr = np.max(np.abs(td_seq) ** 2) / np.mean(np.abs(td_seq) ** 2)
    print(f'PUNC: PAPRlin={papr:.4f}, {10. * np.log10(papr):.4f}dB')
    # plt.plot(np.abs(td_seq))
    plt.plot(np.abs(np.correlate(td_seq, td_seq, 'full')))

    seq1 = generate_zadoff_chu_sequence(nasc, 19)
    # seq1 = np.delete(seq1, nasc // 2)
    td_seq1 = np.zeros_like(seq)
    td_seq1[sc_mask] = seq1

    # td_seq1[1:52 // 2 + 1] = seq1[0: 52 // 2]
    # td_seq1[-52 // 2:] = seq1[52 // 2:]
    td_seq1 = np.fft.ifft(td_seq1)
    # print(td_seq1)
    td_seq1 = np.tile(td_seq1, 2)

    papr = np.max(np.abs(td_seq1) ** 2) / np.mean(np.abs(td_seq1) ** 2)
    print(f'USED: PAPRlin={papr:.4f}, {10. * np.log10(papr):.4f}dB')
    # plt.plot(np.abs(td_seq1))
    # plt.show()
    plt.plot(np.abs(np.correlate(td_seq1, td_seq1, 'full')))

    seq2 = generate_zadoff_chu_sequence(nsc, 19)
    fd_seq2 = np.fft.fft(seq2)
    s = np.zeros_like(seq2)
    s[sc_mask] = fd_seq2[sc_mask]
    td_seq2 = np.fft.ifft(s) / np.sqrt(nsc)
    td_seq2 = np.tile(td_seq2, 2)
    papr = np.max(np.abs(td_seq2) ** 2) / np.mean(np.abs(td_seq2) ** 2)
    print(f'TD:   PAPRlin={papr:.4f}, {10. * np.log10(papr):.4f}dB')
    plt.plot(np.abs(np.correlate(td_seq2, td_seq2, 'full')))

    d = get_random_qpsk(nasc)
    s = np.zeros_like(seq2)
    s[sc_mask] = d

    td_seq3 = np.fft.ifft(s)
    td_seq3 = np.tile(td_seq3, 2)
    papr = np.max(np.abs(td_seq3) ** 2) / np.mean(np.abs(td_seq3) ** 2)
    print(f'QPSK: PAPRlin={papr:.4f}, {10. * np.log10(papr):.4f}dB')
    plt.plot(np.abs(np.correlate(td_seq3, td_seq3, 'full')))

    plt.show()
