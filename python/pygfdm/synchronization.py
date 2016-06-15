#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Andrej Rode.
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
import commpy as cp
import matplotlib.pyplot as plt
import scipy.signal as signal
from modulation import gfdm_tx, gfdm_tx_fft2
from filters import get_frequency_domain_filter
from gfdm_modulation import gfdm_modulate_block
from cyclic_prefix import add_cyclic_prefix, pinch_block, get_raised_cosine_ramp, get_window_len
from mapping import get_data_matrix, get_random_qpsk


def sync_symbol(filtertype, alpha, K, n_mod, N):
    '''
        Generate Schmidls Training Symbols to achieve Receiver Synchronisation
        K: should be an odd number
        Process:
            * First Symbol: Transmit PN-Sequence on all even frequencies while zeros on the odd
            frequencies. Constant signal energy -> Multiply every Symbol with sqrt(2)

            * Second Symbol: Transmit PN-Sequence on all odd frequencies and another PN-Sequence
            on the even frequencies
    '''

    pn_order = 14
    pn_seed = '00101101010010'
    pn_mask = '01001110100111'
    if int(np.floor(K / 2.0)) % 2:
        n_even_freq = int(np.floor(K / 2.0))
    else:
        n_even_freq = int(np.ceil(K / 2.0))
    seq_length = n_even_freq * n_mod
    sym_sequence = np.zeros(K, dtype='complex')
    pn_sequence = cp.pnsequence(pn_order, pn_seed, pn_mask, seq_length)
    qam_mod = cp.modulation.QAMModem(2 ** n_mod)
    qam_sequence = qam_mod.modulate(pn_sequence)
    for i in xrange(len(sym_sequence)):
        if not i % 2:
            sym_sequence[i] = qam_sequence[i / 2]
    ifft_sequence = np.fft.ifftshift(sym_sequence)
    output = gfdm_tx(ifft_sequence, filtertype, alpha, 1, K, N)
    # output = np.fft.ifft(np.sqrt(2)*ifft_sequence)

    return output


def sync_symbol2(filtertype, alpha, K, L, n_mod):
    pn_order = 14
    pn_seed = '01001000111011'
    pn_mask = '01001101001110'
    seq_length = K * n_mod
    pn_sequence = cp.pnsequence(pn_order, pn_seed, pn_mask, seq_length)
    qam_mod = cp.modulation.QAMModem(2 ** n_mod)
    qam_sequence = qam_mod.modulate(pn_sequence)
    output = gfdm_tx_fft2(np.tile(qam_sequence, 2), filtertype, alpha, 2, K, 2, 1)
    return output


def sync_product(x, L):
    '''
    Auto-Korrelation der ersten L Samples mit den naechsten L Samples
    '''
    return np.sum([x[i].conj() * x[i + L] for i in xrange(L)])


def sync_iter(x, L, cp):
    '''
    Schrittweise Iteration ueber alle Samples (zunaechst die ersten 2*L Samples)
    Danach ueber die restlichen len(x)-2*L Samples
    '''
    P_d = np.array([])
    P_d = np.append(P_d, sync_product(x, L))
    for i in xrange(len(x) - 2 * L):
        P_d = np.append(
            P_d, P_d[i] + (x[L + i].conj() * x[2 * L + i]) -
                 (x[i].conj() * x[L + i]))
    P_d_out = P_d
    P_d = np.append(np.zeros(cp, dtype='complex'), P_d)
    # Integrate cp-samples to eliminate cp-plateau
    P_di = np.array([])
    for i in xrange(cp, len(x) - 2 * L):
        P_di = np.append(
            P_di, (1.0 / (cp + 1) * np.sum(np.abs(P_d[i - cp:i]) ** 2)))
    return (P_di, P_d_out)


def sync_energy(x, L):
    '''
    Berechnung der Energie der zweiten Haelfte der Sync-Samples -> normieren
    '''
    R_d = np.array([])
    R_d = np.append(R_d, np.sum([np.abs(x[i + L]) ** 2 for i in xrange(L)]))
    for i in xrange(len(x) - 2 * L):
        R_d = np.append(R_d, R_d[-1] + np.abs(x[2 * L + i]) ** 2 - np.abs(x[L + i]) ** 2)
    return R_d


def sync_perform(x, L, cp):
    (P_di, P_d) = sync_iter(x, L, cp)
    # R_d = sync_energy(x, L)
    # M_d = (np.abs(P_d)**2)/(R_d**2)
    return (P_di, P_d)


def sync_CFO(P_d, P_di):
    '''
    Gewinn von d (Beginn der Pilotsequenz) und d_f Frequenzoffset genormt auf 1/T.
    Kann bei nur einem Pilotsymbol nur +/- 1/T betragen.
    '''
    d = np.argmax(P_di)
    dphi = np.angle(P_d[d])
    d_f = dphi / (np.pi)
    print("P_d:{},df:{})".format(P_d[d], d_f))

    return (d, d_f)


def get_sync_symbol(pn_symbols, H, K, L, cp_len, ramp_len):
    M = 2  # fixed for preamble
    compat_mode = False
    pn_symbols = np.concatenate((pn_symbols, pn_symbols))
    D = get_data_matrix(pn_symbols, K, group_by_subcarrier=compat_mode)
    symbol = gfdm_modulate_block(D, H, M, K, L, compat_mode=compat_mode)
    print np.shape(symbol)
    symbol = add_cyclic_prefix(symbol, cp_len)
    print 'prefixed', np.shape(symbol)
    window_ramp = get_raised_cosine_ramp(ramp_len, get_window_len(cp_len, M, K))
    print 'ramp', np.shape(window_ramp)
    symbol = pinch_block(symbol, window_ramp)
    print np.shape(symbol)
    return symbol


def main():
    print 'Hello World'
    alpha = .5
    M = 2
    K = 32
    L = 2
    cp_len = 4
    ramp_len = 4
    pn_symbols = get_random_qpsk(K)
    H = get_frequency_domain_filter('rrc', alpha, M, K, L)
    symbol = get_sync_symbol(pn_symbols, H, K, L, cp_len, ramp_len)
    symbol *= 1. / np.sqrt(np.sum(symbol ** 2))
    print np.shape(symbol)
    c = signal.correlate(symbol, symbol)
    # c *= 1. / len(symbol)
    print np.shape(np.abs(c))
    nc = c[np.argmax(np.abs(c))]
    print nc
    plt.plot(np.abs(c))
    plt.show()


if __name__ == '__main__':
    main()
