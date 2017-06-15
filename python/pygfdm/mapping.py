#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016 Johannes Demel.
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


def get_data_matrix(data, K, group_by_subcarrier=False):
    # [2] "Bit Error Rate Performance of Generalized Frequency Division Multiplexing"
    # function yields data matrix according to [2]
    if group_by_subcarrier:
        # alternative grouping. Used in other papers.
        return np.reshape(data, (-1, K))
    else:
        # data grouped as described in [2]
        return np.reshape(data, (K, -1)).T


def get_data_stream(data, group_by_subcarrier=False):
    if group_by_subcarrier:
        return np.reshape(data,(-1,1)).T[0]
    else:
        return np.reshape(data.T, (1,-1))[0]


def reshape_input(x, M, K, group_by_subcarrier=True):
    '''
    1. pick every m*Kth symbol and append to output.
    2. Increase counter one time
    3. perform step 1.+2. M times
    '''
    D = get_data_matrix(x, K, group_by_subcarrier=group_by_subcarrier)
    return D.T.flatten()


def map_to_waveform_resources(syms, active_subcarriers, fft_len, subcarrier_map, per_timeslot=True):
    frame = map_to_waveform_resource_grid(syms, active_subcarriers, fft_len, subcarrier_map, per_timeslot)
    return frame.flatten()


def demap_from_waveform_resource_grid(syms, subcarriers, subcarrier_map):
    tm = np.reshape(syms, (subcarriers, -1))
    td = tm[subcarrier_map, :]
    return td.T.flatten()


def map_to_waveform_resource_grid(syms, active_subcarriers, fft_len, subcarrier_map, per_timeslot=True):
    n_data_syms = len(syms)
    ts = int(np.ceil(1. * n_data_syms / active_subcarriers))
    n_frame_syms = active_subcarriers * ts
    s = np.concatenate((syms, np.zeros(n_frame_syms - n_data_syms)))
    if per_timeslot:
        s = np.reshape(s, [-1, active_subcarriers]).T
    else:
        s = np.reshape(s, [active_subcarriers, -1])
    frame = np.zeros((fft_len, ts), dtype=np.complex64)
    frame[subcarrier_map, :] = s
    return frame


def get_subcarrier_map(subcarriers, active_subcarriers, dc_free=False):
    if dc_free:
        return np.concatenate((np.arange(1, active_subcarriers // 2 + 1), np.arange(subcarriers - active_subcarriers // 2, subcarriers)))
    return np.concatenate((np.arange(0, active_subcarriers // 2), np.arange(subcarriers - active_subcarriers // 2, subcarriers)))


def resource_mapping_test():
    subcarriers = 8
    timeslots = 4
    smap = get_subcarrier_map(subcarriers, subcarriers)
    group_tests = [False, True]
    for g in group_tests:
        d = np.arange(subcarriers * timeslots, dtype=np.complex64)
        D = get_data_matrix(d, subcarriers, g)
        m = map_to_waveform_resource_grid(d, subcarriers, subcarriers, smap, g)
        if not np.all(D.T == m):
            raise ValueError('mapping matrices do not match!')

        f = map_to_waveform_resources(d, subcarriers, subcarriers, smap, g)
        F = reshape_input(d, timeslots, subcarriers, g)
        if not np.all(f == F):
            raise ValueError('stacked mapping vectors do not match!')


def main():
    np.set_printoptions(precision=2)
    resource_mapping_test()


if __name__ == '__main__':
    main()
