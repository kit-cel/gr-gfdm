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

'''
This file holds functions necessary for cyclic prefix/suffix insertion. Also, it provides tools for block pinching.

[0] Generalized Frequency Division Multiplexing for 5th Generation Cellular Networks
[1] https://en.wikipedia.org/wiki/Raised-cosine_filter
[2] https://en.wikipedia.org/wiki/Root-raised-cosine_filter
[3] A Survey on Multicarrier Communications: Prototype Filters, Lattice Structures and Implementation Aspects
'''

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def get_window_len(cp_len, n_timeslots, n_subcarriers):
    return n_timeslots * n_subcarriers + cp_len


def fourth_order_polynomial(x):
    return (x ** 4) * (35 - 84 * x + 70 * (x ** 2) - 20 * (x ** 3))


def window_ramp(ramp_len, window_len):
    r = np.arange(0, 1, 1. / ramp_len)
    return np.concatenate((1. - r, np.zeros(window_len - 2 * ramp_len), r))


def calculate_raised_cosine(r):
    return .5 * (1. + np.cos(np.pi * r))


def get_raised_cosine_ramp(ramp_len, window_len):
    r = window_ramp(ramp_len, window_len)
    x = calculate_raised_cosine(r)
    return x


def get_root_raised_cosine_ramp(ramp_len, window_len):
    r = get_raised_cosine_ramp(ramp_len, window_len)
    return np.sqrt(r)


def get_fourth_order_raised_cosine_ramp(ramp_len, window_len):
    r = window_ramp(ramp_len, window_len)
    r = fourth_order_polynomial(r)
    x = calculate_raised_cosine(r)
    return x


def pinch_block(data, window_taps):
    return data * window_taps


def add_cyclic_starfix(data, cp_len, cs_len):
    # b = np.concatenate(data[-cp_len:], data)
    return np.concatenate((data[-cp_len:], data, data[0:cs_len]))


def add_cyclic_prefix(data, cp_len):
    return add_cyclic_starfix(data, cp_len, 0)


def pinch_cp_add_block(data, timeslots, subcarriers, cp_len, ramp_len):
    window_len = get_window_len(cp_len, timeslots, subcarriers)
    window = get_root_raised_cosine_ramp(ramp_len, window_len)
    d = add_cyclic_prefix(data, cp_len)
    return pinch_block(d, window)


def plot_window_ramps():
    n_subcarriers = 16
    n_timeslots = 15
    cp_len = n_subcarriers # / 2
    ramp_len = cp_len #/ 2
    window_len = get_window_len(cp_len, n_timeslots, n_subcarriers)

    d = window_ramp(ramp_len, window_len)
    plt.plot(d, label='window ramp')
    df = fourth_order_polynomial(d)
    plt.plot(df, label='fourth order ramp')
    d = get_raised_cosine_ramp(ramp_len, window_len)
    plt.plot(d, label='RC ramp')
    plt.plot(get_fourth_order_raised_cosine_ramp(ramp_len, window_len), label='4th RC ramp')
    plt.plot(get_root_raised_cosine_ramp(ramp_len, window_len), label='RRC ramp')

    plt.plot(np.arange(-1, 2), np.ones(3) * .5)
    plt.grid()
    plt.legend()
    plt.show()


def main():
    plot_window_ramps()


if __name__ == '__main__':
    main()
