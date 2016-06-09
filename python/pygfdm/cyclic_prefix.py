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


def window_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros=True):
    # returns symmetric ram function for block pinching. edges excluded.
    block_len = n_timeslots * n_subcarriers
    window_len = block_len + cp_len
    alpha = 1. * ramp_len / block_len
    w = alpha + (1. - alpha) / 2.
    if cut_zeros:
        f = np.linspace(-w, w, window_len + 2, dtype=np.float)
        f = f[1:-1]
    else:
        f = np.linspace(-w, w, window_len, dtype=np.float)
    r = truncated_lin(f, alpha)
    # print ramp_len
    # print r
    return f, r


def get_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros=True):
    f, r = window_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros)
    x = .5 * (1. + np.cos(np.pi * r))
    return f, x


def get_fourth_order_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros=True):
    f, r = window_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros)
    r = fourth_order_polynomial(r)
    x = .5 * (1. + np.cos(np.pi * r))
    return f, x


def get_root_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros=True):
    f, r = get_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers, cut_zeros)
    r = np.sqrt(r)
    return f, r


def fourth_order_polynomial(x):
    return (x ** 4) * (35 - 84 * x + 70 * (x ** 2) - 20 * (x ** 3))


def example_raised_cosine(alpha, nsamples):
    w = alpha + (1. - alpha) / 2.
    step_size = 2. * w / (nsamples)
    f = np.arange(-w, w + step_size, step_size)
    np.linspace(-w, w, nsamples, dtype=np.float)
    x = truncated_lin(f, alpha)
    x = .5 * (1. + np.cos(np.pi * x))
    return f, x


def truncated_lin(x, alpha):
    # alpha: roll-off
    # partially described in [0]. Made consistent with [1], [2], [3]
    f = (1. - alpha) / (2. * alpha)
    x = np.abs(x) / alpha - f
    x = np.maximum(np.zeros(len(x)), x)
    x = np.minimum(np.ones(len(x)), x)
    return x


def pinch_block(data, window_taps):
    return data * window_taps


def add_cyclic_starfix(data, cp_len, cs_len):
    # b = np.concatenate(data[-cp_len:], data)
    return np.concatenate((data[-cp_len:], data, data[0:cs_len]))


def add_cyclic_prefix(data, cp_len):
    return add_cyclic_starfix(data, cp_len, 0)


def plot_window_ramps():
    n_subcarriers = 8
    n_timeslots = 16
    cp_len = n_subcarriers # / 2
    ramp_len = cp_len #/ 2

    f, d = window_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers)
    print len(f), len(d)
    plt.plot(f, d, label='window ramp')
    df = fourth_order_polynomial(d)
    plt.plot(f, df, label='fourth order ramp')
    f, d = get_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers)
    plt.plot(f, d, label='RC ramp')
    plt.plot(*get_fourth_order_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers), label='4th RC ramp')
    plt.plot(*get_root_raised_cosine_ramp(ramp_len, cp_len, n_timeslots, n_subcarriers), label='RRC ramp')

    plt.plot(np.arange(-1, 2), np.ones(3) * .5)
    plt.grid()
    plt.legend()
    plt.show()



def main():
    plot_window_ramps()
    # data = np.arange(16)
    #
    # f, d = window_ramp(8, 8, 8, 4)
    # print len(f), len(d)
    # plt.plot(f, d)
    # plt.plot(f, fourth_order_polynomial(d))
    #
    # n_samps = 100
    # for alpha in np.arange(.05, 1., .2):
    #     w = alpha + (1. - alpha) / 2.
    #     plt.axvline(w)
    #
    #     f, d = example_raised_cosine(alpha, n_samps)
    #
    #     plt.plot(f, d, label=alpha)
    #
    # plt.plot(np.arange(-1, 2), np.ones(3) * .5)
    # plt.grid()
    # plt.legend()
    # plt.show()


if __name__ == '__main__':
    main()
