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

import pygfdm
from pygfdm.modulation import gfdm_filter_taps, transmitMatrix_core
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.ticker as ticker
from matplotlib import cm

import commpy as cp


def gfdm_modulation_matrix(p, M, K, oversampling_factor=1):
    # generate GFDM modulation matrix from prototype filter taps p, M time slots and K subcarriers
    # p must hold N = K * M taps!
    N = M * K

    p = np.roll(p, (N * oversampling_factor) // 2)
    A = np.zeros((N * oversampling_factor, N), dtype=np.complex)

    n = np.arange(N * oversampling_factor, dtype=np.complex)
    for m in range(M):
        g = np.roll(p, m * K * oversampling_factor)
        for k in range(K):
            g = g * np.exp(1j * 2 * np.pi * (float(k) / (K * oversampling_factor)) * n)
            A[:, m * K + k] = g
    return A


def plot_gfdm_matrix(A):
    (Nx, Ny) = np.shape(A)

    ax = plt.figure().gca(projection='3d')
    X = np.arange(Nx)
    Y = np.arange(Ny)
    X, Y = np.meshgrid(X, Y)
    Z = A

    ax.plot_trisurf(X.flatten(), Y.flatten(), np.abs(Z).flatten(),  # rstride=1, cstride=1,
                    cmap=cm.viridis,
                    linewidth=0, antialiased=True)

    ax.set_zlim(0.0, 1.0)
    ax.zaxis.set_major_locator(ticker.LinearLocator(10))
    ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.02f'))
    ax.set_zlabel('abs(A)')
    plt.tight_layout()
    # plt.savefig('gfdm_matrix.png')
    # plt.show()


def filter_test():
    M = 8
    K = 4
    N = 1
    alpha = .1
    t, rrc_taps = cp.rcosfilter(M * K * N, alpha, N * K, 2)
    print t
    print rrc_taps
    # if filtertype == "rrc":
    #     time_h, h = cp.rrcosfilter(M * K * N, alpha, N * K, 1)
    # elif filtertype == "rc":
    #     time_h, h = cp.rcosfilter(M * K * N, alpha, N * K, 1)

    plt.plot(t, rrc_taps)
    plt.show()


def main():
    print "hello"
    M = 6
    K = 3
    alpha = 1.0
    oversampling_factor = 2
    t_extract = 2

    taps = p = gfdm_filter_taps('rrc', alpha, M, K, oversampling_factor)
    A0 = gfdm_modulation_matrix(p, M, K, oversampling_factor)
    print 'GFDM shape: ', np.shape(A0)

    indices = np.arange(M * K)
    print indices
    indices = np.reshape(indices, (-1, K)).T.flatten()
    print indices
    A0 = A0[:, indices]
    plot_gfdm_matrix(A0)

    A1 = transmitMatrix_core(taps, M, K, oversampling_factor)
    scaling_factor = abs(A0[0, 0]) / abs(A1[0, 0])
    A1 *= scaling_factor
    print 'GR shape: ', np.shape(A1)
    print 'scaling factor:', scaling_factor

    plot_gfdm_matrix(A1)
    # plot_gfdm_matrix(abs(A0) - abs(A1))
    plot_gfdm_matrix(A0 - A1)

    t0 = A0[:, t_extract]
    t1 = A1[:, t_extract]


    # plt.plot(taps)
    # fig = plt.figure()
    # t0 *= (abs(t1[0]) / abs(t0[0]))
    # plt.plot(np.real(t0))
    # plt.plot(np.real(t1))
    plt.show()


if __name__ == '__main__':
    main()
