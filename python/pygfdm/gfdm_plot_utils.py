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
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.ticker as ticker
from matplotlib import cm


def plot_gfdm_matrix(A):
    '''
    Use this function for a 3D absolute value plot of the GFDM modulation matrix.
    Maybe helpful for debugging or for educational purposes.
    :param A: GFDM modulation matrix
    :return: None
    '''
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
    ax.set_xlabel('MxN samples')
    ax.set_ylabel('MxK samples')
    ax.set_zlabel('abs(A)')
    plt.tight_layout()
    # plt.savefig('gfdm_matrix.png')
    # plt.show()


def plot_complex(data, lc):
    plt.plot(np.real(data), lc)
    plt.plot(np.imag(data), lc + '--')


def plotFFT(x):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x_fft = np.fft.fft(x, n=x.shape[-1])
    x_t = np.fft.fftfreq(x.shape[-1])
    ax.stem(x_t, np.abs(x_fft))
    fig.show()


def plotTransmitMatrix(A):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.meshgrid(np.arange(A.shape[-1]), np.arange(A.shape[-2]))
    ax.plot_wireframe(x, y, np.abs(A), cstride=1, rstride=1)
    #ax.plot_surface(x, y, np.abs(A), rstride=1, cstride=1, cmap=cm.Greys, linewidth=0, antialiased=False)
    fig.show()


def plotSymbols(x):
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(x.real)
    ax2 = fig.add_subplot(212)
    ax2.plot(x.imag)
    fig.show()


def compareRx(*args):
    fig = plt.figure()
    ax = fig.add_subplot(211)
    for x in args:
        ax.plot(x.real)
    ax = fig.add_subplot(212)
    for x in args:
        ax.plot(x.imag)
    fig.show()


def plotScatter(ax,x):
    ax.scatter(x.real,x.imag)


def plotBER(x):
    '''
    nested dict - structure (see testsuite.py)
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for m in x:
        for k in x[m]:
            for qam in x[m][k]:
                ax.plot(x[m][k][qam][0],x[m][k][qam][1], label="M:{},K:{},QAM:{}".format(m,k,qam))
    ax.legend(bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.)
    fig.show()


def plot_sc_int(x,k_range,m_range,qam_range,j):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for k in k_range:
        for m in m_range:
            for qam in qam_range:
                ax.plot([x[k][m][qam][j_it] for j_it in xrange(j)], label="K:{}, M:{},\n QAM:{}".format(k,m,qam))
    ax.legend(bbox_to_anchor=(0.95,1), loc=2, borderaxespad=0.)
    fig.show()


def plot_waterfall(symbols, fft_len, overlap=-1, cmap=plt.get_cmap('gnuplot2'), fftshift=True):
    if overlap < 0:
        overlap = fft_len / 4
    Pxx, freqs, bins, im = plt.specgram(symbols, NFFT=fft_len, Fs=fft_len, noverlap=overlap)

    if fftshift:
        Pxx = np.fft.fftshift(Pxx, axes=0)
    plt.imshow(10 * np.log10(Pxx.T), cmap=cmap)
    plt.xlabel('FFT bins')
    plt.ylabel('time slots')
