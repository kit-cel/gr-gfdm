#!/usr/bin/env python2.7

import matplotlib.pyplot as mp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

__all__ = ['plotFFT', 'plotTransmitMatrix']

def plotFFT(x):
    fig = mp.figure()
    ax = fig.add_subplot(111)
    x_fft = np.fft.fft(x, n=x.shape[-1])
    x_t = np.fft.fftfreq(x.shape[-1])
    ax.stem(x_t, np.abs(x_fft))
    fig.show()


def plotTransmitMatrix(A):
    fig = mp.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.meshgrid(np.arange(A.shape[-1]), np.arange(A.shape[-2]))
    ax.plot_wireframe(x, y, np.abs(A),cstride=1,rstride=1)
    #ax.plot_surface(x, y, np.abs(A), rstride=1, cstride=1, cmap=cm.Greys, linewidth=0, antialiased=False)
    fig.show()


def plotSymbols(x):
    fig = mp.figure()
    ax = fig.add_subplot(211)
    ax.plot(x.real)
    ax2 = fig.add_subplot(212)
    ax2.plot(x.imag)
    fig.show()

def compareRx(*args):
    fig = mp.figure()
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
    fig = mp.figure()
    ax = fig.add_subplot(111)
    for m in x:
        for k in x[m]:
            for qam in x[m][k]:
                ax.plot(x[m][k][qam][0],x[m][k][qam][1], label="M:{},K:{},QAM:{}".format(m,k,qam))
    ax.legend(bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.)
    fig.show()


def plot_sc_int(x,k_range,m_range,qam_range,j):
    fig = mp.figure()
    ax = fig.add_subplot(111)
    for k in k_range:
        for m in m_range:
            for qam in qam_range:
                ax.plot([x[k][m][qam][j_it] for j_it in xrange(j)], label="K:{}, M:{},\n QAM:{}".format(k,m,qam))
    ax.legend(bbox_to_anchor=(0.95,1), loc=2, borderaxespad=0.)
    fig.show()
