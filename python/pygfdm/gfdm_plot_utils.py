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