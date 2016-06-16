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


# [2] "Bit Error Rate Performance of Generalized Frequency Division Multiplexing"
def get_data_matrix(data, K, group_by_subcarrier=False):
    # function yields data matrix according to [2]
    if group_by_subcarrier:
        # alternative grouping. Used in other papers.
        return np.reshape(data, (-1, K))
    else:
        # data grouped as described in [2]
        return np.reshape(data, (K, -1)).T


def reshape_input(x, M, K, group_by_subcarrier=True):
    '''
    1. pick every m*Kth symbol and append to output.
    2. Increase counter one time
    3. perform step 1.+2. M times
    '''
    D = get_data_matrix(x, K, group_by_subcarrier=group_by_subcarrier)
    return D.T.flatten()
