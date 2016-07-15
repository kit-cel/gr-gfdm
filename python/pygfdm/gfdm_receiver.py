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
from modulation import gfdm_modulation_matrix
from mapping import get_data_stream
from filters import get_frequency_domain_filter

def gfdm_transform_input_to_fd(R):
    return np.fft.fft(R.astype(np.complex))

def gfdm_extract_subcarriers(R,K,M,L):
    D = np.empty((K,M*L),np.complex)
    for k in xrange(K):
        for l in xrange(L):
            start_pos = ((k+l+K-1)%K)*M
            copy_pos = ((l+L/2) % L)*M
            D[k][copy_pos:copy_pos+M] = R[start_pos:start_pos+M]
    return D

def gfdm_filter_subcarriers(R,H,K,M,L):
    D = R.flatten()
    F = D * np.tile(H,K)
    return np.reshape(F,(-1,L*M)).T

def gfdm_superposition_subcarriers(R,K,M,L):
    S = np.zeros((K,M),np.complex)
    D = R.T
    for k in xrange(K):
        S[k] = np.fft.ifft(np.sum(np.reshape(D[k],(L,-1)),axis=0))
    return S.T

def gfdm_demodulate_block(R,H,K,M,L):
    D_0 = gfdm_transform_input_to_fd(R)
    D_1 = gfdm_extract_subcarriers(D_0,K,M,L)
    D_2 = gfdm_filter_subcarriers(D_1,H,K,M,L)
    D_3 = gfdm_superposition_subcarriers(D_2,K,M,L)
    return get_data_stream(D_3)

def gfdm_demodulate_fft(data,alpha,M,K,overlap):
    H = get_frequency_domain_filter('rrc', alpha, M, K, overlap)
    return gfdm_demodulate_block(data,H.conj(),K,M,overlap)
