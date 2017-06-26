#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2017 <+YOU OR YOUR COMPANY+>.
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
from gnuradio import gr
from pygfdm import validation_utils

class channel_estimator_cc(gr.interp_block):
    """
    docstring for block channel_estimator_cc
    """
    def __init__(self, preamble, fft_len, timeslots, active_subcarriers):
        gr.interp_block.__init__(self,
            name="channel_estimator_cc",
            in_sig=[np.complex64],
            out_sig=[np.complex64], interp=timeslots)
        self._kernel = validation_utils.frame_estimator(preamble, fft_len, timeslots, active_subcarriers)
        self._preamble_len = 2 * fft_len
        self.set_output_multiple(fft_len * timeslots)
        self._work_count = 0


    def work(self, input_items, output_items):
        in0 = input_items[0]
        out = output_items[0]
        if len(in0) < self._preamble_len:
            return 0
        # print len(in0), len(out)
        H_estimate = self._kernel.estimate_frame(in0[0:self._preamble_len])
        out[0:len(H_estimate)] = H_estimate.astype(dtype=np.complex64)
        print(self._work_count, 'channel estimator:', len(H_estimate))
        self._work_count += 1
        return len(H_estimate)

