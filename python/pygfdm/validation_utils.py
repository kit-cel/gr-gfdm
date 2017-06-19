#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Johannes Demel.
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
import utils
import mapping
import preamble
import gfdm_modulation
import cyclic_prefix
import filters


def generate_reference_frame(timeslots, subcarriers, active_subcarriers, cp_len, cs_len, alpha=.2):
    p_seed = utils.generate_seed('awesome preamble')
    f_seed = utils.generate_seed('awesome frame')
    subcarrier_map = mapping.get_subcarrier_map(subcarriers, active_subcarriers, dc_free=True)
    overlap = 2
    frame_preamble, x_preamble = preamble.mapped_preamble(p_seed, 'rrc', alpha, active_subcarriers, subcarriers, subcarrier_map, overlap, cp_len, cs_len)
    d = utils.get_random_qpsk(timeslots * active_subcarriers, f_seed)
    d_frame = mod_frame = gfdm_modulation.modulate_mapped_gfdm_block(d, timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)
    symbol = cyclic_prefix.add_cyclic_starfix(d_frame, cp_len, cs_len)
    window_ramp = cyclic_prefix.get_raised_cosine_ramp(cs_len, cyclic_prefix.get_window_len(cp_len, timeslots, subcarriers, cs_len))
    d_frame = cyclic_prefix.pinch_block(symbol, window_ramp)

    H = filters.get_frequency_domain_filter('rrc', alpha, timeslots, subcarriers, overlap)
    return np.concatenate((frame_preamble, d_frame)), mod_frame, x_preamble, d, H


def generate_sc_qpsk_frame(timeslots, subcarriers, active_subcarriers, cp_len, cs_len, alpha=.2):
    p_seed = utils.generate_seed('awesome preamble')
    f_seed = utils.generate_seed('awesome frame')
    subcarrier_map = mapping.get_subcarrier_map(subcarriers, active_subcarriers, dc_free=True)
    overlap = 2
    frame_preamble, x_preamble = preamble.mapped_preamble(p_seed, 'rrc', alpha, active_subcarriers, subcarriers, subcarrier_map, overlap, cp_len, cs_len)
    d = utils.get_random_qpsk(timeslots * subcarriers, f_seed)
    # d_frame = mod_frame = gfdm_modulation.modulate_mapped_gfdm_block(d, timeslots, subcarriers, active_subcarriers, overlap, alpha, dc_free=True)
    # symbol = cyclic_prefix.add_cyclic_starfix(d_frame, cp_len, cs_len)

    window_ramp = cyclic_prefix.get_raised_cosine_ramp(cs_len, cyclic_prefix.get_window_len(cp_len, timeslots, subcarriers, cs_len))
    d_frame = cyclic_prefix.pinch_block(d, window_ramp)

    H = filters.get_frequency_domain_filter('rrc', alpha, timeslots, subcarriers, overlap)
    return np.concatenate((frame_preamble, d_frame)), d, x_preamble, d, H