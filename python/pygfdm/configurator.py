#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from collections import namedtuple
import numpy as np
from .mapping import get_subcarrier_map
from .filters import get_frequency_domain_filter
from .preamble import mapped_preamble
from .gfdm_modulation import modulate_mapped_gfdm_block
from .cyclic_prefix import add_cyclic_starfix, get_raised_cosine_ramp, pinch_block, get_window_len


def round_up_power_of_2(value):
    return int(2 ** np.ceil(np.log2(float(value))))


def get_padding_configuration(frame_len):
    padded_frame_len = round_up_power_of_2(frame_len)
    if padded_frame_len - frame_len < 500:
        padded_frame_len *= 2
    padded_len = padded_frame_len - frame_len
    pre_padding_len = 256
    post_padding_len = 128
    while pre_padding_len + post_padding_len < padded_len:
        pre_padding_len += 128
        post_padding_len += 128
    post_padding_len -= pre_padding_len + post_padding_len - padded_len
    return pre_padding_len, post_padding_len


preamble_seed = int(3660365253)


def get_gfdm_configuration(timeslots=9, subcarriers=64, active_subcarriers=52, overlap=2,
                           cp_len=16, cs_len=8, filtertype='rrc', filteralpha=0.2,
                           cyclic_shifts=[0, ]):
    ramp_len = cs_len
    dconf = {
        'timeslots': timeslots,
        'subcarriers': subcarriers,
        'active_subcarriers': active_subcarriers,
        'overlap': overlap,
        'cp_len': cp_len,
        'cs_len': cs_len,
        'ramp_len': ramp_len,
        'cyclic_shifts': cyclic_shifts,
        'seed': preamble_seed,
    }
    dconf['block_len'] = block_len = timeslots * subcarriers
    dconf['window_len'] = window_len = block_len + cp_len + cs_len
    # for k, v in dconf.items():
    #    print('{:16s}:\t{}'.format(k, v))

    dconf['subcarrier_map'] = subcarrier_map = get_subcarrier_map(
        subcarriers, active_subcarriers, dc_free=True)
    preambles = [mapped_preamble(preamble_seed, filtertype, filteralpha, active_subcarriers,
                                 subcarriers, subcarrier_map, overlap, cp_len, ramp_len, use_zadoff_chu=True,
                                 cyclic_shift=cs) for cs in cyclic_shifts]
    dconf['full_preambles'] = full_preambles = [p[0] for p in preambles]
    dconf['core_preamble'] = core_preamble = preambles[0][1]
    dconf['preamble_len'] = preamble_len = full_preambles[0].size
    dconf['core_preamble_len'] = core_preamble_len = core_preamble.size
    dconf['frame_len'] = frame_len = window_len + preamble_len
    pre_padding_len, post_padding_len = get_padding_configuration(frame_len)
    dconf['pre_padding_len'] = pre_padding_len
    dconf['post_padding_len'] = post_padding_len
    dconf['padded_frame_len'] = pre_padding_len + frame_len + post_padding_len

    dconf['window_taps'] = window_taps = get_raised_cosine_ramp(
        ramp_len, get_window_len(cp_len, timeslots, subcarriers, cs_len))

    dconf['tx_filter_taps'] = tx_filter_taps = get_frequency_domain_filter(
        filtertype, filteralpha, timeslots, subcarriers, overlap)
    dconf['rx_filter_taps'] = rx_f_taps = np.conjugate(tx_filter_taps)

    co = namedtuple("configuration", dconf.keys())(*dconf.values())
    return co
