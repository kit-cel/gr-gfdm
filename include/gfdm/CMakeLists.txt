# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file was generated by gr_modtool, a tool from the GNU Radio framework
# This file is a part of gr-gfdm
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Install public header files
########################################################################
install(FILES
    api.h
    modulator_cc.h
    gfdm_utils.h
    receiver_kernel_cc.h
    cyclic_prefixer_cc.h
    remove_prefix_cc.h
    simple_modulator_cc.h
    modulator_kernel_cc.h
    add_cyclic_prefix_cc.h
    simple_receiver_cc.h
    advanced_receiver_sb_cc.h
    resource_mapper_kernel_cc.h
    resource_mapper_cc.h
    resource_demapper_cc.h
    advanced_receiver_kernel_cc.h
    extract_burst_cc.h
    preamble_channel_estimator_cc.h
    gfdm_kernel_utils.h
    channel_estimator_cc.h
    transmitter_kernel.h
    transmitter_cc.h
    short_burst_shaper.h DESTINATION include/gfdm
)
