# Copyright 2011 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
install(FILES
    gfdm_modulator_cc.block.yml
    gfdm_cyclic_prefixer_cc.block.yml
    gfdm_preamble_generator.block.yml
    gfdm_remove_prefix_cc.block.yml
    gfdm_simple_modulator_cc.block.yml
    gfdm_simple_receiver_cc.block.yml
    gfdm_advanced_receiver_sb_cc.block.yml
    gfdm_resource_mapper_cc.block.yml
    gfdm_resource_demapper_cc.block.yml
    gfdm_extract_burst_cc.block.yml
    gfdm_channel_estimator_cc.block.yml
    gfdm_transmitter_cc.block.yml
    gfdm_short_burst_shaper.block.yml DESTINATION share/gnuradio/grc/blocks
)
