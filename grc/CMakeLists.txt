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
    gfdm_framer_cc.xml
    gfdm_modulator_cc.xml
    gfdm_sync_cc.xml
    gfdm_cyclic_prefixer_cc.xml
    gfdm_preamble_generator.xml
    gfdm_remove_prefix_cc.xml
    gfdm_simple_modulator_cc.xml
    gfdm_simple_receiver_cc.xml
    gfdm_advanced_receiver_sb_cc.xml
    gfdm_advanced_receiver_tsb_cc.xml
    gfdm_resource_mapper_cc.xml
    gfdm_frame_energy_detector_cc.xml
    gfdm_simple_preamble_sync_cc.xml
    gfdm_resource_demapper_cc.xml
    gfdm_extract_burst_cc.xml
    gfdm_channel_estimator_cc.xml DESTINATION share/gnuradio/grc/blocks
)
