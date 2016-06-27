# Copyright 2011,2012 Free Software Foundation, Inc.
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

########################################################################
# Install public header files
########################################################################
install(FILES
    api.h
    framer_cc.h
    modulator_cc.h
    gfdm_utils.h
    receiver_kernel_cc.h
    #receiver_cc.h
    advanced_receiver_cc.h
    sync_cc.h
    cyclic_prefixer_cc.h
    preamble_generator.h
    remove_prefix_cc.h
    simple_modulator_cc.h
    modulator_kernel_cc.h
    add_cyclic_prefix_cc.h
    improved_sync_algorithm_kernel_cc.h DESTINATION include/gfdm
)
