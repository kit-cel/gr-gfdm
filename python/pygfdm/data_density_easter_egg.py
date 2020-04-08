#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
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


def data_per_volume():
    '''
    This is a little funtime activity easter egg trying to estimate data content of the human genome.
    '''
    micro_sd_dimensions = [11e-3, 15e-3, 1e-3]  # according to spec
    micro_sd_volume = micro_sd_dimensions[0] * micro_sd_dimensions[1] * micro_sd_dimensions[2]
    micro_sd_bits = 200.e9  # use value for some large known mirco SD card size
    genome_volume = (0.34e-9 ** 3) * 6.e9  # maybe correct? http://hypertextbook.com/facts/1998/StevenChen.shtml
    genome_bits = 4.e6  # losslessly compressed aomunt of data
    micro_sd_density = micro_sd_bits / micro_sd_volume
    genome_density = genome_bits / genome_volume
    print(f'micro SD card volume: {micro_sd_volume:.3e} cubic meter')
    print(f'genome volume: {genome_volume:.3e} cubic meter')
    print(f'micro SD data per volume {micro_sd_density:.4e} bits/(cubic meter)')
    print(f'genome data per volume {genome_density:.4e} bits/(cubic meter)')


def main():
    np.set_printoptions(precision=2, suppress=True)
    data_per_volume()


if __name__ == '__main__':
    main()
