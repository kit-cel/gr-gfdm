/* -*- c++ -*- */
/* 
 * Copyright 2016 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#include <gfdm/resource_demapper_kernel_cc.h>
#include <string.h>
#include <stdexcept>

namespace gr {
  namespace gfdm {

    resource_demapper_kernel_cc::resource_demapper_kernel_cc(int timeslots, int subcarriers, int active_subcarriers, std::vector<int> subcarrier_map, bool per_timeslot)
    {
    }

    resource_demapper_kernel_cc::~resource_demapper_kernel_cc()
    {
    }

  } /* namespace gfdm */
} /* namespace gr */

