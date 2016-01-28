/* -*- c++ -*- */
/*
 * Copyright 2016 Andrej Rode
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_GFDM_RECEIVER_H
#define INCLUDED_GFDM_RECEIVER_H
#endif

#include <gfdm/gfdm_receiver.h>

namespace gr {
  namespace gfdm {
    namespace kernel {
          
      gfdm_receiver::gfdm_receiver()
      {
      }
     
      gfdm_receiver::~gfdm_receiver()
      {

      }
      
      void
      gfdm_receiver::filter_superposition(std::vector< std::vector<gr_complex> > &out,
          gr_complex &in)
      {

      }

      void
      gfdm_receiver::demodulate_subcarrier(gr_complex &out,
          std::vector< std::vector<gr_complex> > &sc_fdomain)
      {

      }

    } /* namespace kernel */
  } /* namespace filter */
} /* namespace gr */


