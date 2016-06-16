/* -*- c++ -*- */
/* 
 * Copyright 2016 Andrej Rode.
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


#ifndef INCLUDED_GFDM_CYCLIC_PREFIXER_CC_H
#define INCLUDED_GFDM_CYCLIC_PREFIXER_CC_H

#include <gfdm/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Add Cyclic Prefix to GFDM block and apply block pinching (W-GFDM).
     * \ingroup gfdm
     *
     */
    class GFDM_API cyclic_prefixer_cc : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<cyclic_prefixer_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of gfdm::cyclic_prefixer_cc.
       *
       * To avoid accidental use of raw pointers, gfdm::cyclic_prefixer_cc's
       * constructor is in a private implementation
       * class. gfdm::cyclic_prefixer_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(int cp_length, int ramp_len, int block_len, std::vector<gr_complex> window_taps);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_CYCLIC_PREFIXER_CC_H */

