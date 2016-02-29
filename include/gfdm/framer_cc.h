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


#ifndef INCLUDED_GFDM_FRAMER_CC_H
#define INCLUDED_GFDM_FRAMER_CC_H

#include <gfdm/api.h>
#include <gfdm/preamble_generator.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief <+description of block+>
     * \ingroup gfdm
     *
     */
    class GFDM_API framer_cc : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<framer_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of gfdm::framer_cc.
       *
       * To avoid accidental use of raw pointers, gfdm::framer_cc's
       * constructor is in a private implementation
       * class. gfdm::framer_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(
          int nsubcarrier,
          int ntimeslots,
          bool sync,
          std::vector<gr_complex> sync_symbols,
          gr::gfdm::preamble_generator_sptr preamble_generator,
          const std::string& len_tag_key = "frame_len");
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_FRAMER_CC_H */

