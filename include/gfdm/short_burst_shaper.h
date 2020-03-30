/* -*- c++ -*- */
/*
 * Copyright 2019 Johannes Demel.
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


#ifndef INCLUDED_GFDM_SHORT_BURST_SHAPER_H
#define INCLUDED_GFDM_SHORT_BURST_SHAPER_H

#include <gfdm/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Add Padding and Scale Burst.
     * \ingroup gfdm
     *
     */
    class GFDM_API short_burst_shaper : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<short_burst_shaper> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of gfdm::short_burst_shaper.
       *
       * To avoid accidental use of raw pointers, gfdm::short_burst_shaper's
       * constructor is in a private implementation
       * class. gfdm::short_burst_shaper::make is the public interface for
       * creating new instances.
       */
      static sptr make(int pre_padding, int post_padding, gr_complex scale,
                       const std::string &length_tag_name="packet_len");

      /*!
      * \brief Return multiplicative constant
      */
      virtual gr_complex scale() const = 0;

      /*!
      * \brief Set multiplicative constant
      */
      virtual void set_scale(gr_complex scale) = 0;
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_SHORT_BURST_SHAPER_H */

