/* -*- c++ -*- */
/* 
 * Copyright 2016 Johannes Demel.
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


#ifndef INCLUDED_GFDM_FRAME_ENERGY_DETECTOR_CC_H
#define INCLUDED_GFDM_FRAME_ENERGY_DETECTOR_CC_H

#include <gfdm/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Detect frames based on energy ramp detection.
     * \ingroup gfdm
     *
     */
    class GFDM_API frame_energy_detector_cc : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<frame_energy_detector_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of gfdm::frame_energy_detector_cc.
       *
       * To avoid accidental use of raw pointers, gfdm::frame_energy_detector_cc's
       * constructor is in a private implementation
       * class. gfdm::frame_energy_detector_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(float alpha, int average_len, int frame_len, int backoff_len, const std::string& tag_key);

      virtual float alpha() = 0;
      virtual void set_alpha(float alpha) = 0;
      virtual int backoff_len() = 0;
      virtual void set_backoff_len(int backoff_len) = 0;
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_FRAME_ENERGY_DETECTOR_CC_H */

