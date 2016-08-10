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


#ifndef INCLUDED_GFDM_ADVANCED_RECEIVER_TSB_CC_H
#define INCLUDED_GFDM_ADVANCED_RECEIVER_TSB_CC_H

#include <gfdm/api.h>
#include <gnuradio/tagged_stream_block.h>
#include <gnuradio/digital/api.h>
#include <gnuradio/digital/constellation.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief GFDM-Demodulator with iterative interference cancellation implemented as tagged stream block.
     * \ingroup gfdm
     *
     */
    class GFDM_API advanced_receiver_tsb_cc : virtual public gr::tagged_stream_block
    {
     public:
      typedef boost::shared_ptr<advanced_receiver_tsb_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of gfdm::advanced_receiver_tsb_cc.
       *
       * To avoid accidental use of raw pointers, gfdm::advanced_receiver_tsb_cc's
       * constructor is in a private implementation
       * class. gfdm::advanced_receiver_tsb_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(int n_timeslots, int n_subcarriers, int overlap, int ic_iter, std::vector< gr_complex > frequency_taps, gr::digital::constellation_sptr constellation, const std::string& len_tag_key = "gfdm_block");
      virtual void set_ic(int ic_iter){};
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_ADVANCED_RECEIVER_TSB_CC_H */

