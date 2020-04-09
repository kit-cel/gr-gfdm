/* -*- c++ -*- */
/*
 * Copyright 2018 Johannes Demel.
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


#ifndef INCLUDED_GFDM_TRANSMITTER_CC_H
#define INCLUDED_GFDM_TRANSMITTER_CC_H

#include <gnuradio/block.h>
#include <gfdm/api.h>

namespace gr {
namespace gfdm {

/*!
 * \brief <+description of block+>
 * \ingroup gfdm
 *
 */
class GFDM_API transmitter_cc : virtual public gr::block
{
public:
    typedef boost::shared_ptr<transmitter_cc> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of gfdm::transmitter_cc.
     *
     * To avoid accidental use of raw pointers, gfdm::transmitter_cc's
     * constructor is in a private implementation
     * class. gfdm::transmitter_cc::make is the public interface for
     * creating new instances.
     */
    static sptr make(int timeslots,
                     int subcarriers,
                     int active_subcarriers,
                     int cp_len,
                     int cs_len,
                     int ramp_len,
                     std::vector<int> subcarrier_map,
                     bool per_timeslot,
                     int overlap,
                     std::vector<gr_complex> frequency_taps,
                     std::vector<gr_complex> window_taps,
                     std::vector<gr_complex> preamble,
                     const std::string& tsb_tag_key = "");
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_TRANSMITTER_CC_H */
