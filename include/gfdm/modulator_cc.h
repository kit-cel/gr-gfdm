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


#ifndef INCLUDED_GFDM_MODULATOR_CC_H
#define INCLUDED_GFDM_MODULATOR_CC_H

#include <gnuradio/tagged_stream_block.h>
#include <gfdm/api.h>

namespace gr {
namespace gfdm {

/*!
 * \brief modulator_cc modulates an input tagged input stream with subcarrier-wise sorting
 * according to "Generalized Frequency Division Multiplexing: Analysis of an Alternative
 * Multi-Carrier Technique for Next Generation Cellular Systems" by N. Michailow
 * Low Complexity, Sparse Frequency Domain Processing
 * \ingroup gfdm
 *
 */
class GFDM_API modulator_cc : virtual public gr::tagged_stream_block
{
public:
    typedef std::shared_ptr<modulator_cc> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of gfdm::modulator_cc.
     *
     * To avoid accidental use of raw pointers, gfdm::modulator_cc's
     * constructor is in a private implementation
     * class. gfdm::modulator_cc::make is the public interface for
     * creating new instances.
     */
    static sptr make(int nsubcarrier,
                     int ntimeslots,
                     double filter_alpha,
                     int fft_len,
                     int sync_fft_len,
                     const std::string& len_tag_key = "frame_len");
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_MODULATOR_CC_H */
