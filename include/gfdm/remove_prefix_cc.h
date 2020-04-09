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


#ifndef INCLUDED_GFDM_REMOVE_PREFIX_CC_H
#define INCLUDED_GFDM_REMOVE_PREFIX_CC_H

#include <gnuradio/block.h>
#include <gfdm/api.h>

namespace gr {
namespace gfdm {

/*!
 * \brief extract block_len items from frame_len chunks of items, marked with a tag plus
 * offset \ingroup gfdm
 *
 */
class GFDM_API remove_prefix_cc : virtual public gr::block
{
public:
    typedef boost::shared_ptr<remove_prefix_cc> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of gfdm::remove_prefix_cc.
     *
     * To avoid accidental use of raw pointers, gfdm::remove_prefix_cc'
     * constructor is in a private implementation
     * class. gfdm::remove_prefix_cc::make is the public interface for
     * creating new instances.
     */
    static sptr
    make(int frame_len, int block_len, int offset, const std::string& gfdm_sync_tag_key);
};

} // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_REMOVE_PREFIX_CC_H */
