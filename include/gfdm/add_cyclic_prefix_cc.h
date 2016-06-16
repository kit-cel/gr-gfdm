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


#ifndef INCLUDED_GFDM_ADD_CYCLIC_PREFIX_CC_H
#define INCLUDED_GFDM_ADD_CYCLIC_PREFIX_CC_H

#include <complex>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Kernel adds cyclic prefix to GFDM frame and applies block pinching window.
     *
     */
//    class GFDM_API add_cyclic_prefix_cc
    class add_cyclic_prefix_cc
    {
    public:
      typedef std::complex<float> gfdm_complex;
      typedef boost::shared_ptr<add_cyclic_prefix_cc> sptr;

      add_cyclic_prefix_cc(int ramp_len, int cp_len, int block_len, std::vector<gfdm_complex> window_taps);
      ~add_cyclic_prefix_cc();
      void generic_work(gfdm_complex* p_out, const gfdm_complex* p_in);
      int block_size(){ return d_block_len;};
      int frame_size(){ return block_size() + d_cp_len;};
    private:
      int d_ramp_len;
      int d_cp_len;
      int d_block_len;
      gfdm_complex* d_front_ramp;
      gfdm_complex* d_back_ramp;
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_ADD_CYCLIC_PREFIX_CC_H */

