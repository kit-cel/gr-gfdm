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


#ifndef INCLUDED_GFDM_DETECT_FRAME_ENERGY_KERNEL_CL_H
#define INCLUDED_GFDM_DETECT_FRAME_ENERGY_KERNEL_CL_H

#include <complex>
#include <boost/shared_ptr.hpp>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Perform rough energy based synchronization for TDD
     * Calculate Energy over average_len samples and put out flag if previous_energy < alpha * current_energy.
     * Flag is returned for the first such block.
     *
     */
    class detect_frame_energy_kernel_cl
    {
    public:
      typedef std::complex<float> gfdm_complex;
      typedef boost::shared_ptr<detect_frame_energy_kernel_cl> sptr;

      detect_frame_energy_kernel_cl(float alpha, int average_len);
      ~detect_frame_energy_kernel_cl();
      long detect_frame(const gfdm_complex* p_in, const int ninput_items);
      int average_len(){ return d_average_len;};
      float alpha(){ return d_alpha;};
      void set_alpha(float alpha){d_alpha = alpha;};
      void set_initial_energy_value(const float e){d_block_energy = e;};
    private:
      float d_alpha;
      int d_average_len;

      float d_block_energy;

    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_DETECT_FRAME_ENERGY_KERNEL_CL_H */

