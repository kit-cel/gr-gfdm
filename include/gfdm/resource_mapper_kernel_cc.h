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


#ifndef INCLUDED_GFDM_RESOURCE_MAPPER_KERNEL_CC_H
#define INCLUDED_GFDM_RESOURCE_MAPPER_KERNEL_CC_H

#include <complex>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace gr {
  namespace gfdm {

    /*!
     * \brief map complex information symbols to GFDM resource grid.
     * Input is a vector with all complex information symbols for one GFDM frame.
     * Result is a vector which is fed to gfdm_modulator.
     *
     */
    class resource_mapper_kernel_cc
    {
    public:
      typedef std::complex<float> gfdm_complex;
      typedef boost::shared_ptr<resource_mapper_kernel_cc> sptr;

      resource_mapper_kernel_cc(int timeslots, int subcarriers, int active_subcarriers, std::vector<int> subcarrier_map, bool per_timeslot = true);
      ~resource_mapper_kernel_cc();
      int input_vector_size(){ return d_active_subcarriers * d_timeslots;};
      int output_vector_size(){ return d_subcarriers * d_timeslots;};
      void generic_work(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size);
    private:
      int d_timeslots;
      int d_subcarriers;
      int d_active_subcarriers;
      std::vector<int> d_subcarrier_map;
      bool d_per_timeslot;

      void map_per_timeslot(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size);
      void map_per_subcarrier(gfdm_complex* p_out, const gfdm_complex* p_in, const int ninput_size);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_RESOURCE_MAPPER_KERNEL_CC_H */

