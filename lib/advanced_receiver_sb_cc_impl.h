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

#ifndef INCLUDED_GFDM_ADVANCED_RECEIVER_SB_CC_IMPL_H
#define INCLUDED_GFDM_ADVANCED_RECEIVER_SB_CC_IMPL_H

#include <gfdm/advanced_receiver_sb_cc.h>
#include <gfdm/receiver_kernel_cc.h>

namespace gr {
  namespace gfdm {

    class advanced_receiver_sb_cc_impl : public advanced_receiver_sb_cc
    {
     private:
       int d_n_subcarriers;
       int d_n_timeslots;
       int d_ic_iter;
       gr::digital::constellation_sptr d_constellation;
       void map_sc_symbols(std::vector< std::vector<gr_complex> > &sc_symbols);
       void map_symbols_to_constellation_points(gr_complex* symbols);
       void remove_sc_interference(
           std::vector< std::vector<gr_complex> > &sc_symbols,
           std::vector< std::vector<gr_complex> > &sc_fdomain);
       receiver_kernel_cc::sptr d_kernel;
       std::vector< std::vector<gr_complex> > d_sc_fdomain;
       std::vector< std::vector<gr_complex> > d_sc_symbols;

     public:
      advanced_receiver_sb_cc_impl(int n_timeslots, int n_subcarriers, int overlap, int ic_iter, std::vector< gr_complex > frequency_taps, gr::digital::constellation_sptr constellation);
      ~advanced_receiver_sb_cc_impl();
      void set_ic(int ic_iter){d_ic_iter = ic_iter;}

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_ADVANCED_RECEIVER_SB_CC_IMPL_H */

