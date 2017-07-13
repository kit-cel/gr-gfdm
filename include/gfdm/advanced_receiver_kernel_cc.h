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


#ifndef INCLUDED_GFDM_ADVANCED_RECEIVER_KERNEL_CC_H
#define INCLUDED_GFDM_ADVANCED_RECEIVER_KERNEL_CC_H

#include <gfdm/api.h>
#include <gfdm/receiver_kernel_cc.h>
#include <gnuradio/digital/constellation.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief Hold config and functions for advanced IC kernel.
     *
     */
    class GFDM_API advanced_receiver_kernel_cc
    {
    public:
      typedef boost::shared_ptr<advanced_receiver_kernel_cc> sptr;
      advanced_receiver_kernel_cc(int timeslots, int subcarriers, int overlap, std::vector<gr_complex> frequency_taps,
                                  std::vector<int> subcarrier_map, int ic_iter,
                                  gr::digital::constellation_sptr constellation);
      ~advanced_receiver_kernel_cc();

      void generic_work(gr_complex *p_out, const gr_complex *p_in);
      void generic_work_equalize(gr_complex *out, const gr_complex *in, const gr_complex* f_eq_in);
      void set_ic(int ic_iter) { d_ic_iter = ic_iter; }
      int block_size() { return d_kernel->block_size();}

    private:
      receiver_kernel_cc::sptr d_kernel;
      std::vector<int> d_subcarrier_map;
      int d_ic_iter;
      gr::digital::constellation_sptr d_constellation;

      void map_symbols_to_constellation_points(gr_complex *p_out, const gr_complex *p_in);
      void perform_ic_iterations(gr_complex *p_out, const gr_complex *p_freq_block);

      gr_complex *d_freq_block;
      gr_complex *d_ic_time_buffer;
      gr_complex *d_ic_freq_buffer;


    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_ADVANCED_RECEIVER_KERNEL_CC_H */

