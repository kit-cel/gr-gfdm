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

#ifndef INCLUDED_GFDM_RECEIVER_CC_IMPL_H
#define INCLUDED_GFDM_RECEIVER_CC_IMPL_H

#include <gfdm/receiver_cc.h>
#include <gnuradio/fft/fft.h>
#include <volk/volk.h>
#include <gnuradio/filter/firdes.h>
#include <gfdm/receiver_kernel_cc.h>

namespace gr {
  namespace gfdm {

    class receiver2_cc_impl :  public receiver_cc
    {
     private:
      receiver_kernel_cc::sptr d_kernel;
     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      receiver2_cc_impl(
          int nsubcarrier,
          int ntimeslots,
          double filter_alpha, 
          int fft_len,
          const std::string& len_tag_key);
      ~receiver2_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_RECEIVER_CC_IMPL_H */

