/* -*- c++ -*- */
/*
 * Copyright 2019 Johannes Demel.
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

#ifndef INCLUDED_GFDM_SHORT_BURST_SHAPER_IMPL_H
#define INCLUDED_GFDM_SHORT_BURST_SHAPER_IMPL_H

#include <gfdm/short_burst_shaper.h>

namespace gr {
  namespace gfdm {

    class short_burst_shaper_impl : public short_burst_shaper
    {
     private:
      const int d_pre_padding;
      const int d_post_padding;
      gr_complex d_scale;

     protected:
      int calculate_output_stream_length(const gr_vector_int &ninput_items);

     public:
      short_burst_shaper_impl(int pre_padding, int post_padding, gr_complex scale,
                              const std::string &length_tag_name);
      ~short_burst_shaper_impl();

      gr_complex scale() const { return d_scale; }
      void set_scale(gr_complex scale) { d_scale = scale; }

      // Where all the action really happens
      int work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_SHORT_BURST_SHAPER_IMPL_H */

