/* -*- c++ -*- */
/*  
 * Copyright 2016 Andrej Rode
 *  
 * This file is part of GNU Radio
 *  
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *  
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */ 

#ifndef INCLUDED_GFDM_UTILS_H
#define INCLUDED_GFDM_UTILS_H
#endif

#include <gnuradio/filter/firdes.h>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace gfdm {

      class rrc_filter_sparse {
        private:
          std::vector<gr_complex> d_filter_taps;

        public:
          rrc_filter_sparse(int ntaps, double alpha, int filter_width, int nsubcarrier, int ntimeslots);
          void get_taps(std::vector<gr_complex> &out);
          ~rrc_filter_sparse();
      };

  } /* namespace gfdm */
} /* namespace gr */
