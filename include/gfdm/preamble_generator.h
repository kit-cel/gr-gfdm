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


#ifndef INCLUDED_GFDM_PREAMBLE_GENERATOR_H
#define INCLUDED_GFDM_PREAMBLE_GENERATOR_H

#include <gfdm/api.h>

namespace gr {
  namespace gfdm {

    /*!
     * \brief <+description+>
     *
     */
    class GFDM_API preamble_generator
    {
    public:
      preamble_generator(int nsubcarrier, double filter_alpha);
      ~preamble_generator();
    private:
    };

  } // namespace gfdm
} // namespace gr

#endif /* INCLUDED_GFDM_PREAMBLE_GENERATOR_H */

