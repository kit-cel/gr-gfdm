/* -*- c++ -*- */
/* 
 * Copyright 2016 Andrej Rode, Johannes Demel.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "advanced_receiver_sb_cc_impl.h"
#include <volk/volk.h>

namespace gr
{
  namespace gfdm
  {

    advanced_receiver_sb_cc::sptr
    advanced_receiver_sb_cc::make(int n_timeslots, int n_subcarriers, int overlap, int ic_iter,
                                  std::vector<gr_complex> frequency_taps, gr::digital::constellation_sptr constellation,
                                  std::vector<int> subcarrier_map)
    {
      return gnuradio::get_initial_sptr
              (new advanced_receiver_sb_cc_impl(n_timeslots, n_subcarriers, overlap, ic_iter, frequency_taps,
                                                constellation, subcarrier_map));
    }

    /*
     * The private constructor
     */
    advanced_receiver_sb_cc_impl::advanced_receiver_sb_cc_impl(int n_timeslots, int n_subcarriers, int overlap,
                                                               int ic_iter, std::vector<gr_complex> frequency_taps,
                                                               gr::digital::constellation_sptr constellation,
                                                               std::vector<int> subcarrier_map)
            : gr::sync_block("advanced_receiver_sb_cc",
                             gr::io_signature::make(1, 2, sizeof(gr_complex)),
                             gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
      d_adv_kernel = advanced_receiver_kernel_cc::sptr(new advanced_receiver_kernel_cc(n_timeslots, n_subcarriers, overlap, frequency_taps, subcarrier_map, ic_iter, constellation));
      set_output_multiple(d_adv_kernel->block_size());
    }

    /*
     * Our virtual destructor.
     */
    advanced_receiver_sb_cc_impl::~advanced_receiver_sb_cc_impl()
    {
    }

    int
    advanced_receiver_sb_cc_impl::work(int noutput_items,
                                       gr_vector_const_void_star &input_items,
                                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      const int n_blocks = noutput_items / d_adv_kernel->block_size();

      if(input_items.size() > 1){
        std::cout << "use EQ\n";
        const gr_complex *in_eq = (const gr_complex *) input_items[1];
        for (int i = 0; i < n_blocks; ++i) {
          d_adv_kernel->generic_work_equalize(out, in, in_eq);

          in += d_adv_kernel->block_size();
          in_eq += d_adv_kernel->block_size();
          out += d_adv_kernel->block_size();
        }
      }
      else{
        for (int i = 0; i < n_blocks; ++i) {
          d_adv_kernel->generic_work(out, in);

          in += d_adv_kernel->block_size();
          out += d_adv_kernel->block_size();
        }
      }

      return n_blocks * d_adv_kernel->block_size();
    }

  } /* namespace gfdm */
} /* namespace gr */

