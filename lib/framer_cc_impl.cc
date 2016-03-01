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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "framer_cc_impl.h"

namespace gr {
  namespace gfdm {

    framer_cc::sptr
    framer_cc::make(
        int nsubcarrier,
        int ntimeslots,
        bool sync,
        std::vector<gr_complex> sync_symbols,
        gr::gfdm::preamble_generator_sptr preamble_generator,
        const std::string& len_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new framer_cc_impl(nsubcarrier,
                            ntimeslots,
                            sync,
                            sync_symbols,
                            preamble_generator,
                            len_tag_key));
    }

    /*
     * The private constructor
     */
    framer_cc_impl::framer_cc_impl(
        int nsubcarrier,
        int ntimeslots,
        bool sync,
        std::vector<gr_complex> sync_symbols,
        gr::gfdm::preamble_generator_sptr preamble_generator,
        const std::string& len_tag_key)
      : gr::tagged_stream_block("framer_cc",
              gr::io_signature::make(1,1, sizeof(gr_complex)),
              gr::io_signature::make(1,1, sizeof(gr_complex)),
              len_tag_key),
      d_nsubcarrier(nsubcarrier),
      d_ntimeslots(ntimeslots),
      d_preamble_generator(preamble_generator),
      d_sync(sync)
    {
      if (d_sync)
      {
        if (d_preamble_generator){
          //d_sync_symbols.resize(d_preamble_generator->get_preamble_len());
          d_sync_symbols = d_preamble_generator->get_preamble();
          //std::memcpy(&d_sync_symbols[0],&d_preamble_generator->get_preamble()[0],sizeof(gr_complex)*d_preamble_generator->get_preamble_len());
        //}
        //else if (sync_symbols.size() < d_nsubcarrier)
        //{
        //  throw std::invalid_argument("number of sync symbols must be equal to or greater than nsubcarrier");
        }else
        {
          d_sync_symbols.resize(2*d_nsubcarrier);
          std::memcpy(&d_sync_symbols[0],&sync_symbols[0],sizeof(gr_complex)*2*nsubcarrier);
        }
      }
    
    }

    /*
     * Our virtual destructor.
     */
    framer_cc_impl::~framer_cc_impl()
    {
    }

    int
    framer_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items = ninput_items[0]; 
      if (d_sync)
      {
        noutput_items = ninput_items[0]+(d_nsubcarrier*2); 
      }
      return noutput_items;
    }

    int
    framer_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0]; 
        
        int sync_offset = 0;
        if (d_sync)
        {
          sync_offset = d_sync_symbols.size();
//          for (int i=0; i<d_nsubcarrier; i++)
//          {
//            out[2*i+1] = d_sync_symbols[i];
//            out[2*i] = d_sync_symbols[i];
//          }
          std::memcpy(&out[0],&d_sync_symbols[0],sizeof(gr_complex)*d_sync_symbols.size());
          add_item_tag(0, nitems_written(0),
              pmt::string_to_symbol("gfdm_sync"),
              pmt::from_uint64(d_sync_symbols.size()));
        }
        add_item_tag(0, nitems_written(0)+sync_offset,
              pmt::string_to_symbol("gfdm_data"),
              pmt::from_uint64(d_ntimeslots*d_nsubcarrier));
        for (int k=0; k<d_nsubcarrier; k++)
        {
          for(int m=0; m<d_ntimeslots; m++)
          {
            out[(k*d_ntimeslots)+m+sync_offset] = in[(m*d_nsubcarrier+k)];
            
          }
        }
        int new_noutput_items = d_nsubcarrier*d_ntimeslots+sync_offset;

        return new_noutput_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

