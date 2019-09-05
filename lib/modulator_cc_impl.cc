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
#include "modulator_cc_impl.h"

namespace gr {
  namespace gfdm {

    modulator_cc::sptr
    modulator_cc::make(
        int nsubcarrier,
        int ntimeslots,
    double filter_alpha,
        int fft_len,
        int sync_fft_len,
        const std::string& len_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new modulator_cc_impl(nsubcarrier,
                               ntimeslots,
                               filter_alpha,
                               fft_len,
                               sync_fft_len,
                               len_tag_key)
         );

    }

    /*
     * The private constructor
     */
    modulator_cc_impl::modulator_cc_impl(
        int nsubcarrier,
        int ntimeslots,
        double filter_alpha,
        int fft_len,
        int sync_fft_len,
        const std::string& len_tag_key)
      : gr::tagged_stream_block("modulator_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              len_tag_key),
      d_ntimeslots(ntimeslots),
      d_nsubcarrier(nsubcarrier),
      d_N(ntimeslots*nsubcarrier),
      d_fft_len(fft_len),
      d_sync_fft_len(sync_fft_len),
      d_len_tag_key(len_tag_key)
    {
      set_tag_propagation_policy(TPP_DONT);
      set_relative_rate(double(fft_len)/double(d_N));
      d_filter_width = 2;
      if (d_fft_len < d_N)
      {
        throw std::invalid_argument("fft_len must be greater than or equal to nsubcarrier*ntimeslots");
      }

      std::vector<gr_complex> filter_taps;
      rrc_filter_sparse *filter_gen = new rrc_filter_sparse(d_N,filter_alpha,d_filter_width,nsubcarrier,ntimeslots);
      filter_gen->get_taps(filter_taps);
      delete filter_gen;
      d_filter_taps = (gr_complex*) volk_malloc(filter_taps.size() * sizeof(gr_complex), volk_get_alignment());
      std::memcpy(d_filter_taps, &filter_taps[0], sizeof(gr_complex) * filter_taps.size());

      //Initialize FFT per subcarrier
      d_sc_fft = new fft::fft_complex(d_ntimeslots,true,1);
      d_sc_fft_in = d_sc_fft->get_inbuf();
      d_sc_fft_out = d_sc_fft->get_outbuf();

      //Initiailize sync FFTs in case there are sync symbols
      d_sync_ifft = new fft::fft_complex(d_sync_fft_len, false, 1);
      d_sync_ifft_in = d_sync_ifft->get_inbuf();
      d_sync_ifft_out = d_sync_ifft->get_outbuf();

      //Initialize resulting FFT
      d_out_ifft = new fft::fft_complex(d_fft_len,false,1);
      d_out_ifft_in = d_out_ifft->get_inbuf();
      d_out_ifft_out = d_out_ifft->get_outbuf();

      // holds intermediate data during frame modulation
      d_sc_tmp = (gr_complex*) volk_malloc(d_filter_width * d_ntimeslots * sizeof(gr_complex), volk_get_alignment());

    }

    /*
     * Our virtual destructor.
     */
    modulator_cc_impl::~modulator_cc_impl()
    {
      delete d_sc_fft;
      delete d_sync_ifft;
      delete d_out_ifft;
      volk_free(d_sc_tmp);
    }

    int
    modulator_cc_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      int noutput_items;
      if (ninput_items[0] == d_N)
      {
        noutput_items = d_fft_len;
      } else if (ninput_items[0] == d_N+d_sync_fft_len)
      {
        noutput_items = d_fft_len+d_sync_fft_len;
      } else
      {
        throw std::runtime_error("wrong number of input_items");
      }
      return noutput_items;
    }

    void
    modulator_cc_impl::update_length_tags(int n_produced, int n_ports)
    {
      return ;
    }

    void
    modulator_cc_impl::modulate_gfdm_frame(gr_complex *out, const gr_complex *in)
    {
      std::memset(d_sc_tmp, 0x00, d_filter_width * d_ntimeslots * sizeof(gr_complex));
      std::memset(d_out_ifft_in, 0x00, sizeof(gr_complex)*d_fft_len);
      for (int k = 0; k < d_nsubcarrier; k++) {
        // 1. FFT on subcarrier
        std::memcpy(d_sc_fft_in, in + k * d_ntimeslots, sizeof(gr_complex) * d_ntimeslots);
        d_sc_fft->execute();
        // 2. Multiply  with filtertaps (times filter_width)
        for (int l = 0; l < d_filter_width; l++) {
          ::volk_32fc_x2_multiply_32fc(d_sc_tmp + l * d_ntimeslots, d_sc_fft_out,
                                       d_filter_taps + l * d_ntimeslots, d_ntimeslots);
        }

        // 3. Add to ifft-vector
        // Calculate ifft offset (not shifted, possibly longer than symbols, some outofband radiation, per subcarrier offset
        int ifft_offset = (d_fft_len / 2 + (d_fft_len - d_N) / 2 - ((d_filter_width - 1) * (d_ntimeslots)) / 2 +
                           k * d_ntimeslots) % d_fft_len;

        // VOLKify me!
        for (int n = 0; n < d_filter_width * d_ntimeslots; n++) {
          d_out_ifft_in[((ifft_offset + n) % d_fft_len)] += d_sc_tmp[(n + (d_filter_width * d_ntimeslots) / 2) %
                                                                   (d_filter_width * d_ntimeslots)];
        }
      }
      d_out_ifft->execute();
      // scale result
      ::volk_32fc_s32fc_multiply_32fc(out, d_out_ifft_out, gr_complex(1.0 / d_N, 0), d_fft_len);
    }

    int
    modulator_cc_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];
        std::vector <gr::tag_t> tags;
        
        get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0)+ninput_items[0]);
        bool sync = false;
        uint64_t sync_offset = 0;
        uint64_t data_offset = 0;
        int sync_length = 0;
        uint64_t nsync_items = 0;
        for (std::vector<gr::tag_t>::iterator it = tags.begin() ; it!= tags.end(); ++it)
        {
          if (pmt::symbol_to_string(it->key) == "gfdm_sync")
          {
            sync = true;
            sync_offset = it->offset - nitems_read(0);
            nsync_items = pmt::to_uint64(it->value);

          }else if (pmt::symbol_to_string(it->key) == "gfdm_data")
          {
            data_offset = it->offset - nitems_read(0);       
          }
        }

        if (sync)
        {
        //do sync stuff
        sync_length = d_sync_fft_len;
        std::memcpy(&out[0],&in[sync_offset],sizeof(gr_complex)*nsync_items);
        add_item_tag(0, nitems_written(0),
            pmt::string_to_symbol(d_len_tag_key),
            pmt::from_long(nsync_items));
        }

//        std::cout << "sync_length = " << sync_length << ", offset = " << data_offset << std::endl;
        // This is where all the action really happens now!
        modulate_gfdm_frame(out + sync_length, in + data_offset);

        add_item_tag(0, nitems_written(0)+sync_length,
            pmt::string_to_symbol(d_len_tag_key),
            pmt::from_long(d_fft_len));
        return sync_length+d_fft_len;
    }

  } /* namespace gfdm */
} /* namespace gr */

