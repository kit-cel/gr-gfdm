/* -*- c++ -*- */
/*
 * Copyright 2017 <+YOU OR YOUR COMPANY+>.
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
#include "extract_burst_cc_impl.h"
#include <cmath>
#include <volk/volk.h>

namespace gr {
  namespace gfdm {

    extract_burst_cc::sptr
    extract_burst_cc::make(int burst_len, int tag_backoff, std::string burst_start_tag,
                           bool activate_cfo_correction)
    {
      return gnuradio::get_initial_sptr
        (new extract_burst_cc_impl(burst_len, tag_backoff, burst_start_tag,
                                   activate_cfo_correction));
    }

    /*
     * The private constructor
     */
    extract_burst_cc_impl::extract_burst_cc_impl(int burst_len, int tag_backoff,
                                                 std::string burst_start_tag,
                                                 bool activate_cfo_correction)
      : gr::block("extract_burst_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
              d_burst_len(burst_len), d_tag_backoff(tag_backoff),
              d_activate_cfo_correction(activate_cfo_correction)
    {
      set_output_multiple(burst_len);
      d_burst_start_tag = pmt::string_to_symbol(burst_start_tag);
      set_tag_propagation_policy(TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    extract_burst_cc_impl::~extract_burst_cc_impl()
    {
    }

    void
    extract_burst_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
        ninput_items_required[0] = noutput_items;
    }

    float
    extract_burst_cc_impl::get_scale_factor(pmt::pmt_t info){
      float scale_factor = 1.0f;
      if(pmt::is_dict(info)){
        pmt::pmt_t scl = pmt::dict_ref(info,
                                      pmt::mp("scale_factor"),
                                      pmt::PMT_NIL);
        if(pmt::is_real(scl)){
          scale_factor = pmt::to_double(scl);
        }
      }
      return scale_factor;
    }

    gr_complex
    extract_burst_cc_impl::get_phase_rotation(pmt::pmt_t info){
      gr_complex fq_comp_rot= 1;
      pmt::pmt_t sc_rot= pmt::dict_ref(info,
                                       pmt::mp("sc_rot"),
                                       pmt::PMT_NIL);

      if(pmt::is_complex(sc_rot)) {
        fq_comp_rot= std::conj(pmt::to_complex(sc_rot));
        fq_comp_rot/= std::abs(fq_comp_rot);
      }
      return fq_comp_rot;
    }

    void
    extract_burst_cc_impl::normalize_power_level(gr_complex* p_out, const gr_complex* p_in, const float norm_factor, const int ninput_size)
    {
      volk_32f_s32f_multiply_32f((float*) p_out, (const float*) p_in,
                                 norm_factor, 2 * ninput_size);
    }

    void
    extract_burst_cc_impl::compensate_cfo(gr_complex* p_out, const gr_complex* p_in, const gr_complex phase_increment, const int ninput_size)
    {
      gr_complex initial_phase = gr_complex(1.0f, 0.0f);
      volk_32fc_s32fc_x2_rotator_32fc(p_out, p_in,
                                      phase_increment,
                                      &initial_phase,
                                      ninput_size);
    }


    int
    extract_burst_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      int n_out_bursts = noutput_items / d_burst_len;
      int avail_items = ninput_items[0];
      int consumed_items = avail_items;
      int produced_items = 0;

      std::vector<tag_t> tags;
      get_tags_in_window(tags, 0, 0, avail_items, d_burst_start_tag);
      const int n_max_bursts = std::min(int(tags.size()), n_out_bursts);
      for (int i = 0; i < n_max_bursts; ++i) {
        int burst_start = tags[i].offset - nitems_read(0);
        int actual_start = burst_start - d_tag_backoff;

        pmt::pmt_t info = tags[i].value;
        const float scale_factor = get_scale_factor(info);

        if(avail_items - burst_start >= d_burst_len){

          if(actual_start < 0){
            int num_prepend_zeros = std::abs(actual_start);
            memset(out, 0, sizeof(gr_complex) * num_prepend_zeros);
            // memcpy(out + num_prepend_zeros, in, sizeof(gr_complex) * d_burst_len);
            normalize_power_level(out + num_prepend_zeros, in, scale_factor, d_burst_len);
          }
          else{
            // memcpy(out, in + actual_start, sizeof(gr_complex) * d_burst_len);
            normalize_power_level(out, in + actual_start, scale_factor, d_burst_len);
          }

          if (d_activate_cfo_correction){
            gr_complex fq_comp_rot= get_phase_rotation(info);
            compensate_cfo(out, out, fq_comp_rot, d_burst_len);
          }

          add_item_tag(0, nitems_written(0) + produced_items, d_burst_start_tag,
                       tags[i].value, pmt::string_to_symbol(name()));

          produced_items += d_burst_len;
          consumed_items = burst_start + d_burst_len;
          out += d_burst_len;
        }
        else{
          consumed_items = std::max(0, actual_start);
          break;
        }
      }

      consume_each (consumed_items);
      return produced_items;
    }

  } /* namespace gfdm */
} /* namespace gr */

