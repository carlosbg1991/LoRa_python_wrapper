/* -*- c++ -*- */
/*
 * Copyright 2017 Pieter Robyns.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifndef CHANNELIZER_H
#define CHANNELIZER_H

#include <string>
#include <vector>
#include <fstream>
#include <stdint.h>
#include <stdexcept>

// #include <gnuradio/filter/freq_xlating_fir_filter_ccf.h>
#include "freq_xlating_fir_filter_ccf.h"
#include "fractional_resampler.h"
#include "firdes.h"

class channelizer
{
 private:
     // gr::filter::freq_xlating_fir_filter_ccf::sptr d_xlating_fir_filter;
     // fractional_resampler d_resampler = new fractional_resampler(10e6, 1e6);
     freq_xlating_fir_filter_ccf * d_xlating_fir_filter;
     fractional_resampler * d_resampler;
     std::vector<float> d_lpf;
     float d_cfo;
     uint32_t d_freq_offset;
     int d_resampler_in_out_ratio;
     gr_complex * bridge_resampler_xlatingFIR;
     int d_batch_size;

 public:
  channelizer(float in_samp_rate, float out_samp_rate, float center_freq, std::vector<float> channel_list);
  ~channelizer();
  void work_resampler(const gr_complex * in, gr_complex * out, int noutput_items);
  void work_freq_xlating_fir_filter(const gr_complex * in, gr_complex * out, int noutput_items);
  void work(const gr_complex * in, gr_complex * out, int noutput_items);
  int get_resampler_in_out_ratio();
  int get_resampler_nTaps();
  void set_batch_size(int batch_size);

  // Where all the action really happens
};

#endif /* CHANNELIZER_H */
