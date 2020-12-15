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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "channelizer.h"
#include <iostream>

/*
 * The private constructor
 */
channelizer::channelizer(float in_samp_rate, float out_samp_rate, float center_freq, std::vector<float> channel_list)
{
	std::cout << "[CHANNELIZER] - Creating Channelizer" << std::endl;
    d_cfo = 0.0;
    d_lpf = firdes::low_pass(1.0, out_samp_rate, 62500+15000, 10000, firdes::WIN_HAMMING, 6.67);
    d_freq_offset = channel_list[0] - center_freq;
    d_xlating_fir_filter = new freq_xlating_fir_filter_ccf(1, d_lpf, d_freq_offset, out_samp_rate);
    d_resampler = new fractional_resampler(0, (float)in_samp_rate / (float)out_samp_rate);

    d_resampler_in_out_ratio = ceil(d_resampler->resamp_ratio());

    std::cout << "[CHANNELIZER] - Expected ratio IN/OUT: " << d_resampler_in_out_ratio << std::endl;

    // pre-allocate memory to inter-mediate variable pointer
    bridge_resampler_xlatingFIR = (gr_complex*) malloc(500*sizeof(gr_complex));

    // connect(self(), 0, d_resampler, 0);
    // connect(d_resampler, 0, d_xlating_fir_filter, 0);
    // connect(d_xlating_fir_filter, 0, self(), 0);
}

/*
 * Our virtual destructor.
 */
channelizer::~channelizer() { }

int channelizer::get_resampler_in_out_ratio() { return d_resampler_in_out_ratio; }
int channelizer::get_resampler_nTaps() { return d_resampler->get_ntaps(); }

void channelizer::set_batch_size(int batch_size) 
{ 
	d_batch_size = batch_size; 
	// Allocate memory
	free(bridge_resampler_xlatingFIR);
	bridge_resampler_xlatingFIR = (gr_complex*) malloc(d_batch_size*sizeof(gr_complex));
}

void channelizer::work_resampler(const gr_complex * in, gr_complex * out, int noutput_items)
{
	if (noutput_items<d_resampler_in_out_ratio)
		throw std::out_of_range("The number of inputs is insufficient to generate a single output sample");
	d_resampler->work(in, out, noutput_items);
}

void channelizer::work_freq_xlating_fir_filter(const gr_complex * in, gr_complex * out, int noutput_items)
{
	d_xlating_fir_filter->work(in, out, noutput_items);
}

void channelizer::work(const gr_complex * in, gr_complex * out, int noutput_items)
{
	if (noutput_items<d_resampler_in_out_ratio)
		throw std::out_of_range("The number of inputs is insufficient to generate a single output sample");

	const gr_complex * in_xlating;
	d_resampler->work(in, bridge_resampler_xlatingFIR, noutput_items);

	// cast into const. This assigns the memory location where the first pointer points to, so no 
	// need for memory pre-allocation for in_xlating
	in_xlating = const_cast<const gr_complex *> (bridge_resampler_xlatingFIR); 

	d_xlating_fir_filter->work(in_xlating, out, noutput_items);
}