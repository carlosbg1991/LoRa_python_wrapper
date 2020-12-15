/* -*- c++ -*- */
/*
 * Copyright 2004,2007,2010,2012-2013 Free Software Foundation, Inc.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "fractional_resampler.h"
#include <stdexcept>

fractional_resampler::fractional_resampler(float phase_shift, float resamp_ratio)
{
	std::cout << "[FRACTIONAL RESAMPLER] - Creating Fractional Resampler" << std::endl;

    d_mu = phase_shift;
    d_mu_inc = resamp_ratio;
    d_interp = new mmse_fir_interpolator();

    if (resamp_ratio <= 0)
        throw std::out_of_range("resampling ratio must be > 0");
    if (phase_shift < 0 || phase_shift > 1)
        throw std::out_of_range("phase shift ratio must be > 0 and < 1");

    std::cout << "[FRACTIONAL RESAMPLER] - To generate 1 output sample we require " << \
    (int) ceil( (1 * d_mu_inc) + d_interp->ntaps() ) <<
    " samples at the input."
    << std::endl;
}

fractional_resampler::~fractional_resampler() { 
    delete d_interp; 
}

void fractional_resampler::work(const gr_complex * in, gr_complex * out, int noutput_items)
{
    int ii = 0; // input index
    int oo = 0; // output index

    while (oo < noutput_items) {
        out[oo++] = d_interp->interpolate(&in[ii], d_mu);

        double s = d_mu + d_mu_inc;
        double f = floor(s);
        int incr = (int)f;
        d_mu = s - f;
        ii += incr;
    }
}

float fractional_resampler::mu() const { return d_mu; }

float fractional_resampler::resamp_ratio() const { return d_mu_inc; }

float fractional_resampler::get_ntaps() const { return d_interp->ntaps(); }

void fractional_resampler::set_mu(float mu) { d_mu = mu; }

void fractional_resampler::set_resamp_ratio(float resamp_ratio)
{
    d_mu_inc = resamp_ratio;
}