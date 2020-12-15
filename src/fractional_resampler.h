/* -*- c++ -*- */
/*
 * Copyright 2004,2007,2012-2013 Free Software Foundation, Inc.
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

#ifndef FRACTIONAL_RESAMPLER_H
#define FRACTIONAL_RESAMPLER_H

#include "mmse_fir_interpolator.h"
#include <iostream>

class fractional_resampler
{
private:
    float d_mu;
    float d_mu_inc;
    mmse_fir_interpolator* d_interp;

public:
    // fractional_resampler();
    fractional_resampler(float phase_shift, float resamp_ratio);
    ~fractional_resampler();

    void work(const gr_complex * input, gr_complex * output, int ninput_items);

    float mu() const;
    float resamp_ratio() const;
    float get_ntaps() const;
    void set_mu(float mu);
    void set_resamp_ratio(float resamp_ratio);
};

#endif /* FRACTIONAL_RESAMPLER_H */
