/* -*- c++ -*- */
/*
 * Copyright 2002,2012 Free Software Foundation, Inc.
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

#include "interpolator_taps.h"
#include "mmse_fir_interpolator.h"
#include <stdexcept>

mmse_fir_interpolator::mmse_fir_interpolator()
{

    std::cout << "[MMSE FIR INTERPOLATOR] - Creating FIR filter" << std::endl;
    std::cout << "[MMSE FIR INTERPOLATOR] - NTAPS:" << NTAPS << std::endl;
    std::cout << "[MMSE FIR INTERPOLATOR] - NSTEPS:" << NSTEPS << std::endl;

    filters.resize(NSTEPS + 1);

    for (int i = 0; i < NSTEPS + 1; i++) {
        std::vector<float> t(&taps[i][0], &taps[i][NTAPS]);
        filters[i] = new fir_filter_ccf(1, t);
    }
}

mmse_fir_interpolator::~mmse_fir_interpolator()
{
    for (int i = 0; i < NSTEPS + 1; i++)
        delete filters[i];
}

unsigned mmse_fir_interpolator::ntaps() const { 
    return NTAPS;
}

unsigned mmse_fir_interpolator::nsteps() const { return NSTEPS;
}

gr_complex mmse_fir_interpolator::interpolate(const gr_complex input[], float mu) const
{
    int imu = (int)rint(mu * NSTEPS);

    if ((imu < 0) || (imu > NSTEPS)) {
        throw std::runtime_error("mmse_fir_interpolator: imu out of bounds.\n");
    }

    gr_complex r = filters[imu]->filter(input);
    return r;
}