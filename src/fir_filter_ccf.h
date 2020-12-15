/* -*- c++ -*- */
/*
 * Copyright 2004,2012 Free Software Foundation, Inc.
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

/* WARNING: this file is machine generated. Edits will be overwritten */

#ifndef FILTER_FIR_FILTER_CCF_H
#define	FILTER_FIR_FILTER_CCF_H

#include <vector>
#include <complex>
#include <fstream>
#include <stdint.h>
#include <algorithm>
#include <cstring>

typedef std::complex<float>  gr_complex;

class fir_filter_ccf
{
public:
    fir_filter_ccf(int decimation, const std::vector<float>& taps);
    ~fir_filter_ccf();

    void set_taps(const std::vector<float>& taps);
    void update_tap(float t, unsigned int index);
    std::vector<float> taps() const;
    unsigned int ntaps() const;

    gr_complex filter(const gr_complex input[]);
    void filterN(gr_complex output[], const gr_complex input[], unsigned long n);
    void filterNdec(gr_complex output[],
                    const gr_complex input[],
                    unsigned long n,
                    unsigned int decimate);

protected:
    std::vector<float> d_taps;
    unsigned int d_ntaps;
    float** d_aligned_taps;
    gr_complex* d_output;
    int d_align;
    int d_naligned;
};

#endif /* FILTER_FIR_FILTER_CCF_H */
