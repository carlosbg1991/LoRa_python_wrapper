/* -*- c++ -*- */
/*
 * Copyright 2004,2010,2012 Free Software Foundation, Inc.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "fir_filter_ccf.h"
#include <volk/volk.h>

typedef std::complex<float>  gr_complex;

fir_filter_ccf::fir_filter_ccf(int decimation, const std::vector<float>& taps)
{
    d_align = volk_get_alignment();
    d_naligned = std::max((size_t)1, d_align / sizeof(gr_complex));

    d_aligned_taps = NULL;
    set_taps(taps);

    // Make sure the output sample is always aligned, too.
    d_output = (gr_complex*)volk_malloc(1 * sizeof(gr_complex), d_align);
}

fir_filter_ccf::~fir_filter_ccf()
{
    // Free all aligned taps
    if (d_aligned_taps != NULL) {
        for (int i = 0; i < d_naligned; i++) {
            volk_free(d_aligned_taps[i]);
        }
        ::free(d_aligned_taps);
        d_aligned_taps = NULL;
    }

    // Free output sample
    volk_free(d_output);
}

void fir_filter_ccf::set_taps(const std::vector<float>& taps)
{
    // Free the taps if already allocated
    if (d_aligned_taps != NULL) {
        for (int i = 0; i < d_naligned; i++) {
            volk_free(d_aligned_taps[i]);
        }
        ::free(d_aligned_taps);
        d_aligned_taps = NULL;
    }

    d_ntaps = (int)taps.size();
    d_taps = taps;
    std::reverse(d_taps.begin(), d_taps.end());

    // Make a set of taps at all possible arch alignments
    d_aligned_taps = (float**)malloc(d_naligned * sizeof(float*));
    for (int i = 0; i < d_naligned; i++) {
        d_aligned_taps[i] =
            (float*)volk_malloc((d_ntaps + d_naligned - 1) * sizeof(float), d_align);
        memset(d_aligned_taps[i], 0, sizeof(float) * (d_ntaps + d_naligned - 1));
        for (unsigned int j = 0; j < d_ntaps; j++)
            d_aligned_taps[i][i + j] = d_taps[j];
    }
}

void fir_filter_ccf::update_tap(float t, unsigned int index)
{
    d_taps[index] = t;
    for (int i = 0; i < d_naligned; i++) {
        d_aligned_taps[i][i + index] = t;
    }
}

std::vector<float> fir_filter_ccf::taps() const
{
    std::vector<float> t = d_taps;
    std::reverse(t.begin(), t.end());
    return t;
}

unsigned int fir_filter_ccf::ntaps() const { return d_ntaps; }

gr_complex fir_filter_ccf::filter(const gr_complex input[])
{
    const gr_complex* ar = (gr_complex*)((size_t)input & ~(d_align - 1));
    unsigned al = input - ar;

    volk_32fc_32f_dot_prod_32fc_a(d_output, ar, d_aligned_taps[al], (d_ntaps + al));
    return *d_output;
}

void fir_filter_ccf::filterN(gr_complex output[],
                             const gr_complex input[],
                             unsigned long n)
{
    for (unsigned long i = 0; i < n; i++)
        output[i] = filter(&input[i]);
}


void fir_filter_ccf::filterNdec(gr_complex output[],
                                const gr_complex input[],
                                unsigned long n,
                                unsigned int decimate)
{
    unsigned long j = 0;
    for (unsigned long i = 0; i < n; i++) {
        output[i] = filter(&input[j]);
        j += decimate;
    }
}