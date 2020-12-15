/*
 * decoder.h
 *
 *  Created on: Nov 17, 2020
 *      Author: Carlos Bocanegra
 */

#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "decoder.h"

#include <liquid/liquid.h>
#include "loratap.h"
#include "utilities.h"
#include "tables.h"

typedef std::complex<float>  gr_complex;

decoder::decoder(float samp_rate, int sf)
// decoder::decoder()
{
    std::cout << "[DECODER] - Creating Decoder"<< std::endl;
    d_bw                 = 125000u;
    d_phdr.cr            = 4;
    d_samples_per_second = samp_rate;
    d_payload_symbols    = 0;
    d_cfo_estimation     = 0.0f;
    d_dt                 = 1.0f / d_samples_per_second;
    d_sf                 = sf;
    d_bits_per_second    = (double)d_sf * (double)(4.0 / (4.0 + d_phdr.cr)) / (1u << d_sf) * d_bw;
    d_symbols_per_second = (double)d_bw / (1u << d_sf);
    d_period             = 1.0f / (double)d_symbols_per_second;
    d_bits_per_symbol    = (double)(d_bits_per_second    / d_symbols_per_second);
    d_samples_per_symbol = (uint32_t)(d_samples_per_second / d_symbols_per_second);
    d_delay_after_sync   = d_samples_per_symbol / 4u;
    d_number_of_bins     = (uint32_t)(1u << d_sf);
    d_number_of_bins_hdr = (uint32_t)(1u << (d_sf-2));
    d_decim_factor       = d_samples_per_symbol / d_number_of_bins;
    d_energy_threshold   = 0.01f;
    d_whitening_sequence = gr::lora::prng_payload;
    d_fine_sync = 0;

    std::cout << "Bits (nominal) per symbol: \t"      << d_bits_per_symbol    << std::endl;
    std::cout << "Bins per symbol: \t"      << d_number_of_bins     << std::endl;
    std::cout << "Samples per symbol: \t"   << d_samples_per_symbol << std::endl;
    std::cout << "Decimation: \t\t"         << d_decim_factor       << std::endl;

    // Locally generated chirps
    build_ideal_chirps();

    d_state = DecoderState::DETECT;

    // FFT decoding preparations
    d_fft.resize(d_samples_per_symbol); //reduced the size of d_fft to the # of samples_per_symbol
    d_mult_hf.resize(d_samples_per_symbol);
    d_tmp.resize(d_number_of_bins);
    d_q  = fft_create_plan(d_samples_per_symbol, &d_mult_hf[0], &d_fft[0],     LIQUID_FFT_FORWARD, 0);
    d_qr = fft_create_plan(d_number_of_bins,     &d_tmp[0],     &d_mult_hf[0], LIQUID_FFT_BACKWARD, 0);
}

float decoder::stddev(const float *values, const uint32_t len, const float mean) {
    float variance = 0.0f;

    for (uint32_t i = 0u; i < len; i++) {
        const float temp = values[i] - mean;
        variance += temp * temp;
    }

    variance /= (float)len;
    //std::cout<< "Standard Dev for given array is calculated"<< std::endl; //comment by me
    return std::sqrt(variance);
    
    //standard deviation for the given array
}

void decoder::samples_to_file(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size) {
    //std::cout << "Given value array has been dumped into text file"<< std::endl;//comment by me
    //#ifdef DEBUG
        std::ofstream out_file;
        out_file.open(path.c_str(), std::ios::out | std::ios::binary);

        //for(std::vector<gr_complex>::const_iterator it = v.begin(); it != v.end(); ++it) {
        for (uint32_t i = 0u; i < length; i++) {
            out_file.write(reinterpret_cast<const char *>(&v[i]), elem_size);
        }

        out_file.close();
        
    //#else
        (void) path;
        (void) v;
        (void) length;
        (void) elem_size;
        //std::cout << "lenth of v is"<< std::type (v) << std::endl;//comment by me

    //#endif
}

inline void decoder::instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window) {
	if (window < 2u) {
	    std::cerr << "[LoRa Decoder] WARNING : window size < 2 !" << std::endl;
	    return;
	}

	/* instantaneous_phase */
	for (uint32_t i = 1u; i < window; i++) {
	    const float iphase_1 = std::arg(in_samples[i - 1]);
	          float iphase_2 = std::arg(in_samples[i]);

	    // Unwrapped loops from liquid_unwrap_phase
	    while ( (iphase_2 - iphase_1) >  M_PI ) iphase_2 -= 2.0f*M_PI;
	    while ( (iphase_2 - iphase_1) < -M_PI ) iphase_2 += 2.0f*M_PI;

	    out_ifreq[i - 1] = iphase_2 - iphase_1;
	}

	// Make sure there is no strong gradient if this value is accessed by mistake
	out_ifreq[window - 1] = out_ifreq[window - 2];
	//std::cout<< "Calculates instantaneous frequency for the given complex symbol"<< std::endl; //comment by me
	//std::cout<< "iphase_1: \t"<< iphase_1  <<std::endl; //comment by me
	//std::cout<< "iphase_2: \t"<<iphase_2 << std::endl; //comment by me
	//std::cout <<"iphase_1: \t\t"<< out_ifreq[window - 1] << std::endl;
}

inline gr_complex decoder::gr_expj(float phase)
{
    float t_imag, t_real;
    sincosf(phase, &t_imag, &t_real);
    return gr_complex(t_real, t_imag);
}

void decoder::build_ideal_chirps() 
{
    d_downchirp.resize(d_samples_per_symbol);
    d_upchirp.resize(d_samples_per_symbol);
    d_downchirp_ifreq.resize(d_samples_per_symbol);
    d_upchirp_ifreq.resize(d_samples_per_symbol);
    d_upchirp_ifreq_v.resize(d_samples_per_symbol*3);
    gr_complex tmp[d_samples_per_symbol*3];

    const double T       = -0.5 * d_bw * d_symbols_per_second;
    const double f0      = (d_bw / 2.0);
    const double pre_dir = 2.0 * M_PI;
    double t;
    gr_complex cmx       = gr_complex(1.0f, 1.0f);

    for (uint32_t i = 0u; i < d_samples_per_symbol; i++) {
        // Width in number of samples = samples_per_symbol
        // See https://en.wikipedia.org/wiki/Chirp#Linear
        t = d_dt * i;
        d_downchirp[i] = cmx * gr_expj(pre_dir * t * (f0 + T * t)); //ideal downshirp is created
        d_upchirp[i]   = cmx * gr_expj(pre_dir * t * (f0 + T * t) * -1.0f); //ideal upchirp is created
    }

    // Store instantaneous frequency
    instantaneous_frequency(&d_downchirp[0], &d_downchirp_ifreq[0], d_samples_per_symbol);
    instantaneous_frequency(&d_upchirp[0],   &d_upchirp_ifreq[0],   d_samples_per_symbol);

    samples_to_file("/tmp/downchirp", &d_downchirp[0], d_downchirp.size(), sizeof(gr_complex));
    samples_to_file("/tmp/upchirp",   &d_upchirp[0],   d_upchirp.size(),   sizeof(gr_complex));

    // Upchirp sequence
    memcpy(tmp, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
    memcpy(tmp+d_samples_per_symbol, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
    memcpy(tmp+d_samples_per_symbol*2, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
    instantaneous_frequency(tmp, &d_upchirp_ifreq_v[0], d_samples_per_symbol*3);
    //std::cout << "Locally generated chirps have been generated"<< std::endl; //comment by me
    // a way to get a visual of these chirps??
}

float decoder::cross_correlate_ifreq_fast(const float *samples_ifreq, const float *ideal_chirp, const uint32_t window) {
    float result = 0;
    volk_32f_x2_dot_prod_32f(&result, samples_ifreq, ideal_chirp, window);
    //std::cout<< "Random To-do_1"<< std::endl; //comment by me
    return result;
    
}

float decoder::cross_correlate_fast(const gr_complex *samples, const gr_complex *ideal_chirp, const uint32_t window) {
    gr_complex result = 0;
    volk_32fc_x2_conjugate_dot_prod_32fc(&result, samples, ideal_chirp, window);
    //std::cout<< "Random To-do_2"<< std::endl; //comment by me
    return abs(result);
    
}

float decoder::cross_correlate(const gr_complex *samples_1, const gr_complex *samples_2, const uint32_t window) {
    float result = 0.0f;

    for (uint32_t i = 0u; i < window; i++) {
        result += std::real(samples_1[i] * std::conj(samples_2[i]));
    }

    result /= (float)window;
    //std::cout<< "The correlation coefficient(when correlating the given complex symbols in the given window) has been returned"<< std::endl; //comment by me
    return result;
    
}

float decoder::cross_correlate_ifreq(const float *samples_ifreq, const std::vector<float>& ideal_chirp, const uint32_t to_idx) {
    float result = 0.0f;

    const float average   = std::accumulate(samples_ifreq  , samples_ifreq + to_idx, 0.0f) / (float)(to_idx);
    const float chirp_avg = std::accumulate(&ideal_chirp[0], &ideal_chirp[to_idx]  , 0.0f) / (float)(to_idx);
    const float sd        =   stddev(samples_ifreq   , to_idx, average)
                            * stddev(&ideal_chirp[0] , to_idx, chirp_avg);

    for (uint32_t i = 0u; i < to_idx; i++) {
        result += (samples_ifreq[i] - average) * (ideal_chirp[i] - chirp_avg) / sd;
    }

    result /= (float)(to_idx);
    //std::cout<< "The correlation coefficient of real signal is returned (in result)"<< std::endl; //comment by me
    return result;
    


}

void decoder::fine_sync(const gr_complex* in_samples, uint32_t bin_idx, int32_t search_space) {
    int32_t shift_ref = (bin_idx+1) * d_decim_factor;
    float samples_ifreq[d_samples_per_symbol];
    float max_correlation = 0.0f;
    int32_t lag = 0;

    instantaneous_frequency(in_samples, samples_ifreq, d_samples_per_symbol);

    for(int32_t i = -search_space+1; i < search_space; i++) {
        //float c = cross_correlate_fast(in_samples, &d_upchirp_v[shift_ref+i+d_samples_per_symbol], d_samples_per_symbol);
        float c = cross_correlate_ifreq_fast(samples_ifreq, &d_upchirp_ifreq_v[shift_ref+i+d_samples_per_symbol], d_samples_per_symbol);
        if(c > max_correlation) {
             max_correlation = c;
             lag = i;
         }
    }

    #ifdef DEBUG
        //d_debug << "FINE: " << -lag << std::endl;
    #endif

    d_fine_sync = -lag;
    //std::cout<< "we are in fine sync"<< std::endl; //comment by me

    //if(abs(d_fine_sync) >= d_decim_factor / 2)
    //    d_fine_sync = 0;
    //d_fine_sync = 0;
}

float decoder::detect_preamble_autocorr(const gr_complex *samples, const uint32_t window) {
    const gr_complex* chirp1 = samples;
    const gr_complex* chirp2 = samples + d_samples_per_symbol;
    float magsq_chirp1[window];
    float magsq_chirp2[window];
    float energy_chirp1 = 0;
    float energy_chirp2 = 0;
    float autocorr = 0;
    //uint32_t array_test[];
    gr_complex dot_product;

    volk_32fc_x2_conjugate_dot_prod_32fc(&dot_product, chirp1, chirp2, window);
    volk_32fc_magnitude_squared_32f(magsq_chirp1, chirp1, window);
    volk_32fc_magnitude_squared_32f(magsq_chirp2, chirp2, window);
    volk_32f_accumulator_s32f(&energy_chirp1, magsq_chirp1, window);
    volk_32f_accumulator_s32f(&energy_chirp2, magsq_chirp2, window);

    autocorr = abs(dot_product / gr_complex(sqrt(energy_chirp1 * energy_chirp2), 0));
    if (autocorr >= 0.9f){
        //array_test = window;
        //std::cout <<  "array test" << window << std::endl;
       // std::cout <<  "array test length" << sizeof(window) << std::endl;
    }
    //std::cout<< "The auto correlation appraoch to approximate the preamble"<< std::endl; //comment by me ??

    return autocorr;
    
}

float decoder::detect_downchirp(const gr_complex *samples, const uint32_t window) {
    float samples_ifreq[window];
    instantaneous_frequency(samples, samples_ifreq, window);
    //std::cout<< "Base method to start downchirp correlation and return the correlation coefficient"<< std::endl; //comment by me

    return cross_correlate_ifreq(samples_ifreq, d_downchirp_ifreq, window - 1u);
    
}

float decoder::detect_upchirp(const gr_complex *samples, const uint32_t window, int32_t *index) {
    float samples_ifreq[window*2];
    instantaneous_frequency(samples, samples_ifreq, window*2);
    //std::cout<< " Base method to start upchirp detection `sliding_norm_cross_correlate_upchirp` is called"<< std::endl; //comment by me

    return sliding_norm_cross_correlate_upchirp(samples_ifreq, window, index);
    
}

float decoder::sliding_norm_cross_correlate_upchirp(const float *samples_ifreq, const uint32_t window, int32_t *index) {
    float max_correlation = 0;
    #ifdef DEBUG
        std::cout<< " Cross correlation for upchirp "<< std::endl; //comment by me
    #endif
    // Cross correlate
    for (uint32_t i = 0; i < window; i++) {
     //const float samp;
     //samp = samples_ifreq + i;
     const float max_corr = cross_correlate_ifreq_fast(samples_ifreq + i, &d_upchirp_ifreq[0], window - 1u);

     if (max_corr > max_correlation) {
         *index = i;
         max_correlation = max_corr;               
     }
     //samples_to_file("/home/dell/Download/lora_test.txt", &samples_ifreq + i[0], window , sizeof(samples_ifreq + i));
     //std::cout<< "upchirp values"<< samples_ifreq + i <<std::endl; 
     //std::cout<< "size"<< sizeof(samples_ifreq + i) <<std::endl;
      //else cout << "Unable to open file";
    //std::cout<< "upchirp array type"<< type (samples_ifreq + i) <<std::endl;
    }


    // std::cout<< "Shift of the symbol has been corrected to match the ideal upchirp (using sliding correlation"<< std::endl; //comment by me
 return max_correlation;
     
}

uint32_t decoder::max_frequency_gradient_idx(const gr_complex *samples) {
    float samples_ifreq[d_samples_per_symbol];
    float samples_ifreq_avg[d_number_of_bins];

    samples_to_file("/tmp/data", &samples[0], d_samples_per_symbol, sizeof(gr_complex));

    instantaneous_frequency(samples, samples_ifreq, d_samples_per_symbol);

    for(uint32_t i = 0; i < d_number_of_bins; i++) {
        volk_32f_accumulator_s32f(&samples_ifreq_avg[i], &samples_ifreq[i*d_decim_factor], d_decim_factor);
        samples_ifreq_avg[i] /= d_decim_factor;
    }

    float max_gradient = 0.1f;
    float gradient = 0.0f;
    uint32_t max_index = 0;
    for (uint32_t i = 1u; i < d_number_of_bins; i++) {
        gradient = samples_ifreq_avg[i - 1] - samples_ifreq_avg[i];
        if (gradient > max_gradient) {
            max_gradient = gradient;
            max_index = i+1;
        }
    }
    //std::cout<< "the index of the bin containing the frequency change is returned"<< std::endl; //comment by me
    return (d_number_of_bins - max_index) % d_number_of_bins;
   
}


/**
 *  Correct the interleaving by extracting each column of bits after rotating to the left.
 *  <br/>(The words were interleaved diagonally, by rotating we make them straight into columns.)
 */
void decoder::deinterleave(const uint32_t ppm) {
    const uint32_t bits_per_word = d_words.size();
    const uint32_t offset_start  = ppm - 1u;

    std::vector<uint8_t> words_deinterleaved(ppm, 0u);

    if (bits_per_word > 8u) {
        // Not sure if this can ever occur. It would imply coding rate high than 4/8 e.g. 4/9.
        std::cerr << "[LoRa Decoder] WARNING : Deinterleaver: More than 8 bits per word. uint8_t will not be sufficient!\nBytes need to be stored in intermediate array and then packed into words_deinterleaved!" << std::endl;
        exit(1);
    }

    for (uint32_t i = 0u; i < bits_per_word; i++) {
        const uint32_t word = gr::lora::rotl(d_words[i], i, ppm);

        for (uint32_t j = (1u << offset_start), x = offset_start; j; j >>= 1u, x--) {
            words_deinterleaved[x] |= !!(word & j) << i;
        }
    }

    #ifdef DEBUG
        print_vector_bin(d_debug, words_deinterleaved, "D", sizeof(uint8_t) * 8u);
        //print_interleave_matrix(d_debug, d_words, ppm);
    #endif

    // Add to demodulated data
    d_demodulated.insert(d_demodulated.end(), words_deinterleaved.begin(), words_deinterleaved.end());

    // Cleanup
    d_words.clear();
    //std::cout<< "raw modulated words are dinterleaved. No output"<< std::endl; //comment by me
}


void decoder::decode(const bool is_header) {
    static const uint8_t shuffle_pattern[] = {5, 0, 1, 2, 4, 3, 6, 7};

    // For determining shuffle pattern
    //if (!is_header)
    //    values_to_file("/tmp/before_deshuffle", &d_demodulated[0], d_demodulated.size(), 8);

    deshuffle(shuffle_pattern, is_header);

    // For determining whitening sequence
    //if (!is_header)
    //    values_to_file("/tmp/after_deshuffle", &d_words_deshuffled[0], d_words_deshuffled.size(), 8);

    dewhiten(is_header ? gr::lora::prng_header : d_whitening_sequence);

    //if (!is_header)
    //    values_to_file("/tmp/after_dewhiten", &d_words_dewhitened[0], d_words_dewhitened.size(), 8);

    hamming_decode(is_header);
    //std::cout<< "Actual Payload extracted"<< std::endl; //comment by me
    //std::cout<< "THE ENd"<< std::endl; //comment by me
}



bool decoder::demodulate(const gr_complex *samples, const bool is_header) {
    // DBGR_TIME_MEASUREMENT_TO_FILE("SFxx_method");

    // DBGR_START_TIME_MEASUREMENT(false, "only");
    uint32_t bin_idx = max_frequency_gradient_idx(samples);
    //uint32_t bin_idx = get_shift_fft(samples);
    fine_sync(samples, bin_idx, std::max(d_decim_factor / 4u, 2u));

    // DBGR_INTERMEDIATE_TIME_MEASUREMENT();

    // Header has additional redundancy
    if (is_header || d_sf > 10) {
        bin_idx = std::lround(bin_idx / 4.0f) % d_number_of_bins_hdr;
    }

    // Decode (actually gray encode) the bin to get the symbol value
    const uint32_t word = bin_idx ^ (bin_idx >> 1u);

    #ifdef DEBUG
        d_debug << gr::lora::to_bin(word, is_header ? d_sf - 2u : d_sf) << " " << bin_idx  << std::endl;
    #endif
    d_words.push_back(word);

    // Look for 4+cr symbols and stop
    if (d_words.size() == (4u + d_phdr.cr)) {
        // Deinterleave
        deinterleave((is_header || d_sf > 10) ? d_sf - 2u : d_sf);

        //std::cout<< "A block is ready for decoding"<< std::endl; //comment by me
        return true; // Signal that a block is ready for decoding
        //std::cout<< "A block is ready for decoding"<< std::endl; //comment by me
        //std::cout<< "raw modulated words are dinterleaved. No output"<< std::endl; //comment by me
    }
    return false; // We need more words in order to decode a block
    
}


void decoder::deshuffle(const uint8_t *shuffle_pattern, const bool is_header) {
    const uint32_t to_decode = is_header ? 5u : d_demodulated.size();
    const uint32_t len       = sizeof(shuffle_pattern) / sizeof(uint8_t);
    uint8_t result;

    for (uint32_t i = 0u; i < to_decode; i++) {
        result = 0u;

        for (uint32_t j = 0u; j < len; j++) {
            result |= !!(d_demodulated[i] & (1u << shuffle_pattern[j])) << j;
        }

        d_words_deshuffled.push_back(result);
    }

    #ifdef DEBUG
        print_vector_bin(d_debug, d_words_deshuffled, "S", sizeof(uint8_t)*8);
        d_debug << std::endl;
    #endif

    // We're done with these words
    if (is_header){
        d_demodulated.erase(d_demodulated.begin(), d_demodulated.begin() + 5u);
    } else {
        d_demodulated.clear();
    }
    //std::cout<< "1/3: Demodulated words have been deshuffled"<< std::endl; //comment by me
}


void decoder::dewhiten(const uint8_t *prng) {
    const uint32_t len = d_words_deshuffled.size();

    for (uint32_t i = 0u; i < len; i++) {
        uint8_t xor_b = d_words_deshuffled[i] ^ prng[i];
        d_words_dewhitened.push_back(xor_b);
    }

    #ifdef DEBUG
        print_vector_bin(d_debug, d_words_dewhitened, "W", sizeof(uint8_t) * 8);
    #endif

    d_words_deshuffled.clear();
    //std::cout<< "2/3: deshuffled words are dewhitened using XOR"<< std::endl; //comment by me
}

void decoder::hamming_decode(bool is_header) {
    switch(d_phdr.cr) {
        case 4: case 3: // Hamming(8,4) or Hamming(7,4)
            hamming_decode_soft(is_header);
            break;
        case 2: case 1: // Hamming(6,4) or Hamming(5,4)
            // TODO: Report parity error to the user
            extract_data_only(is_header);
            break;
    }

    d_words_dewhitened.clear();
    //std::cout<< "3/3:Data in relevant bytes extracted"<< std::endl; //comment by me

    /*
    fec_scheme fs = LIQUID_FEC_HAMMING84;
    unsigned int n = ceil(d_words_dewhitened.size() * 4.0f / (4.0f + d_phdr.cr));

    unsigned int k = fec_get_enc_msg_length(fs, n);
    fec hamming = fec_create(fs, NULL);

    fec_decode(hamming, n, &d_words_dewhitened[0], out_data);

    d_words_dewhitened.clear();
    fec_destroy(hamming);*/
}


void decoder::hamming_decode_soft(bool is_header) {
    uint32_t len = d_words_dewhitened.size();
    for (uint32_t i = 0u; i < len; i += 2u) {
        const uint8_t d2 = (i + 1u < len) ? gr::lora::hamming_decode_soft_byte(d_words_dewhitened[i + 1u]) : 0u;
        const uint8_t d1 = gr::lora::hamming_decode_soft_byte(d_words_dewhitened[i]);

        if(is_header)
            d_decoded.push_back((d1 << 4u) | d2);
        else
            d_decoded.push_back((d2 << 4u) | d1);
    }
    //std::cout<< "Hamming(8,4) decoding performed on each byte"<< std::endl; //comment by me
}

void decoder::extract_data_only(bool is_header) {
    static const uint8_t data_indices[4] = {1, 2, 3, 5};
    uint32_t len = d_words_dewhitened.size();

    for (uint32_t i = 0u; i < len; i += 2u) {
        const uint8_t d2 = (i + 1u < len) ? gr::lora::select_bits(d_words_dewhitened[i + 1u], data_indices, 4u) & 0xFF : 0u;
        const uint8_t d1 = (gr::lora::select_bits(d_words_dewhitened[i], data_indices, 4u) & 0xFF);

        if(is_header)
            d_decoded.push_back((d1 << 4u) | d2);
        else
            d_decoded.push_back((d2 << 4u) | d1);
    }
    //std::cout<< "Data extracted in given bytes"<< std::endl; //comment by me
}

void decoder::set_abs_threshold(const float threshold) {
    d_energy_threshold = gr::lora::clamp(threshold, 0.0f, 20.0f);
    //std::cout<< "Absolute Threshold set"<< std::endl; //comment by me
}

void decoder::print_complex_pointers(const gr_complex* input, int size)
{
    for (int i=0; i<size; i++)
        std::cout << "input["<<i<<"] = " << input[i] << std::endl;
}

uint32_t decoder::get_d_samples_per_symbol()
{
    return d_samples_per_symbol;
}

int32_t decoder::get_d_fine_sync()
{
    return d_fine_sync;
}

uint32_t decoder::get_d_delay_after_sync()
{
    return d_delay_after_sync;
}

uint32_t decoder::get_d_decim_factor()
{
    return d_decim_factor;
}

uint32_t decoder::get_d_number_of_bins()
{
    return d_number_of_bins;
}

std::vector<uint8_t> decoder::get_d_decoded()
{
    return d_decoded;
}

void decoder::set_d_phdr(std::vector<uint8_t> d_decoded)
{
    d_phdr.length = d_decoded[0];
    d_phdr.crc_msn = d_decoded[1];
    d_phdr.has_mac_crc = d_decoded[2];
    d_phdr.cr = d_decoded[3];
    d_phdr.crc_lsn = d_decoded[4];
    d_phdr.reserved = d_decoded[5];
}

void decoder::clear_d_decoded()
{
    d_decoded.clear();
}

void decoder::set_d_payload_length(uint32_t value)
{
    d_payload_length = value;
}

loraphy_header_t decoder::get_d_phdr()
{
    return d_phdr;
}


std::string decoder::logic_DETECT(const gr_complex * input)
{
    float correlation = detect_preamble_autocorr(input, d_samples_per_symbol);
    //std::cout<< "trying to detect preamble "<< std::endl; //comment by me

    std::cout << correlation << std::endl;

    if (correlation >= 0.90f) {
        std::cout << correlation << std::endl;
        std::cout<< "corr >= 0.90f preamble detected. Lora Frame starts .. "<< std::endl; //comment by me
        #ifdef DEBUG
            d_debug << "Ca: " << correlation << std::endl;
        #endif
        d_corr_fails = 0u;
        std::cout<< "going to decoder state SYNC "<< std::endl; //comment by me
        
        // *d_state = DecoderState::SYNC;
        return "SYNC";
    }
    else {
        // *d_state = DecoderState::DETECT;
        return "DETECT";
    }
}

std::string decoder::logic_SYNC(const gr_complex * input, int * i)
{
    std::cout<< "going to detect upchirp "<< std::endl; //comment by me
    detect_upchirp(input, d_samples_per_symbol, i);
    //float cfo = experimental_determine_cfo(&input[i], d_samples_per_symbol);
    //pmt::pmt_t kv = pmt::cons(pmt::intern(std::string("cfo")), pmt::from_double(cfo));
    //message_port_pub(pmt::mp("control"), kv);

    std::cout<< "going to decoder state FIND SFD "<< std::endl; //comment by me
    return "FIND_SFD";
}

std::string decoder::logic_FIND_SFD(const gr_complex * input)
{
    std::cout << "going to detect downchirp " << std::endl; //comment by me
    const float c = detect_downchirp(input, d_samples_per_symbol);
    
    #ifdef DEBUG
        d_debug << "Cd: " << c << std::endl;
    #endif

    if (c > 0.96f) {
        #ifdef DEBUG
            std::cout << " Sync success. Correlation result > 0.96f " << std::endl; //comment by me
            d_debug << "SYNC: " << c << std::endl;
            samples_to_file("/tmp/sync", input, d_samples_per_symbol, sizeof(gr_complex));
            std::cout <<"going to decode HEADER " << std::endl; //comment by me
        #endif
        
        return "PAUSE";
    } else {
        if(c < -0.97f) {
            std::cout << "correlation result < -0.97f" << std::endl; //comment by me
            std::cout << "going to fine_sync" <<std::endl; //comment by me
            fine_sync(input, d_number_of_bins-1, d_decim_factor * 4);
        } else {
            d_corr_fails++;
            std::cout << "correlation failed " <<std::endl; //comment by me
        }

        if (d_corr_fails > 4u) {
            std::cout << "Lost Sync. Going to detect next preamble " << std::endl; //comment by me
            #ifdef DEBUG
                d_debug << "Lost sync" << std::endl;
            #endif
            
            return "DETECT";
        } else{
            return "FIND_SFD";
        }
    }
}

std::string decoder::logic_PAUSE()
{
    //std::cout << "pausing before going to decode HEADER " << std::endl; //comment by me
    return "DECODE_HEADER";
    //std::cout<< "case:pause end"<< std::endl; //comment by me
}


std::string decoder::logic_DECODE_HEADER(const gr_complex * input)
{
    d_phdr.cr = 4u;
    std::cout << "Trying to decode header " << std::endl; //comment by me
    if (demodulate(input, true)) {
        decode(true);
        std::cout << "Demod and decode of header success. Printing header " << std::endl; //comment by me
        gr::lora::print_vector_hex(std::cout, &d_decoded[0], d_decoded.size(), false);
        uint32_t len_tap = sizeof(loratap_header_t); //by me 
                                uint32_t len_phy = sizeof(loraphy_header_t); //by me
                                //samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/interdigital2b", &input[i-500], 1000, sizeof(gr_complex));//by me
                                samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/header/data-09-2", &input[len_tap], len_phy, sizeof(gr_complex));//by me
        std::cout << "Tap Header Length ="<< len_tap << std::endl;
        std::cout << "Phy Header Length ="<< len_phy << std::endl;
        std::cout << "Tap Header Length ="<< len_tap << std::endl;
        //step1 = &d_decoded; by me
        //int step2 = d_decoded; by me
        //std::cout << "d_decoded" << d_decoded[0] << std::endl; //by me 
        //std::cout << "step2" << step2 << std::endl; //by me
        memcpy(&d_phdr, &d_decoded[0], sizeof(loraphy_header_t));
        std::cout << "d_phdr" << &d_phdr << std::endl; //by me 
        d_decoded.clear();

        d_payload_length = d_phdr.length + MAC_CRC_SIZE * d_phdr.has_mac_crc;
        //std::cout << "corresponding payload length " << d_payload_length << std::endl; //comment by me
        //d_phy_crc = SM(decoded[1], 4, 0xf0) | MS(decoded[2], 0xf0, 4);
        uint32_t len1 = sizeof(loratap_header_t) + sizeof(loraphy_header_t); //+ d_payload_length;
        std::cout << "Length of Lora packet: "<< len1 << std::endl; //comment by me

        // Calculate number of payload symbols needed
        uint8_t redundancy = (d_sf > 10 ? 2 : 0);
        const int symbols_per_block = d_phdr.cr + 4u;
        const float bits_needed     = float(d_payload_length) * 8.0f;
        const float symbols_needed  = bits_needed * (symbols_per_block / 4.0f) / float(d_sf - redundancy);
        const int blocks_needed     = (int)std::ceil(symbols_needed / symbols_per_block);
        d_payload_symbols     = blocks_needed * symbols_per_block;
        int size_payload = sizeof(d_payload_length);
        //std::cout << "Size in Bytes of payload = " << size_payload << "bytes" << std::endl; //comment by me
        std::cout << "  Payload length and payload symbols " << std::endl; //comment by me
        std::cout << "LEN: " << d_payload_length << " (" << d_payload_symbols << " symbols)" << std::endl; //comment by me
        #ifdef DEBUG
            d_debug << "LEN: " << d_payload_length << " (" << d_payload_symbols << " symbols)" << std::endl;
        #endif
        std::cout << "going to decode PAYLOAD " << std::endl; //comment by me
        
        return "DECODE_PAYLOAD";
    } 
    else { 
        std::cout << "Could not decode header " << std::endl; //comment by me
        return "DECODE_HEADER";
    }

    //std::cout<< "case:decode_header end"<< std::endl; //comment by me
}


std::string decoder::logic_DECODE_PAYLOAD(const gr_complex * input)
{
    //**************************************************************************
    // Failsafe if decoding length reaches end of actual data == noise reached?
    // Could be replaced be rejecting packets with CRC mismatch...
    if (std::abs(input[0]) < d_energy_threshold) {
        std::cout << " input has energy lower than threshold. End of data. payload length in HDR is wrong " << std::endl; //comment by me 
        //printf("\n*** Decode payload reached end of data! (payload length in HDR is wrong)\n");
        d_payload_symbols = 0;
    }
    //**************************************************************************

    if (demodulate(input, false)) {
        d_payload_symbols -= (4u + d_phdr.cr);

        if (d_payload_symbols <= 0) {
            decode(false);
            //std::cout << "demod and decode of payload success. Printing .. " << std::endl; //comment by me
            d_payload_length = d_phdr.length + MAC_CRC_SIZE * d_phdr.has_mac_crc;//by me
            uint32_t len1 = sizeof(loratap_header_t) + sizeof(loraphy_header_t) + d_payload_length;//by me
            uint32_t duration = len1 - 500; //by me
            samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/eop/data-09-2", &input[len1-21], 1000, sizeof(gr_complex));//by me
            gr::lora::print_vector_hex(std::cout, &d_decoded[0], d_payload_length, true);
            gr::lora::dev_address(std::cout, &d_decoded[0],true);//DevAddress
            // msg_lora_frame();  // msg is for gnuradio
            std::cout << "d_payload_length:" << d_payload_length<< std::endl; //comment by me
            //gr::lora::print_vector_hex(std::cout << "Array of words deshuffled" << d_words_deshuffled[1]);//<< std::endl;
            //std::cout<< "DevAddress above"<< std::endl;
            std::cout << "End of Lora Packet. going back to detecting preamble " << std::endl; //comment by me
            //samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/interdigital3d", &input[0], 1000, sizeof(gr_complex));//by me doesn't quit terminal
            d_decoded.clear();

            return "DETECT";
        }
    }

    return "DECODE_PAYLOAD";
    //std::cout<< "case:decode_payload end"<< std::endl; //comment by me
}



int decoder::work(const gr_complex * input)
{
    d_fine_sync = 0; // Always reset fine sync

    switch (d_state) {
        case DecoderState::DETECT: {
            float correlation = detect_preamble_autocorr(input, d_samples_per_symbol);
            #ifdef DEBUG
                //std::cout<< "trying to detect preamble "<< std::endl; //comment by me
                std::cout << correlation << std::endl;
            #endif

            if (correlation >= 0.90f) {
                #ifdef DEBUG
                    std::cout << correlation << std::endl;
                    std::cout<< "corr >= 0.90f preamble detected. Lora Frame starts .. "<< std::endl; //comment by me
                    d_debug << "Ca: " << correlation << std::endl;
                    std::cout<< "going to decoder state SYNC "<< std::endl; //comment by me
                #endif
                d_corr_fails = 0u;
                d_state = DecoderState::SYNC;  
                break;
            }

            return d_samples_per_symbol;
            // std::cout<< "Preamble Detected"<< std::endl; //comment by me

        }

        case DecoderState::SYNC: {
            int i = 0;
            #ifdef DEBUG
                std::cout<< "going to detect upchirp "<< std::endl; //comment by me
            #endif
            detect_upchirp(input, d_samples_per_symbol, &i);

            #ifdef DEBUG
                std::cout<< "going to decoder state FIND SFD "<< std::endl; //comment by me
            #endif
            d_state = DecoderState::FIND_SFD;

            return i;
        }

        case DecoderState::FIND_SFD: {
            #ifdef DEBUG
                std::cout << "going to detect downchirp " << std::endl; //comment by me
            #endif
            
            const float c = detect_downchirp(input, d_samples_per_symbol);
            
            #ifdef DEBUG
                d_debug << "Cd: " << c << std::endl;
            #endif

            if (c > 0.96f) {
                #ifdef DEBUG
                    std::cout << " Sync success. Correlation result > 0.96f " << std::endl; //comment by me
                    d_debug << "SYNC: " << c << std::endl;
                    std::cout <<"going to decode HEADER " << std::endl; //comment by me
                #endif
                d_state = DecoderState::PAUSE;
            } else {
                if(c < -0.97f) {
                    #ifdef DEBUG
                        std::cout << "correlation result < -0.97f" << std::endl; //comment by me
                        std::cout << "going to fine_sync" <<std::endl; //comment by me
                    #endif
                    fine_sync(input, d_number_of_bins-1, d_decim_factor * 4);
                } else {
                    d_corr_fails++;
                    #ifdef DEBUG
                        std::cout << "correlation failed " <<std::endl; //comment by me
                    #endif
                }

                if (d_corr_fails > 4u) {
                    d_state = DecoderState::DETECT;
                    #ifdef DEBUG
                        d_debug << "Lost sync" << std::endl;
                    #endif
                }
            }

            return d_samples_per_symbol + d_fine_sync;
        }

        case DecoderState::PAUSE: {
            d_state = DecoderState::DECODE_HEADER;
            return (int32_t) d_samples_per_symbol + d_delay_after_sync;
        }

        case DecoderState::DECODE_HEADER: {
            d_phdr.cr = 4u;
            #ifdef DEBUG
                std::cout << "Trying to decode header " << std::endl; //comment by me
            #endif
            if (demodulate(input, true)) {
                decode(true);

                std::cout << "Demod and decode of header success. Printing header " << std::endl; //comment by me
                gr::lora::print_vector_hex(std::cout, &d_decoded[0], d_decoded.size(), true);
                
                uint32_t len_tap = sizeof(loratap_header_t); //by me 
                                        uint32_t len_phy = sizeof(loraphy_header_t); //by me
                                        //samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/interdigital2b", &input[i-500], 1000, sizeof(gr_complex));//by me
                                        samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/header/data-09-2", &input[len_tap], len_phy, sizeof(gr_complex));//by me
                #ifdef DEBUG
                    std::cout << "Tap Header Length ="<< len_tap << std::endl;
                    std::cout << "Phy Header Length ="<< len_phy << std::endl;
                    std::cout << "Tap Header Length ="<< len_tap << std::endl;
                #endif
                //step1 = &d_decoded; by me
                //int step2 = d_decoded; by me
                //std::cout << "d_decoded" << d_decoded[0] << std::endl; //by me 
                //std::cout << "step2" << step2 << std::endl; //by me
                memcpy(&d_phdr, &d_decoded[0], sizeof(loraphy_header_t));
                #ifdef DEBUG
                    std::cout << "d_phdr" << &d_phdr << std::endl; //by me 
                #endif
                d_decoded.clear();

                d_payload_length = d_phdr.length + MAC_CRC_SIZE * d_phdr.has_mac_crc;
                //std::cout << "corresponding payload length " << d_payload_length << std::endl; //comment by me
                //d_phy_crc = SM(decoded[1], 4, 0xf0) | MS(decoded[2], 0xf0, 4);
                uint32_t len1 = sizeof(loratap_header_t) + sizeof(loraphy_header_t); //+ d_payload_length;
                #ifdef DEBUG
                    std::cout << "Length of Lora packet: "<< len1 << std::endl; //comment by me
                #endif

                // Calculate number of payload symbols needed
                uint8_t redundancy = (d_sf > 10 ? 2 : 0);
                const int symbols_per_block = d_phdr.cr + 4u;
                const float bits_needed     = float(d_payload_length) * 8.0f;
                const float symbols_needed  = bits_needed * (symbols_per_block / 4.0f) / float(d_sf - redundancy);
                const int blocks_needed     = (int)std::ceil(symbols_needed / symbols_per_block);
                d_payload_symbols     = blocks_needed * symbols_per_block;
                int size_payload = sizeof(d_payload_length);
                #ifdef DEBUG
                    std::cout << "Size in Bytes of payload = " << size_payload << "bytes" << std::endl; //comment by me
                    std::cout << "  Payload length and payload symbols " << std::endl; //comment by me
                    d_debug << "LEN: " << d_payload_length << " (" << d_payload_symbols << " symbols)" << std::endl;
                    std::cout << "going to decode PAYLOAD " << std::endl; //comment by me
                #endif
                d_state = DecoderState::DECODE_PAYLOAD;
            } 
            else 
            { 
                #ifdef DEBUG
                    std::cout << "Could not decode header " << std::endl; //comment by me
                #endif
            }

            return (int32_t) d_samples_per_symbol + d_fine_sync;
        }

        case DecoderState::DECODE_PAYLOAD: {
            //**************************************************************************
            // Failsafe if decoding length reaches end of actual data == noise reached?
            // Could be replaced be rejecting packets with CRC mismatch...
            if (std::abs(input[0]) < d_energy_threshold) {
                #ifdef DEBUG
                    std::cout << " input has energy lower than threshold. End of data. payload length in HDR is wrong " << std::endl; //comment by me 
                #endif
                //printf("\n*** Decode payload reached end of data! (payload length in HDR is wrong)\n");
                d_payload_symbols = 0;
            }
            //**************************************************************************

            if (demodulate(input, false)) {
                d_payload_symbols -= (4u + d_phdr.cr);

                if (d_payload_symbols <= 0) {
                    decode(false);
                    //std::cout << "demod and decode of payload success. Printing .. " << std::endl; //comment by me
                    d_payload_length = d_phdr.length + MAC_CRC_SIZE * d_phdr.has_mac_crc;//by me
                    uint32_t len1 = sizeof(loratap_header_t) + sizeof(loraphy_header_t) + d_payload_length;//by me
                    uint32_t duration = len1 - 500; //by me
                    // samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/eop/data-09-2", &input[len1-21], 1000, sizeof(gr_complex));//by me

                    std::cout << "Demod and decode of Payload success. Printing payload " << std::endl; //comment by me
                    gr::lora::print_vector_hex(std::cout, &d_decoded[0], d_payload_length, true);
                    
                    #ifdef DEBUG    
                        // gr::lora::dev_address(std::cout, &d_decoded[0],true);//DevAddress
                        // msg_lora_frame();
                        std::cout << "d_payload_length:" << d_payload_length<< std::endl; //comment by me
                        //gr::lora::print_vector_hex(std::cout << "Array of words deshuffled" << d_words_deshuffled[1]);//<< std::endl;
                        //std::cout<< "DevAddress above"<< std::endl;
                        std::cout << "End of Lora Packet. going back to detecting preamble " << std::endl; //comment by me
                        //samples_to_file("/home/dell/Downloads/gr-lora-0.6.2/lib/interdigital3d", &input[0], 1000, sizeof(gr_complex));//by me doesn't quit terminal
                    #endif

                    d_state = DecoderState::DETECT;
                    d_decoded.clear();
                }
            }

            return (int32_t) d_samples_per_symbol+d_fine_sync;
        }

        case DecoderState::STOP: {
            return d_samples_per_symbol;
        }

        default: {
            std::cerr << "[LoRa Decoder] WARNING : No state! Shouldn't happen\n";
            break;
        }
    }

    return 0;
}