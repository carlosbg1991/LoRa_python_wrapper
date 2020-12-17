/* -*- c++ -*- */
/*
 * Copyright 2017 Pieter Robyns, William Thenaers.
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


#ifndef DECODER_H_
#define DECODER_H_

#include <string>
#include <vector>
#include <fstream>
#include <stdint.h>
#include <liquid/liquid.h>
#include <volk/volk.h>
#include "loraphy.h"


typedef std::complex<float>  gr_complex;

class decoder
{

	private:
		/**
		*  \brief  Return the DecoderState as string for debugging purposes.
		*
		*  \param  s
		*          The state to return to string.
		*/
		/*static std::string state_to_string(DecoderState s) {
		static std::string DecoderStateLUT[] = { "DETECT", "SYNC", "FIND_SFD", "PAUSE", "DECODE_HEADER", "DECODE_PAYLOAD", "STOP" };
		return DecoderStateLUT[ (size_t)s ];
		}*/

		/**
		*  \brief  **LoRa Decoder**
		*          <br/>The main class for the LoRa decoder.
		*          Contains all variables and methods necessary for succesfully decoding LoRa PHY.
		*          <br/>Only the sample rate and spreading factor are needed.
		*          The other settings, like packet length and coding rate, are extracted from the (explicit) HDR.
		*/
		std::vector<gr_complex> d_downchirp;        ///< The complex ideal downchirp.
		std::vector<float>      d_downchirp_ifreq;  ///< The instantaneous frequency of the ideal downchirp.

		std::vector<gr_complex> d_upchirp;          ///< The complex ideal upchirp.
		std::vector<float>      d_upchirp_ifreq;    ///< The instantaneous frequency of the ideal upchirp.
		std::vector<float>      d_upchirp_ifreq_v;  ///< The instantaneous frequency of the ideal upchirp.

		std::vector<liquid_float_complex> d_fft;              ///< Vector containing the FFT resuls.
		std::vector<liquid_float_complex> d_mult_hf;          ///< Vector containing the FFT decimation.
		std::vector<liquid_float_complex> d_tmp;              ///< Vector containing the FFT decimation.

		// std::vector<gr_complex> d_fft;              ///< Vector containing the FFT resuls.
		// std::vector<gr_complex> d_mult_hf;          ///< Vector containing the FFT decimation.
		// std::vector<gr_complex> d_tmp;              ///< Vector containing the FFT decimation.

		uint8_t          d_sf;                      ///< The Spreading Factor.
		uint32_t         d_bw;                      ///< The receiver bandwidth (fixed to `125kHz`).
		loraphy_header_t d_phdr;                    ///< LoRa PHY header.
		uint16_t         d_mac_crc;                 ///< The MAC CRC.
		double           d_bits_per_second;         ///< Indicator of how many bits are transferred each second.
		uint32_t         d_delay_after_sync;        ///< The number of samples to skip in `DecoderState::PAUSE`.
		uint32_t         d_samples_per_second;      ///< The number of samples taken per second by GNU Radio.
		double           d_symbols_per_second;      ///< Indicator of how many symbols (read: chirps) are transferred each second.
		double           d_bits_per_symbol;         ///< The number of bits each of the symbols contain.
		uint32_t         d_samples_per_symbol;      ///< The number of samples in one symbol.
		double           d_period;                  ///< Period of the symbol.
		uint32_t         d_number_of_bins;          ///< Indicates in how many parts or bins a symbol is decimated, i.e. the max value to decode out of one payload symbol.
		uint32_t         d_number_of_bins_hdr;      ///< Indicates in how many parts or bins a HDR symbol is decimated, i.e. the max value to decode out of one HDR symbol.
		int32_t          d_payload_symbols;         ///< The number of symbols needed to decode the payload. Calculated from an indicator in the HDR.
		uint32_t         d_payload_length;          ///< The number of words after decoding the HDR or payload. Calculated from an indicator in the HDR.
		uint32_t         d_corr_fails;              ///< Indicates how many times the correlation failed. After some tries, the state will revert to `DecoderState::DETECT`.
		float            d_energy_threshold;        ///< The absolute threshold to distinguish signal from noise.
		const uint8_t*   d_whitening_sequence;      ///< A pointer to the whitening sequence to be used in decoding. Determined by the SF in the ctor.

		std::vector<uint32_t> d_words;              ///< Vector containing the demodulated words.
		std::vector<uint8_t>  d_demodulated;        ///< Vector containing the words after deinterleaving.
		std::vector<uint8_t>  d_words_deshuffled;   ///< Vector containing the words after deshuffling.
		std::vector<uint8_t>  d_words_dewhitened;   ///< Vector containing the words after dewhitening.
		std::vector<uint8_t>  d_decoded;            ///< Vector containing the words after Hamming decode or the final decoded words.

		fftplan d_q;                                ///< The LiquidDSP::FFT_Plan.
		fftplan d_qr;                               ///< The LiquidDSP::FFT_Plan in reverse.

		uint32_t      d_decim_factor;               ///< The number of samples (data points) in each bin.
		float         d_cfo_estimation;             ///< An estimation for the current Center Frequency Offset.
		double        d_dt;                         ///< Indicates how fast the frequency changes in a symbol (chirp).
		int32_t 	  d_fine_sync;

	public:
		/**
		*  \brief  **DecoderState** : Each state the LoRa decoder can be in.
		*/
		enum class DecoderState {
		    DETECT,
		    SYNC,
		    FIND_SFD,
		    PAUSE,
		    DECODE_HEADER,
		    DECODE_PAYLOAD,
		    STOP
		};
		DecoderState            d_state;            ///< Holds the current state of the decoder (state machine).

		// extern void build_ideal_chirps(void);
		decoder(float samp_rate, int sf);

		inline gr_complex gr_expj(float phase);
		float stddev(const float *values, const uint32_t len, const float mean);

		void samples_to_file(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size);
		inline void instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window);
		void build_ideal_chirps();

		float cross_correlate_ifreq_fast(const float *samples_ifreq, const float *ideal_chirp, const uint32_t window);
		float cross_correlate_fast(const gr_complex *samples, const gr_complex *ideal_chirp, const uint32_t window);
		float cross_correlate(const gr_complex *samples_1, const gr_complex *samples_2, const uint32_t window);
		float cross_correlate_ifreq(const float *samples_ifreq, const std::vector<float>& ideal_chirp, const uint32_t to_idx);

		void fine_sync(const gr_complex* in_samples, uint32_t bin_idx, int32_t search_space);
		float detect_preamble_autocorr(const gr_complex *samples, const uint32_t window);
		float detect_downchirp(const gr_complex *samples, const uint32_t window);
		float detect_upchirp(const gr_complex *samples, const uint32_t window, int32_t *index);
		float sliding_norm_cross_correlate_upchirp(const float *samples_ifreq, const uint32_t window, int32_t *index);
		uint32_t max_frequency_gradient_idx(const gr_complex *samples);
		void deinterleave(const uint32_t ppm);
		void decode(const bool is_header);
		bool demodulate(const gr_complex *samples, const bool is_header);
		void deshuffle(const uint8_t *shuffle_pattern, const bool is_header);
		void dewhiten(const uint8_t *prng);
		void hamming_decode(bool is_header);
		void hamming_decode_soft(bool is_header);
		void extract_data_only(bool is_header);
		void set_abs_threshold(const float threshold);

		void print_complex_pointers(const gr_complex* input, int size);
		uint32_t get_d_samples_per_symbol();
		int32_t get_d_fine_sync();
		uint32_t get_d_delay_after_sync();
		uint32_t get_d_decim_factor();
		uint32_t get_d_number_of_bins();
		std::vector<uint8_t> get_d_decoded();
		void set_d_phdr(std::vector<uint8_t> d_decoded);
		void clear_d_decoded();
		void set_d_payload_length(uint32_t value);
		loraphy_header_t get_d_phdr();


		std::string logic_DETECT(const gr_complex * input);
		std::string logic_SYNC(const gr_complex * input, int * i);
		std::string logic_FIND_SFD(const gr_complex * input);
		std::string logic_PAUSE();
		std::string logic_DECODE_HEADER(const gr_complex * input);
		std::string logic_DECODE_PAYLOAD(const gr_complex * input);
		// void logic_STOP(const gr_complex * input, DecoderState * d_state);

		int work(const gr_complex * input);
};

#endif /* DECODER_H_ */ 