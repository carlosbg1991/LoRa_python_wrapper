# -*- Python -*- #
#
# Copyright 2020 Carlos Bocanegra.
# The Genesys Lab, Northeastern University, Boston, MA
#
#

import sys
import setproctitle as setPT
import lora as lora_gr  # from GNURadio library
from lora_id import *  # from InterDigital customized library
from utils import *
import optparse


def run_channelizer_resampler(channelizer, filename_in, filename_out):
    print "Running channelizer Resampler..."
    dataset = read_complex_array(filename_in)  # after channelizer
    len_dataset = len(dataset)

    init_idx = 0
    idx_chann_out = 0
    d_resampler_in_out_ratio = channelizer.get_resampler_in_out_ratio()
    d_batch = len_dataset - 1  # one round and done...

    my_gr_in = new_grcomplex(d_batch)  # Store samples in gr_complex* (or std::vector<float>*)
    my_gr_out = new_grcomplex(d_batch / d_resampler_in_out_ratio)  # Store samples in gr_complex* (or std::vector<float>*)
    py_chann_out = np.empty(len_dataset / d_resampler_in_out_ratio, dtype=complex)  # samples in python to file

    while init_idx + d_batch < len_dataset:
        for idx_small in range(d_batch):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr_in, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        channelizer.work_resampler(my_gr_in, my_gr_out, d_batch / d_resampler_in_out_ratio)
        init_idx = init_idx + d_batch

        for idx_small in range(d_batch / d_resampler_in_out_ratio):
            py_chann_out[idx_chann_out] = get_grcomplex(my_gr_out, idx_small)
            idx_chann_out = idx_chann_out + 1

    print "Writing output samples from the PyLoRa channelizer/resampler block in {0}".format(filename_out)
    write_complex_array(py_chann_out, filename_out)


def run_channelizer_xlating(channelizer, filename_in, filename_out):
    print "Running channelizer Xlating FIR Filter..."

    dataset = read_complex_array(filename_in)
    len_dataset = len(dataset)

    init_idx = 0
    idx_chann_out = 0
    d_batch = len_dataset - 1  # one round and done...

    my_gr_in = new_grcomplex(d_batch)  # Store samples in gr_complex* (or std::vector<float>*)
    my_gr_out = new_grcomplex(d_batch)  # Store samples in gr_complex* (or std::vector<float>*)
    py_xlating_out = np.empty(len_dataset, dtype=complex)

    while init_idx + d_batch < len_dataset:

        for idx_small in range(d_batch):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr_in, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        channelizer.work_freq_xlating_fir_filter(my_gr_in, my_gr_out, d_batch)

        init_idx = init_idx + d_batch


        for idx_small in range(d_batch):
            py_xlating_out[idx_chann_out] = get_grcomplex(my_gr_out, idx_small)
            idx_chann_out = idx_chann_out + 1

    print "Writing output samples from the PyLoRa channelizer/Xlating_FIR block in {0}".format(filename_out)
    write_complex_array(py_xlating_out, filename_out)


def run_decoder(decoder, filename_in):
    print "Running decoder..."

    dataset = read_complex_array(filename_in)
    len_dataset = len(dataset)

    init_idx = 0
    d_samples_per_symbol = decoder.get_d_samples_per_symbol()

    my_gr = new_grcomplex(2 * d_samples_per_symbol)  # Store samples in gr_complex* (or std::vector<float>*)
    while init_idx + 2 * d_samples_per_symbol < len_dataset:

        for idx_small in range(2 * d_samples_per_symbol):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        idx_update = decoder.work(my_gr)

        init_idx = init_idx + idx_update


if __name__ == '__main__':

    print "Running script: {0}".format(" ".join([x for x in sys.argv]))
    setPT.setproctitle('lora-receiver')
    print "Process name:  " + str(get_proc_name())
    print "Process PID:  " + str(get_proc_pid(get_proc_name()))
    print "*****************************************"

    print "To assess the performance, run the following command on a side terminal"
    print ">> PID={0}".format(get_proc_pid(get_proc_name()))
    print ">> htop -p `pstree -p $PID | perl -ne 'push @t, /\((\d+)\)/g; END { print join \",\", @t }'`"
    print "*****************************************"

    parser = optparse.OptionParser()
    parser.add_option('-n', '--file', type="string", dest="FILENAME", help=" Input IQ file", default="../data/lora-99-100.sigmf-data")
    parser.add_option('-r', '--adcin', type="float", dest="ADC_IN", help="Sampling Rate in Hz", default=10e6)
    parser.add_option('-s', '--adcout', type="float", dest="ADC_OUT", help="Sampling rate after resampler in Hz", default=1e6)
    parser.add_option('-f', '--freq', type="float", dest="FREQ", help="Center Frequency in Hz",default=902.5e6)
    parser.add_option('-c', '--captfreq', type="string", dest="CAPTFREQ", help="Capture Frequency in Hz as a list, e.g., \"902.5e6,907.1e6\"", default="902.5e6")
    parser.add_option('-k', '--sf', type="int", dest="SPREADFACT", help="Spreading Factor", default=7)
    parser.add_option('-d', '--threshold', type="float", dest="THRESHOLD", help="Decision threshold for frame sync", default=0.01)
    options, args = parser.parse_args()
    file_in_resampler = options.FILENAME
    ADC_IN = options.ADC_IN
    ADC_OUT = options.ADC_OUT
    FREQ = options.FREQ
    CAPTFREQ = options.CAPTFREQ.split(",")
    CAPTFREQ = [float(x) for x in CAPTFREQ]
    SPREADFACT = options.SPREADFACT
    THRESHOLD = options.THRESHOLD

    print "Running LoRa receiver with configuration:"
    for k in options.__dict__.keys():
        print "\t{0}: {1}".format(k, options.__dict__[k])
    print "*****************************************"

    channelizer = channelizer(ADC_IN, ADC_OUT, FREQ, CAPTFREQ)
    decoder = decoder(ADC_OUT, SPREADFACT)
    decoder.set_abs_threshold(THRESHOLD)
    print "*****************************************"

    file_out_resampler = "../data/py_lora_output_resampler"
    file_out_xlating = "../data/py_lora_output_xlating"

    run_channelizer_resampler(channelizer, file_in_resampler, file_out_resampler)
    run_channelizer_xlating(channelizer, file_out_resampler, file_out_xlating)
    run_decoder(decoder, file_out_xlating)

    plot_complex_array(file_in_resampler, "INPUT DATASET")
    plot_complex_array("../data/py_lora_output_resampler", "AFTER RESAMPLER")
    plot_complex_array("../data/py_lora_output_xlating", "AFTER XLATING FIR FILTER")
