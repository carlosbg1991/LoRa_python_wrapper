# -*- Python -*- #
#
# Copyright 2020 Carlos Bocanegra.
# The Genesys Lab, Northeastern University, Boston, MA
#
#

import ctypes
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import setproctitle as setPT
import lora as lora_gr  # from GNURadio library
from lora_id import *  # from InterDigital customized library
from utils import *


if __name__ == '__main__':
    setPT.setproctitle('lora-tester-channelizer-py')
    print "Process name:  " + str(get_proc_name())
    print "***********************************"

    # ----------------------------------------------
    # GDB ATTACH (DEBUGGING or performance monitoring)
    # ----------------------------------------------
    GDB_ATTACH = 0
    if (GDB_ATTACH):
        print ('Blocked waiting for GDB attach (pid = %d) ' % (os.getpid(),) + '. Press ENTER after GDB is attached.')
        sys.stdout.flush()
        raw_input()

    # Read from dataset (returns numpy)
    dataset = read_complex_array("../data/lora-99-100.sigmf-data")  # before channelizer
    fileName_out = '../data/py_lora_output_resampler'
    len_dataset = len(dataset)
    type_dataset = type(dataset)
    print(len(dataset))
    print(type(dataset))

    my_channelizer = channelizer(10e6, 1e6, 902.5e6, [902.5e6])

    MAC_CRC_SIZE = 2  # from utilities function
    init_idx = 0
    idx_chann_out = 0
    d_corr_fails = 0
    d_resampler_in_out_ratio = my_channelizer.get_resampler_in_out_ratio()
    d_resampler_nTaps = my_channelizer.get_resampler_nTaps()
    d_resampler_samples_per_batch = 1000000
    d_batch = d_resampler_samples_per_batch*d_resampler_in_out_ratio + d_resampler_nTaps
    d_batch = len_dataset - 1

    my_gr_in = new_grcomplex(d_batch)  # Store samples in gr_complex* (or std::vector<float>*)
    my_gr_out = new_grcomplex(d_batch / d_resampler_in_out_ratio)  # Store samples in gr_complex* (or std::vector<float>*)

    my_gr_chann_out = new_grcomplex(len(dataset)/d_resampler_in_out_ratio)  # Store samples in gr_complex* (or std::vector<float>*)
    py_chann_out = np.empty(len_dataset/d_resampler_in_out_ratio, dtype=complex)

    while (init_idx + d_batch < len_dataset):
        # plt.clf()
        # plt.plot(np.real(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.plot(np.imag(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.show()

        for idx_small in range(d_batch):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr_in, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        my_channelizer.work_resampler(my_gr_in, my_gr_out, d_batch / d_resampler_in_out_ratio)
        init_idx = init_idx + d_batch

        for idx_small in range(d_batch / d_resampler_in_out_ratio):
            py_chann_out[idx_chann_out] = get_grcomplex(my_gr_out, idx_small)
            idx_chann_out = idx_chann_out + 1

    print "Length input (from USRP) samples {0}".format(len(dataset))
    print "Actual Length output (of Channelizer/resampler) samples {0}".format(len(py_chann_out))
    print "Writing output samples from the PyLoRa channelizer/resampler block in {0}".format(fileName_out)
    write_complex_array(py_chann_out, fileName_out)

    plt.figure()
    plt.plot(np.real(py_chann_out), linewidth=0.5)
    plt.plot(np.imag(py_chann_out), linewidth=0.5)
    plt.show()