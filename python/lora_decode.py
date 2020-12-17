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
    setPT.setproctitle('lora-tester-py')
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

    my_decoder = decoder(1e6, 7)
    my_decoder.set_abs_threshold(0.01)

    # Read from dataset (returns numpy)
    # dataset = read_complex_array("../data/lora_output_chann")  # after channelizer
    dataset = read_complex_array("../data/py_lora_output_xlating")  # after channelizer
    # dataset = read_complex_array("/home/carlos/git/gr-lora-0.6.2/misc/data/gr_record_out_channelizer")  # after channelizer
    # dataset = read_complex_array("/home/carlos/git/gr-lora-0.6.2/misc/data/gr_record_out_resampler")  # after channelizer
    print(len(dataset))
    print(type(dataset))

    MAC_CRC_SIZE = 2  # from utilities function
    init_idx = 0
    sync_idx = 0
    d_corr_fails = 0
    d_samples_per_symbol = my_decoder.get_d_samples_per_symbol()
    d_decim_factor = my_decoder.get_d_decim_factor()
    d_number_of_bins = my_decoder.get_d_number_of_bins()
    d_delay_after_sync = my_decoder.get_d_delay_after_sync()

    my_gr = new_grcomplex(2 * d_samples_per_symbol)  # Store samples in gr_complex* (or std::vector<float>*)

    while (init_idx + 2 * d_samples_per_symbol < len(dataset)):
        # plt.clf()
        # plt.plot(np.real(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.plot(np.imag(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.show()

        for idx_small in range(2 * d_samples_per_symbol):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        idx_update = my_decoder.work(my_gr)
        init_idx = init_idx + idx_update

    plt.plot(np.real(dataset), linewidth=0.5)
    plt.plot(np.imag(dataset), linewidth=0.5)
    plt.show()
