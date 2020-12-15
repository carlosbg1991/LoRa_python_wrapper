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
    dataset = read_complex_array("../data/py_lora_output_resampler")  # before before XLATING , after RESAMPLER (Python)
    # dataset = read_complex_array("/home/carlos/git/gr-lora-0.6.2/misc/data/gr_record_out_resampler")  # before before XLATING , after RESAMPLER (GR)
    fileName_out = '../data/py_lora_output_xlating'
    len_dataset = len(dataset)
    type_dataset = type(dataset)
    print(len(dataset))
    print(type(dataset))

    my_channelizer = channelizer(10e6, 1e6, 902.5e6, [902.5e6])

    init_idx = 0
    idx_chann_out = 0
    d_corr_fails = 0
    d_batch = len_dataset-1

    my_gr_in = new_grcomplex(d_batch)  # Store samples in gr_complex* (or std::vector<float>*)
    my_gr_out = new_grcomplex(d_batch)  # Store samples in gr_complex* (or std::vector<float>*)
    py_xlating_out = np.empty(len_dataset, dtype=complex)

    while (init_idx + d_batch < len_dataset):
        # plt.clf()
        # plt.plot(np.real(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.plot(np.imag(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.show()

        for idx_small in range(d_batch):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr_in, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        my_channelizer.work_freq_xlating_fir_filter(my_gr_in, my_gr_out, d_batch)

        init_idx = init_idx + d_batch

        for idx_small in range(d_batch):
            py_xlating_out[idx_chann_out] = get_grcomplex(my_gr_out, idx_small)
            idx_chann_out = idx_chann_out + 1

    print "Length input (from USRP) samples {0}".format(len(dataset))
    print "Actual Length output (of Channelizer) samples {0}".format(len(py_xlating_out))
    print "Writing output samples from the PyLoRa channelizer block in {0}".format(fileName_out)
    write_complex_array(py_xlating_out, fileName_out)

    plt.figure()
    plt.plot(np.real(py_xlating_out), linewidth=0.5)
    plt.plot(np.imag(py_xlating_out), linewidth=0.5)