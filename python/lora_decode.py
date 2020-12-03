import ctypes
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import setproctitle as setPT
import lora as lora_gr  # from GNURadio library
from lora_id import *  # from InterDigital customized library


def read_complex_array(filePath):
    f = open(filePath, 'rb')  # create file descriptor
    x = np.fromfile(f, dtype=np.float32, count=-1)  # read data into long array
    f.close()  # close file descriptor

    y_real = x[0::2]  # real values
    y_imag = x[1::2]  # imag values

    # Reconstruct the original complex array
    res = y_real + 1j * y_imag
    return res


if __name__ == '__main__':
    def set_procname(newname):
        from ctypes import cdll, byref, create_string_buffer
        libc = cdll.LoadLibrary('libc.so.6')  # Loading a 3rd party library C
        buff = create_string_buffer(len(newname) + 1)  # Note: One larger than the name (man prctl says that)
        buff.value = newname  # Null terminated string as it should be
        libc.prctl(15, byref(buff), 0, 0,
                   0)  # Refer to "#define" of "/usr/include/linux/prctl.h" for the misterious value 16 & arg[3..5] are zero as the man page says.

    def get_proc_name():
        from ctypes import cdll, byref, create_string_buffer
        libc = cdll.LoadLibrary('libc.so.6')
        buff = create_string_buffer(128)
        # 16 == PR_GET_NAME from <linux/prctl.h>
        libc.prctl(16, byref(buff), 0, 0, 0)
        return buff.value


    setPT.setproctitle('lora-tester-py')
    print "Process name:  " + str(get_proc_name())
    print "***********************************"

    # ----------------------------------------------
    # GDB ATTACH (DEBUGGING or performance monitoring)
    # ----------------------------------------------
    GDB_ATTACH = 1
    if (GDB_ATTACH):
        print ('Blocked waiting for GDB attach (pid = %d) ' % (os.getpid(),) + '. Press ENTER after GDB is attached.')
        sys.stdout.flush()
        raw_input()

    my_decoder = decoder(1e6, 7)
    my_decoder.set_abs_threshold(0.01)

    # Read from dataset (returns numpy)
    dataset = read_complex_array("../data/lora_output_chann")  # after channelizer
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
    STATUS = 'DETECT'  # DETECT, SYNC, FIND_SFD, DECODE_HEADER, DECODE_PAYLOAD

    while (init_idx + 2 * d_samples_per_symbol < len(dataset)):
        # plt.clf()
        # plt.plot(np.real(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.plot(np.imag(dataset[init_idx:init_idx + d_samples_per_symbol]), linewidth=0.5)
        # plt.show()

        for idx_small in range(2 * d_samples_per_symbol):
            idx_dataset = init_idx + idx_small
            set_grcomplex(my_gr, idx_small, complex(dataset[idx_dataset].real, dataset[idx_dataset].imag))

        if STATUS == 'DETECT':
            STATUS = my_decoder.logic_DETECT(my_gr)
            if STATUS == 'DETECT':
                init_idx = init_idx + d_samples_per_symbol
        elif STATUS == 'SYNC':
            i_ptr = create_int_ptr(0)
            STATUS = my_decoder.logic_SYNC(my_gr, i_ptr)
            sync_idx = get_int(i_ptr)
            init_idx = init_idx + sync_idx
        elif STATUS == 'FIND_SFD':
            STATUS = my_decoder.logic_FIND_SFD(my_gr)
            d_fine_sync = my_decoder.get_d_fine_sync()
            init_idx = init_idx + d_samples_per_symbol + d_fine_sync
        elif STATUS == "PAUSE":
            STATUS = my_decoder.logic_PAUSE()
            init_idx = init_idx + d_samples_per_symbol + d_delay_after_sync
        elif STATUS == "DECODE_HEADER":
            STATUS = my_decoder.logic_DECODE_HEADER(my_gr)
            init_idx = init_idx + d_samples_per_symbol + d_fine_sync
        elif STATUS == "DECODE_PAYLOAD":
            STATUS = my_decoder.logic_DECODE_PAYLOAD(my_gr)
            init_idx = init_idx + d_samples_per_symbol + d_fine_sync

    plt.plot(np.real(dataset), linewidth=0.5)
    plt.plot(np.imag(dataset), linewidth=0.5)
    plt.show()
