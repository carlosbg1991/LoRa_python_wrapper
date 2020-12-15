import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_complex_array(filePath):
    f = open(filePath, 'rb')  # create file descriptor
    x = np.fromfile(f, dtype=np.float32, count=-1)  # read data into long array
    f.close()  # close file descriptor
    y_real = x[0::2]  # real values
    y_imag = x[1::2]  # imag values
    # Reconstruct the original complex array
    res = y_real + 1j * y_imag
    return res

def write_complex_array(data, filePath):
    y = np.empty(2*len(data), dtype=np.float32)
    y[0::2] = data.real
    y[1::2] = data.imag
    y.astype(np.float32).tofile(filePath)

def plot_complex_array(filePath,title = "my plot"):
    data = read_complex_array(filePath)
    plt.figure()
    plt.plot(np.real(data), linewidth=0.5)
    plt.plot(np.imag(data), linewidth=0.5)
    plt.title(title)
    plt.show()

def plot_complex_array_compare(filePath1, filePath2):
    data1 = read_complex_array(filePath1)
    data2 = read_complex_array(filePath2)
    plt.figure()
    plt.subplot(211)
    plt.plot(np.real(data1), linewidth=0.5)
    plt.plot(np.real(data2), linewidth=0.5)
    plt.subplot(212)
    plt.plot(np.imag(data1), linewidth=0.5)
    plt.plot(np.imag(data2), linewidth=0.5)
    plt.show()

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