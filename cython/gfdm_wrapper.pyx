cimport gfdm_interface
import numpy as np
cimport numpy as np


cdef class py_modulator_kernel_cc:
    cdef gfdm_interface.modulator_kernel_cc* kernel

    def __cinit__(self, int timeslots, int subcarriers, int overlap, np.ndarray taps):
        self.kernel = new gfdm_interface.modulator_kernel_cc(timeslots, subcarriers, overlap, taps)

    def __del__(self):
        del self.kernel

    def block_size(self):
        return self.kernel.block_size()

    def generic_work(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf):
        self.kernel.generic_work(<float complex*> outbuf.data, <float complex*> inbuf.data)

    def modulate(self, np.ndarray[np.complex64_t, ndim=1] samples):
        if samples.shape[0] != self.block_size():
            raise ValueError("Size of input array MUST equal timeslots * subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.generic_work(res, samples)
        return res


cdef class py_add_cyclic_prefix_cc:
    cdef gfdm_interface.add_cyclic_prefix_cc* kernel

    def __cinit__(self, int ramp_len, int cp_len, int block_len, np.ndarray window_taps):
        self.kernel = new gfdm_interface.add_cyclic_prefix_cc(ramp_len, cp_len, block_len, window_taps)

    def __del__(self):
        del self.kernel

    def block_size(self):
        return self.kernel.block_size()

    def frame_size(self):
            return self.kernel.frame_size()

    def generic_work(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf):
        self.kernel.generic_work(<float complex*> outbuf.data, <float complex*> inbuf.data)

    def add_cyclic_prefix(self, np.ndarray[np.complex64_t, ndim=1] samples):
        if samples.shape[0] != self.block_size():
            raise ValueError("Size of input array MUST equal timeslots * subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.frame_size(),), dtype=np.complex64)
        self.generic_work(res, samples)
        return res




