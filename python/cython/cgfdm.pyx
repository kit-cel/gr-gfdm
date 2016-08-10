cimport gfdm_interface
import numpy as np
import time
cimport numpy as np
from libcpp cimport bool


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


cdef class py_resource_mapper_kernel_cc:
    cdef gfdm_interface.resource_mapper_kernel_cc* kernel
    def __cinit__(self, int active_subcarriers, int fft_len, int timeslots, np.ndarray subcarrier_map, bool per_timeslot=True):
        self.kernel = new gfdm_interface.resource_mapper_kernel_cc(active_subcarriers, fft_len, timeslots, subcarrier_map, per_timeslot)

    def __del__(self):
        del self.kernel

    def input_vector_size(self):
        return self.kernel.input_vector_size()

    def output_vector_size(self):
        return self.kernel.output_vector_size()

    def generic_work(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf, int input_size):
        self.kernel.generic_work(<float complex*> outbuf.data, <float complex*> inbuf.data, input_size)

    def map_to_resources(self, np.ndarray[np.complex64_t, ndim=1] samples):
        if samples.shape[0] > self.input_vector_size():
            raise ValueError("Size of input array MUST be smaller or equal to timeslots * active_subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.output_vector_size(),), dtype=np.complex64)
        self.generic_work(res, samples, samples.shape[0])

        return res




