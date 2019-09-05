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

    def filter_taps(self):
        return np.array(self.kernel.filter_taps())

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

    def __cinit__(self, int block_len, int cp_len, int cs_len, int ramp_len, np.ndarray window_taps):
        self.kernel = new gfdm_interface.add_cyclic_prefix_cc(block_len, cp_len, cs_len, ramp_len, window_taps)

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
    def __cinit__(self, int timeslots, int subcarriers, int active_subcarriers, np.ndarray subcarrier_map, bool per_timeslot=True):
        self.kernel = new gfdm_interface.resource_mapper_kernel_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, per_timeslot)

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


cdef class py_resource_demapper_kernel_cc:
    cdef gfdm_interface.resource_demapper_kernel_cc* kernel
    def __cinit__(self, int timeslots, int subcarriers, int active_subcarriers, np.ndarray subcarrier_map, bool per_timeslot=True):
        self.kernel = new gfdm_interface.resource_demapper_kernel_cc(timeslots, subcarriers, active_subcarriers, subcarrier_map, per_timeslot)

    def __del__(self):
        del self.kernel

    def input_vector_size(self):
        return self.kernel.input_vector_size()

    def output_vector_size(self):
        return self.kernel.output_vector_size()

    def generic_work(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf, int noutput_size):
        self.kernel.generic_work(<float complex*> outbuf.data, <float complex*> inbuf.data, noutput_size)

    def demap_from_resources(self, np.ndarray[np.complex64_t, ndim=1] samples, noutput_size):
        if samples.shape[0] != self.input_vector_size():
            raise ValueError("Size of input array MUST be equal to timeslots * subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((noutput_size,), dtype=np.complex64)
        self.generic_work(res, samples, noutput_size)
        return res


cdef class py_receiver_kernel_cc:
    cdef gfdm_interface.receiver_kernel_cc* kernel
    def __cinit__(self, int timeslots, int subcarriers, int overlap, np.ndarray freq_taps):
        self.kernel = new gfdm_interface.receiver_kernel_cc(timeslots, subcarriers, overlap, freq_taps)

    def __del__(self):
        del self.kernel

    def block_size(self):
        return self.kernel.block_size()

    def filter_taps(self):
        return np.array(self.kernel.filter_taps())

    def ic_filter_taps(self):
        return np.array(self.kernel.ic_filter_taps())

    def generic_work(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf):
        self.kernel.generic_work(<float complex*> outbuf.data, <float complex*> inbuf.data)

    def generic_work_equalize(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf, np.ndarray[np.complex64_t, ndim=1] channel_estimate):
        self.kernel.generic_work_equalize(<float complex*> outbuf.data, <float complex*> inbuf.data, <float complex*> channel_estimate.data)

    def cpp_fft_filter_downsample(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf):
        self.kernel.fft_filter_downsample(<float complex*> outbuf.data, <float complex*> inbuf.data)

    def cpp_transform_subcarriers_to_td(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] inbuf):
        self.kernel.transform_subcarriers_to_td(<float complex*> outbuf.data, <float complex*> inbuf.data)

    def cpp_cancel_sc_interference(self, np.ndarray[np.complex64_t, ndim=1] outbuf, np.ndarray[np.complex64_t, ndim=1] td_inbuf, np.ndarray[np.complex64_t, ndim=1] fd_inbuf):
        self.kernel.cancel_sc_interference(<float complex*> outbuf.data, <float complex*> td_inbuf.data, <float complex*> fd_inbuf.data);

    def demodulate(self, np.ndarray[np.complex64_t, ndim=1] samples):
        if samples.shape[0] != self.block_size():
            raise ValueError("CGFDM: Size of input array MUST be equal to timeslots * active_subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.generic_work(res, samples)
        return res

    def demodulate_equalize(self, np.ndarray[np.complex64_t, ndim=1] samples, np.ndarray[np.complex64_t, ndim=1] channel_estimate):
        if samples.shape[0] != self.block_size() or channel_estimate.shape[0] != self.block_size():
            raise ValueError("CGFDM: Size of input array MUST be equal to timeslots * active_subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.generic_work_equalize(res, samples, channel_estimate)
        return res

    def fft_filter_downsample(self, np.ndarray[np.complex64_t, ndim=1] samples):
        if samples.shape[0] != self.block_size():
            raise ValueError("CGFDM: Size of input array MUST be equal to timeslots * active_subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.cpp_fft_filter_downsample(res, samples)
        return res

    def fft_equalize_filter_downsample(self, np.ndarray[np.complex64_t, ndim=1] samples, np.ndarray[np.complex64_t, ndim=1] fd_equalizer_taps):
        if samples.size != self.block_size():
            raise ValueError("CGFDM: Size of input array ({}) MUST be equal to timeslots * active_subcarriers ({})!".format(samples.size, self.block_size()))
        if fd_equalizer_taps.size != self.block_size():
            raise ValueError("CGFDM: Size of equalizer array ({}) MUST be equal to timeslots * active_subcarriers ({})!".format(fd_equalizer_taps.size, self.block_size()))
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.kernel.fft_equalize_filter_downsample(<float complex*> res.data, <float complex*> samples.data, <float complex*> fd_equalizer_taps.data);
        return res

    def transform_subcarriers_to_td(self, np.ndarray[np.complex64_t, ndim=1] samples):
        if samples.shape[0] != self.block_size():
            raise ValueError("CGFDM: Size of input array MUST be equal to timeslots * active_subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.cpp_transform_subcarriers_to_td(res, samples)
        return res

    def cancel_sc_interference(self, np.ndarray[np.complex64_t, ndim=1] td_in, np.ndarray[np.complex64_t, ndim=1] fd_in):
        if td_in.shape[0] != self.block_size():
            raise ValueError("CGFDM: Size of input array MUST be equal to timeslots * active_subcarriers!")
        if td_in.shape[0] != fd_in.shape[0]:
            raise ValueError("CGFDM: Size of BOTH input arrays MUST be equal to timeslots * active_subcarriers!")
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.block_size(),), dtype=np.complex64)
        self.cpp_cancel_sc_interference(res, td_in, fd_in)
        return res


cdef class py_preamble_channel_estimator_cc:
    cdef gfdm_interface.preamble_channel_estimator_cc* kernel

    def __cinit__(self, int timeslots, int fft_len, int active_subcarriers, is_dc_free, which_estimator, np.ndarray preamble):
        self.kernel = new gfdm_interface.preamble_channel_estimator_cc(timeslots, fft_len, active_subcarriers, is_dc_free, which_estimator, preamble.astype(dtype=np.complex64))

    def __del__(self):
        del self.kernel

    def preamble_filter_taps(self):
        return np.array(self.kernel.preamble_filter_taps())

    def estimate_preamble_channel(self, np.ndarray[np.complex64_t, ndim=1] rx_preamble):
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.kernel.fft_len(),), dtype=np.complex64)
        self.kernel.estimate_preamble_channel(<float complex*> res.data, <float complex*> rx_preamble.data)
        return res

    def filter_preamble_estimate(self, np.ndarray[np.complex64_t, ndim=1] estimate):
        res_len = self.kernel.active_subcarriers()
        if self.kernel.is_dc_free():
            res_len += 1
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((res_len,), dtype=np.complex64)
        self.kernel.filter_preamble_estimate(<float complex*> res.data, <float complex*> estimate.data)
        return res

    def interpolate_frame(self, np.ndarray[np.complex64_t, ndim=1] estimate):
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.kernel.fft_len() * self.kernel.timeslots(),), dtype=np.complex64)
        self.kernel.interpolate_frame(<float complex*> res.data, <float complex*> estimate.data)
        return res

    def estimate_frame(self, np.ndarray[np.complex64_t, ndim=1] rx_preamble):
        cdef np.ndarray[np.complex64_t, ndim=1] res = np.zeros((self.kernel.fft_len() * self.kernel.timeslots(),), dtype=np.complex64)
        self.kernel.estimate_frame(<float complex*> res.data, <float complex*> rx_preamble.data)
        return res
