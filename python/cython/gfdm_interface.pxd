from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "modulator_kernel_cc.h" namespace "gr::gfdm":
    cdef cppclass modulator_kernel_cc:
        modulator_kernel_cc(int, int, int, vector[float complex]) except +
        int block_size()
        vector[float complex] filter_taps()
        void generic_work(float complex* p_out, const float complex* p_in);

cdef extern from "add_cyclic_prefix_cc.h" namespace "gr::gfdm":
    cdef cppclass add_cyclic_prefix_cc:
        add_cyclic_prefix_cc(int, int, int, int, vector[float complex]) except +
        int block_size()
        int frame_size()
        void generic_work(float complex* p_out, const float complex* p_in);

cdef extern from "resource_mapper_kernel_cc.h" namespace "gr::gfdm":
    cdef cppclass resource_mapper_kernel_cc:
        resource_mapper_kernel_cc(int, int, int, vector[int], bool) except +
        int input_vector_size()
        int output_vector_size()
        void generic_work(float complex* p_out, const float complex* p_in, const int ninput_size);

cdef extern from "resource_demapper_kernel_cc.h" namespace "gr::gfdm":
    cdef cppclass resource_demapper_kernel_cc:
        resource_demapper_kernel_cc(int, int, int, vector[int], bool) except +
        int input_vector_size()
        int output_vector_size()
        void generic_work(float complex* p_out, const float complex* p_in, const int noutput_size);

cdef extern from "receiver_kernel_cc.h" namespace "gr::gfdm":
    cdef cppclass receiver_kernel_cc:
        receiver_kernel_cc(int, int, int, vector[float complex]) except +
        int block_size()
        vector[float complex] filter_taps()
        vector[float complex] ic_filter_taps()
        void generic_work(float complex* p_out, const float complex* p_in);
        void generic_work_equalize(float complex* p_out, const float complex* p_in, const float complex* f_eq_in)
        void fft_filter_downsample(float complex* p_out, const float complex* p_in);
        void fft_equalize_filter_downsample(float complex* p_out, const float complex* p_in, const float complex* f_eq_in);
        void transform_subcarriers_to_td(float complex* p_out, const float complex* p_in);
        void cancel_sc_interference(float complex* p_out, const float complex* p_td_in, const float complex* p_fd_in);

cdef extern from "preamble_channel_estimator_cc.h" namespace "gr::gfdm":
    cdef cppclass preamble_channel_estimator_cc:
        preamble_channel_estimator_cc(int, int, int, bool, int, vector[float complex]) except +
        int fft_len()
        int timeslots()
        int active_subcarriers()
        bool is_dc_free()
        vector[float] preamble_filter_taps()
        void estimate_preamble_channel(float complex* fd_preamble_channel, const float complex* rx_preamble)
        void filter_preamble_estimate(float complex* filtered, const float complex* estimate)
        void interpolate_frame(float complex* frame_estimate, const float complex* estimate)
        void estimate_frame(float complex* frame_estimate, const float complex* rx_preamble)
