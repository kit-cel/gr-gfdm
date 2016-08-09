from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "modulator_kernel_cc.h" namespace "gr::gfdm":
    cdef cppclass modulator_kernel_cc:
        modulator_kernel_cc(int, int, int, vector[float complex]) except +
        int block_size()
        void generic_work(float complex* p_out, const float complex* p_in);

cdef extern from "add_cyclic_prefix_cc.h" namespace "gr::gfdm":
    cdef cppclass add_cyclic_prefix_cc:
        add_cyclic_prefix_cc(int, int, int, vector[float complex]) except +
        int block_size()
        int frame_size()
        void generic_work(float complex* p_out, const float complex* p_in);

cdef extern from "resource_mapper_kernel_cc.h" namespace "gr::gfdm":
    cdef cppclass resource_mapper_kernel_cc:
        resource_mapper_kernel_cc(int, int, int, vector[int], bool) except +
        int input_vector_size()
        int output_vector_size()
        void generic_work(float complex* p_out, const float complex* p_in, const int ninput_size);
