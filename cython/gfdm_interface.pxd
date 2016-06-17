from libcpp.vector cimport vector

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
