import sys, os
sys.path.insert(0, os.path.abspath('/home/demel/src/gr-gfdm/python/build/lib.linux-x86_64-2.7/'))
print sys.path

import cgfdm
import numpy as np
from gfdm.pygfdm.gfdm_modulation import gfdm_modulate_block
from gfdm.pygfdm.mapping import get_data_matrix
from gfdm.pygfdm.utils import get_random_qpsk
from gfdm.pygfdm.filters import get_frequency_domain_filter
from gfdm.pygfdm.cyclic_prefix import get_raised_cosine_ramp


def modulator_test():
    fft_len = 128
    timeslots = 205
    overlap = 2
    taps = get_frequency_domain_filter('rrc', .5, timeslots, fft_len, overlap)
    kernel = cgfdm.py_modulator_kernel_cc(timeslots, fft_len, overlap, taps)

    N = 10
    for i in range(N):
        data = np.random.randn(2, fft_len * timeslots)
        data = data[0] + 1j * data[1]
        data = data.astype(dtype=np.complex64)
        kernel.modulate(data)


def resource_mapping_test():
    active = 110
    fft_len = 128
    timeslots = 205
    smap = np.arange(active) + (active - active) // 2
    mapper = cgfdm.py_resource_mapper_kernel_cc(active, fft_len, timeslots, smap, True)

    N = 10
    for i in range(N):
        data = np.random.randn(2, active * timeslots)
        data = data[0] + 1j * data[1]
        data = data.astype(dtype=np.complex64)
        mapper.map_to_resources(data)



def main():
    np.set_printoptions(precision=2, suppress=True)
    err_margin = 1e-5
    M = 5
    K = 16
    L = 2
    taps = get_frequency_domain_filter('rrc', .5, M, K, L)
    kernel = cgfdm.py_modulator_kernel_cc(M, K, L, taps)
    print kernel.block_size()
    in_buf = get_random_qpsk(M * K).astype(dtype=np.complex64)
    out_buf = np.zeros(kernel.block_size(), dtype=np.complex64)
    print type(in_buf[0])
    kernel.generic_work(out_buf, in_buf)
    res = kernel.modulate(in_buf)
    # print np.reshape(res, (M, -1))

    D = get_data_matrix(in_buf, K, group_by_subcarrier=False)
    ref = gfdm_modulate_block(D, taps, M, K, L, False)
    # print np.reshape(ref, (M, -1))

    print np.all(np.abs(out_buf - res) < err_margin)
    print np.all(np.abs(ref - res) < err_margin)

    window_taps = get_raised_cosine_ramp(4, M * K + 4)
    cpler = cgfdm.py_add_cyclic_prefix_cc(4, 4, M * K, window_taps)
    print cpler.block_size()
    print cpler.frame_size()
    block = cpler.add_cyclic_prefix(in_buf)
    print np.shape(block)

    # resource_mapping_test()
    modulator_test()


if __name__ == '__main__':
    main()
