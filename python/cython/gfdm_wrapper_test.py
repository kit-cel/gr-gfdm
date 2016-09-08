import sys, os
sys.path.insert(0, os.path.abspath('/home/demel/src/gr-gfdm/python/build/lib.linux-x86_64-2.7/'))
print sys.path

import cgfdm
import numpy as np
from pygfdm.gfdm_modulation import gfdm_modulate_block
from pygfdm.mapping import get_data_matrix
from pygfdm.utils import get_random_qpsk
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.cyclic_prefix import get_raised_cosine_ramp


def data_per_volume():
    '''
    This is a little funtime activity easter egg trying to estimate data content of the human genome.
    '''
    micro_sd_dimensions = [11e-3, 15e-3, 1e-3]  # according to spec
    micro_sd_volume = micro_sd_dimensions[0] * micro_sd_dimensions[1] * micro_sd_dimensions[2]
    micro_sd_bits = 200.e9  # use value for some large known mirco SD card size
    genome_volume = (0.34e-9 ** 3) * 6.e9  # maybe correct? http://hypertextbook.com/facts/1998/StevenChen.shtml
    genome_bits = 4.e6  # lossless compressed aomunt of data
    print 'micro SD card volume: ', micro_sd_volume, 'cubic meter'
    print 'genome volume: ', genome_volume, 'cubic meter'
    print 'micro SD data per volume', micro_sd_bits / micro_sd_volume, 'bits/(cubic meter)'
    print 'genome data per volume', genome_bits / genome_volume, 'bits/(cubic meter)'


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
    cpler = cgfdm.py_add_cyclic_prefix_cc(M * K, 4, 4, window_taps)
    print cpler.block_size()
    print cpler.frame_size()
    block = cpler.add_cyclic_prefix(in_buf)
    print np.shape(block)

    energy_detector = cgfdm.py_detect_frame_energy_kernel_cl(50.3, 8)
    print 'energy_detector: alpha', energy_detector.alpha(), ', average_len', energy_detector.average_len()
    n_alpha = 47.11
    energy_detector.set_alpha(n_alpha)
    print 'energy_detector: alpha', energy_detector.alpha(), 'expected after reset', n_alpha

    syms = np.concatenate((np.ones(20), np.ones(8) * 2 * n_alpha)).astype(dtype=np.complex64)
    print syms
    pos = energy_detector.detect_frame(syms)
    print 'frame pos:', pos

    # resource_mapping_test()
    modulator_test()
    data_per_volume()


if __name__ == '__main__':
    main()
