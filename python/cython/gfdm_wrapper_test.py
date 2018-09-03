import sys, os
sys.path.insert(0, os.path.abspath('/home/demel/src/gr-gfdm/python/build/lib.linux-x86_64-2.7/'))
print sys.path

import cgfdm
import numpy as np
from pygfdm.gfdm_modulation import gfdm_modulate_block
from pygfdm.mapping import get_data_matrix
from pygfdm.utils import get_random_qpsk, calculate_awgn_noise_variance, get_complex_noise_vector
from pygfdm.filters import get_frequency_domain_filter
from pygfdm.cyclic_prefix import get_raised_cosine_ramp
from pygfdm.synchronization import find_frame_start, simplified_sync_algo
from pygfdm.preamble import mapped_preamble


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
    per_timeslot = True
    mapper = cgfdm.py_resource_mapper_kernel_cc(timeslots, fft_len, active, smap, per_timeslot)
    demapper = cgfdm.py_resource_demapper_kernel_cc(timeslots, fft_len, active, smap, per_timeslot)

    assert mapper.input_vector_size() == demapper.output_vector_size()
    assert mapper.output_vector_size() == demapper.input_vector_size()

    N = 10
    for i in range(N):
        data = np.random.randn(2, active * timeslots)
        data = data[0] + 1j * data[1]
        data = data.astype(dtype=np.complex64)
        tx = mapper.map_to_resources(data)
        rd = demapper.demap_from_resources(tx, len(data))
        assert len(data) == len(rd)
        assert np.all(np.abs(data - rd) < 1e-10)


def cp_test():
    M = 5
    K = 16

    window_taps = get_raised_cosine_ramp(4, M * K + 4)
    cpler = cgfdm.py_add_cyclic_prefix_cc(M * K, 4, 0, 4, window_taps)
    print cpler.block_size()
    print cpler.frame_size()
    in_buf = get_random_qpsk(M * K, dtype=np.complex64)
    block = cpler.add_cyclic_prefix(in_buf)
    print np.shape(block)


def energy_detector_test():
    energy_detector = cgfdm.py_detect_frame_energy_kernel_cl(50.3, 8)
    print 'energy_detector: alpha', energy_detector.alpha(), ', average_len', energy_detector.average_len()
    n_alpha = 47.11
    energy_detector.set_alpha(n_alpha)
    print 'energy_detector: alpha', energy_detector.alpha(), 'expected after reset', n_alpha
    syms = np.concatenate((np.ones(20), np.ones(8) * 2 * n_alpha)).astype(dtype=np.complex64)
    print syms
    pos = energy_detector.detect_frame(syms)
    print 'frame pos:', pos


def modulate_gfdm_frame(tx_symbols, params):
    mapper = cgfdm.py_resource_mapper_kernel_cc(params['timeslots'], params['fft_len'], params['active_subcarriers'],
                                               params['subcarrier_map'])
    taps = get_frequency_domain_filter(params['prototype_type'], params['prototype_alpha'], params['timeslots'],
                                       params['fft_len'], params['overlap'])
    kernel = cgfdm.py_modulator_kernel_cc(params['timeslots'], params['fft_len'], params['overlap'], taps)
    win_filt = get_raised_cosine_ramp(params['filter_len'], params['timeslots'] * params['fft_len'] + params['cp_len'])
    cpler = cgfdm.py_add_cyclic_prefix_cc(params['timeslots'] * params['fft_len'], params['cp_len'], params['filter_len'],
                                         win_filt)

    syms = mapper.map_to_resources(tx_symbols)
    frame = kernel.modulate(syms.flatten())
    frame = cpler.add_cyclic_prefix(frame)
    return frame


def preamble_sync_test():
    seed = 4711
    timeslots = 9
    subcarriers = 32
    active_subcarriers = 20
    cp_len = 16
    filter_len = 8
    overlap = 2
    subcarrier_map = np.concatenate((np.arange(0, active_subcarriers // 2), np.arange(subcarriers - active_subcarriers // 2, subcarriers)))
    print subcarrier_map

    snr_db = 5.
    offset = 1000
    preamble, x_preamble = mapped_preamble(seed, 'rrc', .5, active_subcarriers, subcarriers, subcarrier_map, overlap,
                           cp_len, filter_len)

    noise_variance = calculate_awgn_noise_variance(preamble, snr_db)
    rx = get_complex_noise_vector(len(preamble) + 2 * offset, noise_variance)
    rx[offset:offset + len(preamble)] += preamble
    print('frame offset', offset, 'with cp_len', cp_len)
    xc = np.correlate(rx, x_preamble, 'valid')
    # xc *= np.max(np.abs(rx)) / np.max(np.abs(xc))
    # xc /= calculate_average_signal_energy(tx)
    peak = np.argmax(np.abs(xc))
    print('frame located @{}'.format(peak))

    nc, cfo, abs_corr_vals, corr_vals, napcc, apcc = find_frame_start(rx, x_preamble, subcarriers,
                                                                      cp_len)
    print('find frame start res @{} and cfo: {}'.format(nc, cfo))
    s_nm, s_cfo, simple_ac, phase = simplified_sync_algo(rx, x_preamble, subcarriers, cp_len)
    print('simple algo res @{} and cfo: {}'.format(s_nm, s_cfo))
    sync_kernel = cgfdm.py_auto_cross_corr_multicarrier_sync_cc(subcarriers, cp_len, x_preamble)
    kpos, kcfo = sync_kernel.detect_frame(rx.astype(dtype=np.complex64))

    print('kernel: nc={}, cfo={}'.format(kpos, kcfo))

    assert(kpos == s_nm)
    assert(np.abs(s_cfo - kcfo) < 1e-7)


def equalize(rx_frame, H):
    F = np.fft.fft(rx_frame)
    E = F / H
    return np.fft.ifft(E)


def estimate_preamble_aided_channel(rx, preamble, smap):
    P = np.fft.fft(preamble)
    return np.fft.fft(rx)[smap] / P[smap]


def interpolate_subcarriers(H, fft_len, smap):
    t = np.zeros(fft_len, dtype=np.complex)
    t[smap] = H
    return t


def interpolate_channel(est_H, frame_len, fft_len, cp_len, smap):
    t = interpolate_subcarriers(est_H, fft_len, smap)
    est_time = np.fft.ifft(t)
    est_time = est_time[0:cp_len]
    H_frame = np.fft.fft(est_time, frame_len)
    return H_frame



def main():
    np.set_printoptions(precision=2, suppress=True)
    # data_per_volume()
    err_margin = 1e-5
    M = 5
    K = 16
    L = 2

    cp_test()
    resource_mapping_test()
    modulator_test()
    energy_detector_test()
    preamble_sync_test()





if __name__ == '__main__':
    main()
