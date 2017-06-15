#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
import scipy.signal as signal
from fractions import gcd

import synchronization as sync
import preamble
import mapping
import converter
import utils
import gfdm_modulation as gmod
import cyclic_prefix
import gfdm_receiver
import filters
import validation_utils

import matplotlib.pyplot as plt

import cgfdm


def synchronize(frame, ref_frame, x_preamble, fft_len, cp_len):
    ac = sync.auto_correlate_signal(frame, fft_len)
    nm = np.argmax(np.abs(ac))
    # cfo = 2. * np.angle(ac[nm-5:nm+5]) / (2. * np.pi)
    # print(cfo)
    # cfo = np.sum(cfo) / len(cfo)
    # cfo = 2 * np.angle(ac[nm]) / (2. * np.pi)
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    print('CFO:', cfo, cfo * 3.125e6 / fft_len)

    phase_inc = sync.cfo_to_phase_increment(-cfo, fft_len)
    wave = sync.complex_sine(phase_inc, len(frame), 0.0)
    frame *= wave

    ac = sync.auto_correlate_signal(frame, fft_len)
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    print('CFO:', cfo, cfo * 3.125e6 / fft_len)

    xc = np.correlate(frame, x_preamble, 'valid')
    cc = sync.multiply_valid(np.abs(ac), np.abs(xc))
    nc = np.argmax(np.abs(cc))
    print('correlation frame start:', nc)
    sample_nc = nc - fft_len
    print('sample frame start:     ', sample_nc)
    p_len = cp_len + 2 * fft_len + 16 + cp_len
    print('data frame start:       ', sample_nc + p_len)
    phase = np.angle(xc[nc])
    print('phase:', phase)
    frame *= np.exp(-1j * phase)

    ref_e = utils.calculate_signal_energy(x_preamble)
    rx_e = utils.calculate_signal_energy(frame[nc:nc + len(x_preamble)])
    agc_factor = np.sqrt(ref_e / rx_e)
    print('AGC values:', ref_e, rx_e, agc_factor)
    frame *= agc_factor

    # plt.plot(np.abs(ref_frame))
    # plt.plot(np.abs(frame))
    # plt.plot(np.abs(ac))
    # # # plt.plot(np.abs(xc))
    # plt.plot(cc)
    # # plt.axvline(nc, color='y')
    # plt.show()
    # return
    sframe = frame[sample_nc:sample_nc + len(ref_frame)]
    return sframe


def preamble_estimate(rx_preamble, x_preamble, fft_len):
    e0 = np.fft.fft(rx_preamble[0:fft_len]) / np.fft.fft(x_preamble[0:fft_len])
    e1 = np.fft.fft(rx_preamble[fft_len:]) / np.fft.fft(x_preamble[fft_len:])
    H = (e0 + e1) / 2
    return H, e0, e1


def estimate_frame(r_frame, ref_frame, fft_len, timeslots):
    t = np.arange(10)
    print(t)
    t= np.reshape(t, (-1, 2))
    print(t)
    rx = np.reshape(r_frame, (timeslots, fft_len))
    tx = np.reshape(ref_frame, (timeslots, fft_len))
    e = np.zeros(np.shape(rx), dtype=np.complex)
    for i in range(timeslots):
        r = rx[i]
        t = tx[i]
        e[i] = np.fft.fft(r) / np.fft.fft(t)
    return e


def main():
    np.set_printoptions(precision=2, linewidth=150)
    alpha = .2
    active_subcarriers = 52
    timeslots = 9
    fft_len = 64
    cp_len = fft_len // 2
    cs_len = cp_len // 2
    subcarrier_map = mapping.get_subcarrier_map(fft_len, active_subcarriers, dc_free=True)
    print(subcarrier_map)
    # filename = '/home/demel/iq_samples/gfdm_replay.dat'
    # frame = converter.load_gr_iq_file(filename)[34000:40000]
    # plt.semilogy(*signal.welch(frame))

    ref_frame, modulated_frame, x_preamble, data, freq_filter_taps = validation_utils.generate_reference_frame(timeslots, fft_len, active_subcarriers, cp_len, cs_len, alpha)
    rx_kernel = cgfdm.py_receiver_kernel_cc(timeslots, fft_len, 2, freq_filter_taps)
    demapper = cgfdm.py_resource_demapper_kernel_cc(timeslots, fft_len, active_subcarriers, subcarrier_map, True)

    ref_syms64 = modulated_frame.astype(dtype=np.complex64)
    # fd_syms = rx_kernel.fft_filter_downsample(ref_syms64)
    # txd2 = rx_kernel.transform_subcarriers_to_td(fd_syms)

    txd = rx_kernel.demodulate(ref_syms64)
    txd = demapper.demap_from_resources(txd, len(data))
    plt.scatter(txd.real, txd.imag, color='r')
    # plt.semilogy(*signal.welch(ref_frame))
    plt.show()
    return


    sframe = synchronize(frame, ref_frame, x_preamble, fft_len, cp_len)
    # xc = np.correlate(frame, ref_frame)
    # nc = np.argmax(np.abs(xc))
    # print(np.abs(xc[nc]), np.angle(xc[nc]))
    # plt.plot(np.angle(ref_frame / sframe))
    # plt.show()
    # return

    rx_preamble = sframe[cp_len:cp_len + 2 * fft_len]
    H, e0, e1 = preamble_estimate(rx_preamble, x_preamble, fft_len)
    H_est = np.repeat(H, timeslots)
    # H = np.fft.fftshift(H)
    # e0 = np.fft.fftshift(e0)
    # e1 = np.fft.fftshift(e1)
    # plt.plot(np.angle(e0))
    # plt.plot(np.angle(e1))
    # plt.plot(np.angle(H))
    # plt.grid()
    # plt.show()
    # return

    f_start = cp_len + 2 * fft_len + cs_len
    d_start = f_start + cp_len

    data_syms = sframe[d_start:d_start + fft_len * timeslots]
    ref_syms = ref_frame[d_start:d_start + fft_len * timeslots]


    rx_syms64 = data_syms.astype(dtype=np.complex64)
    fd_syms = rx_kernel.fft_filter_downsample(rx_syms64)
    fd_syms = np.reshape(fd_syms, (fft_len, timeslots))
    fd_eq = np.zeros(np.shape(fd_syms), dtype=np.complex64)
    for i in range(fft_len):
        fd_eq[i] = fd_syms[i] * np.conj(H[i])
    fd_eq = fd_eq.flatten()


    r_data = rx_kernel.transform_subcarriers_to_td(fd_eq)


    # f_estimate = estimate_frame(data_syms, ref_syms, fft_len, timeslots)

    ref_syms64 = ref_syms.astype(dtype=np.complex64)
    fd_syms = rx_kernel.fft_filter_downsample(ref_syms64)
    txd2 = rx_kernel.transform_subcarriers_to_td(fd_syms)

    txd = rx_kernel.demodulate(ref_syms64)
    t_data = gfdm_receiver.gfdm_demodulate_fft(ref_syms, .5, timeslots, fft_len, 2, sic_rounds=0)

    sc0d = t_data[timeslots:2*timeslots]
    t_data_m = np.reshape(t_data, (fft_len, timeslots))
    print(np.abs(sc0d - t_data_m[1]) < 1e-8)


    # r_data = gfdm_receiver.gfdm_demodulate_fft(data_syms, .5, timeslots, fft_len, 2, sic_rounds=0)
    r_data_m = np.reshape(r_data, (fft_len, timeslots))

    # for i in subcarrier_map[0:10]:
    #     # plt.scatter(t_data_m[i].real, t_data_m[i].imag)
    #     # plt.scatter(r_data_m[i].real, r_data_m[i].imag, color='r')
    #     a = np.angle(r_data_m[i] / t_data_m[i])
    #     plt.plot(a)
    #     print(a)
    # plt.show()
    # return
    # tt = rx_kernel.demodulate(data_syms.astype(dtype=np.complex64))



    d = mapping.demap_from_waveform_resource_grid(t_data, fft_len, subcarrier_map)

    r_data = mapping.demap_from_waveform_resource_grid(r_data, fft_len, subcarrier_map)
    plt.scatter(r_data.real, r_data.imag)
    plt.scatter(t_data.real, t_data.imag, color='r')
    # plt.scatter(d.real, d.imag, color='g')
    # plt.scatter(tt.real, tt.imag, color='m')

    print(d[0:5])
    b = utils.demodulate_qpsk(d)
    print(b[0:10])
    br = utils.demodulate_qpsk(r_data)
    print(br[0:10])
    print(np.sum(br == b) / len(b))

    # ref_frame = converter.convert_to_cf64(ref_frame)
    # ref_frame.tofile('/home/demel/iq_samples/gfdm_reference_frame.dat')

    plt.show()


if __name__ == '__main__':
    main()
