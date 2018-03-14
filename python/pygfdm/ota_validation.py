#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
import scipy.signal as signal
from fractions import gcd
import sys, os
import time

sys.path.insert(0, os.path.abspath('/home/demel/src/gr-gfdm/examples/'))
import gfdm_file_sync

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

try:
import matplotlib.pyplot as plt
except ImportError:
    pass

import cgfdm


def kernel_synchronize_time(frame, x_preamble, fft_len, cp_len, frame_len):
    sync_kernel = cgfdm.py_auto_cross_corr_multicarrier_sync_cc(fft_len, cp_len, x_preamble)
    kac = sync_kernel.fixed_lag_auto_correlate(frame.astype(dtype=np.complex64))
    knm = sync_kernel.find_peak(np.abs(kac).astype(dtype=np.float32))
    xstart = np.maximum(knm - sync_kernel.subcarriers(), 0)
    xc_range = np.arange(xstart, xstart + 4 * sync_kernel.subcarriers())
    kxc = sync_kernel.cross_correlate_preamble(frame[xc_range].astype(dtype=np.complex64))
    kcc = np.abs(kac[xstart:xstart + len(kxc)]) * np.abs(kxc)
    knc = xstart + sync_kernel.find_peak(kcc.astype(dtype=np.float32))
    # kpos, kcfo = sync_kernel.detect_frame(frame.astype(dtype=np.complex64))
    sample_nc = knc - cp_len
    sframe = frame[sample_nc:sample_nc + frame_len]
    print(knc, knm, sample_nc, len(sframe))
    p_att = sync_kernel.calculate_preamble_attenuation(sframe[cp_len:cp_len + 2 * fft_len].astype(dtype=np.complex64))
    sframe *= (1. / p_att)
    return sframe


def synchronize_time(frame, ref_frame, x_preamble, fft_len, cp_len, samp_rate=12.e6):
    kframe = kernel_synchronize_time(frame, x_preamble, fft_len, cp_len, len(ref_frame))

    ac = sync.auto_correlate_signal(frame, fft_len)
    nm = np.argmax(np.abs(ac))
    print('AC start: ', nm)
    # cfo = 2 * np.angle(ac[nm]) / (2. * np.pi)
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    print('CFO:', cfo, cfo * samp_rate / fft_len)

    phase_inc = sync.cfo_to_phase_increment(-cfo, fft_len)
    wave = sync.complex_sine(phase_inc, len(frame), 0.0)
    # print(len(wave), len(frame))
    # frame *= wave

    ac = sync.auto_correlate_signal(frame, fft_len)
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    print('CFO:', cfo, cfo * samp_rate / fft_len)

    xc = np.correlate(frame, x_preamble, 'valid')
    cc = sync.multiply_valid(np.abs(ac), np.abs(xc))
    nc = np.argmax(np.abs(cc))
    print('correlation frame start:', nc)
    cfo = np.angle(ac[nc]) / (2. * np.pi)
    print('CFO:', cfo, cfo * samp_rate / fft_len)
    sample_nc = nc - cp_len
    print('sample frame start:     ', sample_nc)
    p_len = cp_len + 2 * fft_len + cp_len // 2 + cp_len
    print('data frame start:       ', sample_nc + p_len)
    phase = np.angle(xc[nc])
    # phase = 0.0
    print('phase:', phase)
    # frame *= np.exp(-1j * phase)

    ref_e = utils.calculate_signal_energy(x_preamble)
    rx_e = utils.calculate_signal_energy(frame[nc:nc + len(x_preamble)])
    agc_factor = np.sqrt(ref_e / rx_e)
    print('AGC values:', ref_e, rx_e, agc_factor)
    frame *= agc_factor
    sframe = frame[sample_nc:sample_nc + len(ref_frame)]
    # plt.plot(np.abs(ref_frame))
    # plt.plot(np.abs(frame))
    # plt.plot(np.abs(ac))
    # plt.plot(np.abs(xc))
    # plt.plot(cc)
    # # # plt.axvline(sample_nc, color='y')
    # print(np.abs(kframe - sframe) < 1e-4)
    # plt.plot(np.abs(kframe))
    # plt.plot(np.abs(sframe))
    # plt.show()
    return sframe


def estimate_frame_channel(H, fft_len, frame_len):
    used_sc = 52
    # subcarrier_map = mapping.get_subcarrier_map(fft_len, used_sc, dc_free=True)
    # subcarrier_map = np.roll(subcarrier_map, len(subcarrier_map) // 2)
    # ts = time.time()
    H[0] = (H[1] + H[-1]) / 2.
    # plt.plot(subcarrier_map, np.angle(H[subcarrier_map]))
    H = np.fft.fftshift(H)
    active_sc = np.arange((fft_len - used_sc)//2, (fft_len + used_sc)//2+1)
    active_sc = active_sc[3:-3]
    g = signal.gaussian(9, 1.0)
    g_factor = 1.
    g /= np.sqrt(g_factor * g.dot(g))
    Ha = H[active_sc]
    Hb = np.concatenate((np.repeat(Ha[0], 4), Ha, np.repeat(Ha[-1], 4)))
    Hg = np.correlate(Hb, g)

    # print('channel averaging error:', utils.calculate_signal_energy(Ha), utils.calculate_signal_energy(Hg), utils.calculate_signal_energy(Ha) / utils.calculate_signal_energy(Hg))
    Hg *= np.sqrt(utils.calculate_signal_energy(Ha) / utils.calculate_signal_energy(Hg))

    freqs = np.fft.fftfreq(len(H))
    freqs = np.fft.fftshift(freqs)

    fr_freqs = np.fft.fftfreq(frame_len)
    fr_freqs = np.fft.fftshift(fr_freqs)
    H_frame = np.interp(fr_freqs, freqs[active_sc], Hg.real) + 1j * np.interp(fr_freqs, freqs[active_sc], Hg.imag)
    # te = time.time()
    # print('interpolation timing:', te - ts)

    A = np.array([freqs[active_sc], np.ones(len(active_sc))])
    # print(np.shape(A))

    m, c = np.linalg.lstsq(A.T, H[active_sc])[0]
    Hl = m * freqs[active_sc] + c


    # plt.plot(freqs[active_sc], Ha.real)
    # plt.plot(freqs[active_sc], Ha.imag)
    # plt.plot(freqs[active_sc], Hg.real, marker='x', ms=10)
    # plt.plot(freqs[active_sc], Hl.real)

    # Hg = Ha

    fr_freqs = np.fft.fftfreq(frame_len)
    fr_freqs = np.fft.fftshift(fr_freqs)
    H_frame = np.interp(fr_freqs, freqs[active_sc], Hg.real) + 1j * np.interp(fr_freqs, freqs[active_sc], Hg.imag)
    # plt.plot(fr_freqs, H_frame.real)

    p_freqs = np.fft.fftfreq(fft_len * 9)
    p_freqs = np.fft.fftshift(p_freqs)
    H_p = np.interp(p_freqs, freqs[active_sc], Hg.real) + 1j * np.interp(p_freqs, freqs[active_sc], Hg.imag)
    # plt.plot(p_freqs, H_p.real)

    Hlf = m * fr_freqs + c
    # plt.plot(fr_freqs, Hlf.real)

    # plt.show()

    # plt.plot(xp, Hg.real, marker='x')
    # plt.plot(x, H_frame.real)
    # plt.plot(xf, H_p.real)

    # phase_m = m * fft_len / len(sframe)
    # p = np.arange(-frame_len, frame_len)
    # plt.plot(np.angle(H_p))
    # plt.plot(np.angle(H_frame))
    #
    # plt.plot(np.angle(Ha))
    # plt.plot(np.angle(Hg))
    # plt.plot(np.abs(Ha))
    # plt.plot(np.abs(Hg))
    # plt.show()
    return H_frame


def calculate_agc_factor(rx_preamble, x_preamble):
    ref_e = utils.calculate_signal_energy(x_preamble)
    rx_e = utils.calculate_signal_energy(rx_preamble)
    # print('AGC factor', np.sqrt(ref_e / rx_e))
    return np.sqrt(ref_e / rx_e)


def equalize_frame(sframe, x_preamble, fft_len, cp_len, cs_len):
    rx_preamble = sframe[cp_len:cp_len + len(x_preamble)]
    agc_factor = calculate_agc_factor(rx_preamble, x_preamble)
    # print('AGC values:', ref_e, rx_e, agc_factor)
    sframe *= agc_factor

    frame_start = cp_len + 2 * fft_len + cs_len + cp_len
    H, e0, e1 = preamble_estimate(rx_preamble, x_preamble, fft_len)
    H_estimate = estimate_frame_channel(H, fft_len, len(sframe))

    H_p = estimate_frame_channel(H, fft_len, fft_len * 9)
    p = sframe[frame_start:frame_start + fft_len * 9]

    P = np.fft.fft(p)
    P *= np.fft.fftshift(np.conj(H_p))
    p = np.fft.ifft(P)
    print('equalize p:', utils.calculate_signal_energy(p))

    F = np.fft.fft(sframe)
    F *= np.fft.fftshift(np.conj(H_estimate))
    sframe = np.fft.ifft(F)

    s = sframe[frame_start:frame_start + fft_len * 9]
    print('equalize s:', utils.calculate_signal_energy(s))
    # plt.plot(s.real)
    # plt.plot(p.real)
    # plt.plot(s.imag)
    # plt.plot(p.imag)
    # plt.plot(np.abs(P))
    # plt.plot(np.abs(F))
    # # plt.plot(np.abs(P - F))
    # plt.show()

    return sframe


def synchronize_freq_offsets(sframe, modulated_frame, x_preamble, fft_len, cp_len, samp_rate=12.5e6):
    rx_preamble = sframe[cp_len:cp_len + len(x_preamble)]
    frame_start = cp_len + 2 * fft_len + 16 + cp_len
    # rx_frame = sframe[frame_start:frame_start + len(modulated_frame)]
    # subcarrier_map = mapping.get_subcarrier_map(fft_len, 52, dc_free=True)
    H, e0, e1 = preamble_estimate(rx_preamble, x_preamble, fft_len)
    H_estimate = estimate_frame_channel(H, fft_len, len(sframe))

    H = np.fft.fftshift(H)

    used_sc = 52
    active_sc = np.concatenate((np.arange((fft_len - used_sc)//2, fft_len//2), np.arange(fft_len//2+1, (fft_len + used_sc)//2+1)))
    # print(active_sc)
    A = np.array([active_sc, np.ones(len(active_sc))])
    # print(np.shape(A))

    m, c = np.linalg.lstsq(A.T, np.unwrap(np.angle(H[active_sc])))[0]
    phase_m = m * fft_len / len(sframe)
    p = np.arange(-len(sframe)//2, len(sframe)//2)
    eq = np.exp(-1j * (p * phase_m + 0.0))
    plt.plot(np.angle(eq))
    plt.plot(np.angle(np.conj(H_estimate)))
    plt.show()

    eq = np.fft.fftshift(eq)
    F = np.fft.fft(sframe)
    # F *= eq
    F *= np.fft.fftshift(np.conj(H_estimate))
    sframe = np.fft.ifft(F)

    rx_frame = sframe[frame_start:frame_start + len(modulated_frame)]

    fine_cfo = m / (2 * np.pi)
    print('estimated fine freq offset: ', fine_cfo * samp_rate / fft_len)
    # phase_inc = sync.cfo_to_phase_increment(-fine_cfo, fft_len)
    # wave = sync.complex_sine(phase_inc, len(frame), 0.0)
    # frame *= wave
    #
    # sframe = frame[sample_nc:sample_nc + len(ref_frame)]
    # rx_preamble = sframe[cp_len:cp_len + len(x_preamble)]
    # H, e0, e1 = preamble_estimate(rx_preamble, x_preamble, fft_len)
    # H = np.fft.fftshift(H)
    # m, c = np.linalg.lstsq(A.T, np.unwrap(np.angle(H[active_sc])))[0]
    # fine_cfo = m / (2 * np.pi)
    # print('estimated fine freq offset: ', fine_cfo * samp_rate / fft_len)
    H_frame = np.fft.fft(rx_frame) / np.fft.fft(modulated_frame)
    H_frame = np.fft.fftshift(H_frame)[50:-50]
    # frame_sc = np.arange(-len(H_frame) // 2, len(H_frame) // 2)
    frame_sc = np.concatenate((np.arange(-len(H_frame) // 2, 0), np.arange(1, len(H_frame) // 2)))
    H_frame = H_frame[frame_sc]
    B = np.array([frame_sc, np.ones(len(frame_sc))])
    mf, cf = np.linalg.lstsq(B.T, np.unwrap(np.angle(H_frame)))[0]
    plt.plot(frame_sc, np.unwrap(np.angle(H_frame)))
    plt.plot(frame_sc, mf * frame_sc + cf, marker='x')
    print('m: ', m, mf, m * fft_len / len(rx_frame))
    print('c: ', c, cf)
    m *= fft_len / len(rx_frame)


    # e0 = np.fft.fftshift(e0)
    # e1 = np.fft.fftshift(e1)
    # plt.plot(np.angle(e0))
    # plt.plot(np.angle(e1))
    # plt.plot(np.angle(H))
    active_sc -= fft_len//2
    # plt.plot(active_sc, np.unwrap(np.angle(H[active_sc])), marker='x')
    plt.plot(active_sc, m * active_sc + c, marker='o')
    plt.grid()
    plt.show()
    return sframe


def synchronize_integrated(frame, ref_frame, x_preamble, fft_len, cp_len):
    samp_rate = 12.5e6
    ac = sync.auto_correlate_signal(frame, fft_len)

    nm = np.argmax(np.abs(ac[0:len(ac) // 2]))
    print('AC start: ', nm)
    # cfo = 2 * np.angle(ac[nm]) / (2. * np.pi)
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    print('CFO:', cfo, cfo * samp_rate / fft_len)

    phase_inc = sync.cfo_to_phase_increment(-cfo, fft_len)
    wave = sync.complex_sine(phase_inc, len(frame), 0.0)
    # frame *= wave

    ac = sync.auto_correlate_signal(frame, fft_len)
    cfo = np.angle(ac[nm]) / (2. * np.pi)
    print('CFO:', cfo, cfo * samp_rate / fft_len)
    ac = np.roll(ac, cp_len)

    xc = np.correlate(frame, x_preamble, 'valid')
    cc = sync.multiply_valid(np.abs(ac), np.abs(xc))
    nc = np.argmax(np.abs(cc[0:len(cc)//2]))
    print('correlation frame start:', nc)
    sample_nc = nc - cp_len
    print('sample frame start:     ', sample_nc)
    p_len = cp_len + 2 * fft_len + cp_len // 2 + cp_len
    print('data frame start:       ', sample_nc + p_len)
    phase = np.angle(xc[nc])
    # phase = 0.0
    print('phase:', phase)
    # frame *= np.exp(-1j * phase)

    ref_e = utils.calculate_signal_energy(x_preamble)
    rx_e = utils.calculate_signal_energy(frame[nc:nc + len(x_preamble)])
    agc_factor = np.sqrt(ref_e / rx_e)
    print('AGC values:', ref_e, rx_e, agc_factor)
    frame *= agc_factor
    sframe = frame[sample_nc:sample_nc + len(ref_frame)]
    # plt.plot(np.abs(ref_frame))
    plt.plot(np.abs(frame))
    plt.plot(np.abs(ac))
    plt.plot(np.abs(xc))
    plt.plot(cc)
    # # plt.axvline(sample_nc, color='y')
    plt.show()
    # return

    return sframe


def preamble_estimate(rx_preamble, x_preamble, fft_len):
    # e = np.fft.fft(rx_preamble) / np.fft.fft(x_preamble)
    # e = np.concatenate((e[0:fft_len//2], e[-fft_len//2:]))
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


def demodulate_data_frame(frame, rx_kernel, demapper, n_data_syms):
    ref_syms64 = frame.astype(dtype=np.complex64)
    txd = rx_kernel.demodulate(ref_syms64)
    txd = demapper.demap_from_resources(txd, n_data_syms)
    return txd


def demodulate_equalize_frame(frame, rx_kernel, demapper, H, fft_len, n_data_syms):
    timeslots = len(frame) // fft_len
    syms64 = frame.astype(dtype=np.complex64)

    fd_syms = rx_kernel.fft_filter_downsample(syms64)
    fd_syms = np.reshape(fd_syms, (fft_len, -1))
    fd_eq = np.zeros(np.shape(fd_syms), dtype=np.complex64)
    for i in range(fft_len):
        fd_eq[i] = fd_syms[i] * np.conj(H[i])
    fd_eq_vec = fd_eq.flatten()
    print(fd_eq_vec[0:timeslots] == fd_eq[0])

    r_data = rx_kernel.transform_subcarriers_to_td(fd_eq_vec)
    txd = demapper.demap_from_resources(r_data, n_data_syms)
    return txd


def plot_constellation(ref_data, rx_data, rx_eq_data, start, end):
    r = ref_data[start:end]
    x = rx_data[start:end]
    e = rx_eq_data[start:end]

    for i in range(len(r)):
        plt.plot([r[i].real, x[i].real], [r[i].imag, x[i].imag], color='g')
        # plt.arrow(r[i].real, r[i].imag, x[i].real - r[i].real, x[i].imag - r[i].imag) #, head_width=0.05, head_length=0.1, fc='k', ec='k')

    plt.scatter(r.real, r.imag, color='r', label='reference')
    plt.scatter(x.real, x.imag, color='b', label='RX')
    plt.scatter(e.real, e.imag, color='g', label='RX EQ')
    consti = np.array([1+1j, 1-1j, -1+1j, -1-1j, ])
    consti /= np.sqrt(2.)
    plt.scatter(consti.real, consti.imag, color='m', marker='x')
    plt.grid()
    my_lims = 1.5
    lim_span = [-my_lims, my_lims]
    plt.xlim(lim_span)
    plt.ylim(lim_span)
    plt.legend()


def calculate_frame_ber(ref_symbols, rx_symbols):
    b = utils.demodulate_qpsk(ref_symbols)
    # print(b[0:10])
    br = utils.demodulate_qpsk(rx_symbols)
    # print(br[0:10])
    return (len(b) - np.sum(b == br)) / len(b)


def gr_load_frame():
    # gfdm_file_sync.main()

    tb = gfdm_file_sync.gfdm_file_sync()
    tb.start()
    import time
    time.sleep(1)
    # while not len(tb.sync_sink.data()):
    #     print('wait')
    tb.stop()
    tb.wait()
    rx_frames = np.array(tb.sync_sink.data())
    bursts = np.array(tb.burst_sink.data())
    # plt.plot(np.abs(rx_frames))
    # plt.show()

    tags = tb.sync_sink.tags()
    rx_len = len(rx_frames)
    print(rx_len, rx_len / 3200)
    print(len(tags))
    return rx_frames


def calculate_avg_phase(rx_data, ref_data):
    A = np.array([np.arange(len(rx_data)), np.ones(len(rx_data))])

    phase_err = np.angle(rx_data) - np.angle(ref_data)
    phase_err = np.unwrap(phase_err)
    # plt.plot(phase_err)
    # plt.plot(np.unwrap(phase_err))

    m, c = np.linalg.lstsq(A.T, phase_err)[0]
    # plt.plot(np.arange(len(phase_err)), m * np.arange(len(phase_err)))
    # pm = np.reshape(phases, (-1, active_subcarriers))
    # for i in range(active_subcarriers):
    #     p = pm[:, i]
    #     plt.plot(p)

    avg_phase = np.sum(phase_err) / len(phase_err)
    print('AVG phase shift: ', avg_phase)
    print('return val:      ', c + m * len(phase_err))
    # print('lin reg: ', m, c)
    # plt.plot(phase_err - avg_phase)
    # plt.title('phase error')
    # plt.show()
    # return avg_phase
    return c + m * len(phase_err)


def demodulate_frame(rx_data_frame, modulated_frame, rx_kernel, demapper, data, timeslots, fft_len, H_estimate=None):
    ref_data = demodulate_data_frame(modulated_frame, rx_kernel, demapper, len(data))
    rx_data = demodulate_data_frame(rx_data_frame, rx_kernel, demapper, len(data))

    if H_estimate is None:
        d = np.array([1+1j, -1+1j, -1-1j, 1-1j]) / np.sqrt(2.)
    else:
        d = rx_kernel.demodulate_equalize(rx_data_frame.astype(dtype=np.complex64), H_estimate.astype(dtype=np.complex64))
        d = demapper.demap_from_resources(d.astype(dtype=np.complex64), len(data))

    mse = np.sum(np.abs(rx_data - ref_data) ** 2) / len(ref_data)
    print('frame EQ demodulated energy', utils.calculate_signal_energy(rx_data), utils.calculate_signal_energy(rx_data) / len(rx_data), utils.calculate_signal_energy(ref_data) / len(ref_data))

    # calculate_avg_phase(rx_data, ref_data)

    fber = calculate_frame_ber(ref_data, rx_data)
    print('Frame BER: ', fber, 'with MSE: ', mse)

    plot_constellation(ref_data, rx_data, d, 0, timeslots * fft_len)
    plt.show()


def rx_oversampled(frames, ref_frame, modulated_frame, x_preamble, data, rx_kernel, demapper, timeslots, fft_len, cp_len, cs_len):
    ref_frame_os = signal.resample(ref_frame, 2 * len(ref_frame))
    x_preamble_os = signal.resample(x_preamble, 2 * len(x_preamble))

    nyquist_frame_len = cp_len + 2 * fft_len + cs_len + cp_len + timeslots * fft_len + cs_len
    n_frames = np.shape(frames)[0]
    sync_frames = np.zeros((n_frames, nyquist_frame_len), dtype=np.complex)
    print('nyquist sampled frame len', nyquist_frame_len, 'with n_frames', n_frames)
    f_start = cp_len + 2 * fft_len + cs_len
    d_start = f_start + cp_len
    print('data start: ', d_start)
    for i, f in enumerate(frames[0:2]):
        tf = np.roll(f, 1)
        tf[0] = 0
        ff = signal.resample(tf, len(f) // 2)
        sframe = synchronize_time(ff, ref_frame_os, x_preamble_os, 2 * fft_len, 2 * cp_len)
        sframe = signal.resample(sframe, len(sframe) // 2)
        sframe = synchronize_freq_offsets(sframe, modulated_frame, x_preamble, fft_len, cp_len, samp_rate=3.125e6)
        print(len(sframe), len(ref_frame))
        rx_preamble = sframe[cp_len:cp_len + 2 * fft_len]
        avg_phase = calculate_avg_phase(rx_preamble, x_preamble)
        # m, c = calculate_avg_phase(rx_preamble, x_preamble)
        # avg_phase = calculate_avg_phase(sframe, ref_frame)
        # phase_eqs = m * np.arange(-cp_len, len(sframe) - cp_len) + c
        # sframe *= np.exp(-1j * phase_eqs)
        # sframe *= np.exp(-1j * avg_phase)
        sync_frames[i] = sframe
        rx_data_frame = sframe[d_start:d_start + fft_len * timeslots]
        # # rx_data_frame *= np.exp(-1j * avg_phase)
        #
        demodulate_frame(rx_data_frame, modulated_frame, rx_kernel, demapper, data, timeslots, fft_len)

    for i, f in enumerate(sync_frames[0:3]):
        rx_data_frame = f[d_start:d_start + fft_len * timeslots]
        demodulate_frame(rx_data_frame, modulated_frame, rx_kernel, demapper, data, timeslots, fft_len)


def rx_nyquist_sampled(frames, ref_frame, modulated_frame, x_preamble, data, rx_kernel, demapper, timeslots, fft_len, cp_len, cs_len):
    f_start = cp_len + 2 * fft_len + cs_len
    d_start = f_start + cp_len
    print('data start: ', d_start)
    for f in frames[0:3]:
        for offset in range(3):
            print('\n\n OFFSET: ', offset)
            ff = f[offset:-(4-offset)]
            print('offset frame len', len(ff))
            df = signal.resample(ff, len(ff) // 4)
            sframe = synchronize_time(df, ref_frame, x_preamble, fft_len, cp_len, samp_rate=3.125e6)
            sframe = synchronize_freq_offsets(sframe, modulated_frame, x_preamble, fft_len, cp_len, samp_rate=3.125e6)
            # print(len(sframe), len(ref_frame))
            rx_preamble = sframe[cp_len:cp_len + 2 * fft_len]
            avg_phase = calculate_avg_phase(rx_preamble, x_preamble)
            # m, c = calculate_avg_phase(rx_preamble, x_preamble)
            # avg_phase = calculate_avg_phase(sframe, ref_frame)
            # phase_eqs = m * np.arange(-cp_len, len(sframe) - cp_len) + c
            # sframe *= np.exp(-1j * phase_eqs)
            sframe *= np.exp(-1j * avg_phase)
            rx_data_frame = sframe[d_start:d_start + fft_len * timeslots]

            demodulate_frame(rx_data_frame, modulated_frame, rx_kernel, demapper, data, timeslots, fft_len)

            # plt.plot(np.abs(sframe))
            # plt.plot(np.abs(ref_frame))
            # plt.plot(np.abs(sframe - ref_frame))
            # plt.show()


def rx_demodulate(frames, ref_frame, modulated_frame, x_preamble, data, rx_kernel, demapper, timeslots, fft_len, cp_len, cs_len):
    f_start = cp_len + 2 * fft_len + cs_len
    d_start = f_start + cp_len
    print('data start: ', d_start)
    sync_kernel = cgfdm.py_auto_cross_corr_multicarrier_sync_cc(64, 32, x_preamble)

    estimator = validation_utils.frame_estimator(x_preamble, fft_len, timeslots, 52)
    c_est = cgfdm.py_preamble_channel_estimator_cc(9, 64, 52, True, x_preamble)
    p_filter_taps = c_est.preamble_filter_taps()
    print(estimator._p_filter)
    print(p_filter_taps)
    print(np.abs(p_filter_taps - estimator._p_filter) < 1e-6)

    for f in frames[0:30]:
        rxs = time.time()
        rx_preamble = f[cp_len:cp_len + 2 * fft_len]
        agc_factor = sync_kernel.calculate_preamble_attenuation(rx_preamble.astype(dtype=np.complex64))
        f = sync_kernel.normalize_power_level(f, agc_factor)

        rx_preamble = f[cp_len:cp_len + 2 * fft_len]
        H_estimate = estimator.estimate_frame(rx_preamble)

        rx_t_frame = f[d_start:d_start + fft_len * timeslots]

        kde = rx_kernel.demodulate_equalize(rx_t_frame.astype(dtype=np.complex64), H_estimate.astype(dtype=np.complex64))
        de = demapper.demap_from_resources(kde.astype(dtype=np.complex64), len(data))
        rxe = time.time()
        print('receiver chain time: ', 1e6 * (rxe - rxs), 'us')

        # p_active = np.concatenate((np.arange(1, 27), np.arange(37, 64)))
        # cp_est = c_est.estimate_preamble_channel(rx_preamble)
        # p_est = estimator._estimate_preamble(rx_preamble)
        # print('Python vs C Preamble channel correct:', np.all(np.abs(cp_est[p_active] - p_est[p_active]) < 1e-4), np.all(np.angle(cp_est[p_active]) - np.angle(p_est[p_active]) < 1e-6))

        # pf = estimator._filter_preamble_estimate(p_est)
        # cf = c_est.filter_preamble_estimate(p_est.astype(dtype=np.complex64))
        # print(np.all(np.angle(pf) - np.angle(cf) < 1e-6))

        # frame_estimate = c_est.interpolate_frame(pf.astype(dtype=np.complex64))
        rxs = time.time()
        frame_estimate = c_est.estimate_frame(rx_preamble)
        rxe = time.time()
        print('CPP channel estimator duration: ', 1e6 * (rxe - rxs), 'us')
        print('Python vs C Preamble channel correct:', np.all(np.abs(frame_estimate - H_estimate) < 1e-6))
        H_estimate = frame_estimate
        # plt.plot(H_estimate.real)
        # plt.plot(H_estimate.imag)
        # plt.plot(frame_estimate.real)
        # plt.plot(frame_estimate.imag)
        #
        # # plt.plot(pf.real)
        # # plt.plot(pf.imag)
        # # plt.plot(cf.real)
        # # plt.plot(cf.imag)
        # plt.show()


        P = np.fft.fft(rx_t_frame)
        P *= np.conj(H_estimate)
        rx_t_frame = np.fft.ifft(P)
        # kd = rx_kernel.demodulate(rx_t_frame.astype(dtype=np.complex64))
        # d = demapper.demap_from_resources(kd.astype(dtype=np.complex64), len(data))

        print('kernel FER: ', calculate_frame_ber(data, de))

        sframe = equalize_frame(f, x_preamble, fft_len, cp_len, cs_len)

        # rx_preamble = sframe[cp_len:cp_len + 2 * fft_len]
        # avg_phase = calculate_avg_phase(rx_preamble, x_preamble)
        # sframe *= np.exp(-1j * avg_phase)
        rx_data_frame = sframe[d_start:d_start + fft_len * timeslots]

        ekd = utils.calculate_signal_energy(rx_t_frame - rx_data_frame)
        print('MSE for kernel vs Python demod', ekd, ekd / len(rx_data_frame))
        plt.show()
        # plt.scatter(d.real, d.imag, label='kernel')
        plt.scatter(de.real, de.imag, label='EQ kern')
        H_estimate = np.ones(len(H_estimate))

        # demodulate_frame(rx_data_frame, modulated_frame, rx_kernel, demapper, data, timeslots, fft_len, H_estimate)
        demodulate_frame(rx_t_frame, modulated_frame, rx_kernel, demapper, data, timeslots, fft_len, H_estimate)



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
    filename = '/lhome/records/gfdm_replay_ref_frame_time_synced.dat'
    # filename = '/lhome/records/gfdm_gr_fg_synced_frames.dat'
    # filename = '/lhome/records/gfdm_ref_frame_50ms_slice.dat'
    slice_len = 800
    offset = 0
    n_frames = 200
    frame_start = 0
    frame_end = 800
    frame = converter.load_gr_iq_file(filename)[offset:]
    n_max_frames = int(len(frame) // slice_len)
    print('max number of frames:', n_max_frames)
    frame = frame[2000 * slice_len:2000 * slice_len + slice_len * n_frames]
    frames = np.reshape(frame, (-1, slice_len))
    # frames = frames[:, frame_start:frame_end]
    # frame = converter.load_gr_iq_file(filename)
    print('num samples', len(frame))
    # f_frame = np.fft.fft(frame)
    # # plt.semilogy(np.abs(f_frame))
    # plt.plot(np.abs(frame))
    # plt.show()
    # for f in frames:
    #     plt.plot(np.abs(f))
    # # # # # # plt.semilogy(*signal.welch(frame))
    # plt.show()
    # return


    ref_frame, modulated_frame, x_preamble, data, freq_filter_taps = validation_utils.generate_reference_frame(timeslots, fft_len, active_subcarriers, cp_len, cs_len, alpha)
    #ref_frame, modulated_frame, x_preamble, data, freq_filter_taps = validation_utils.generate_sc_qpsk_frame(timeslots, fft_len, active_subcarriers, cp_len, cs_len, alpha)

    rx_kernel = cgfdm.py_receiver_kernel_cc(timeslots, fft_len, 2, freq_filter_taps)
    demapper = cgfdm.py_resource_demapper_kernel_cc(timeslots, fft_len, active_subcarriers, subcarrier_map, True)

    # rx_oversampled(frames, ref_frame, modulated_frame, x_preamble, data, rx_kernel, demapper, timeslots, fft_len, cp_len, cs_len)
    # rx_nyquist_sampled(frames, ref_frame, modulated_frame, x_preamble, data, rx_kernel, demapper, timeslots, fft_len, cp_len, cs_len)
    rx_demodulate(frames, ref_frame, modulated_frame, x_preamble, data, rx_kernel, demapper, timeslots, fft_len, cp_len, cs_len)
    return


if __name__ == '__main__':
    main()
