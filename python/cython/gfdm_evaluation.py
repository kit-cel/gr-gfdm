#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
from future_builtins import *
import numpy as np
import scipy.signal as signal
import cgfdm as gfdm
import matplotlib.pyplot as plt


def map_to_waveform_resources(syms, subcarrier_map, active_subcarriers, fft_len):
    n_data_syms = len(syms)
    ts = int(np.ceil(n_data_syms / active_subcarriers))
    n_frame_syms = active_subcarriers * ts
    # print(n_data_syms, ts)
    s = np.concatenate((syms, np.zeros(n_frame_syms - n_data_syms)))
    # print(syms)
    s = np.reshape(s, [-1, active_subcarriers]).T
    # print(s)
    nt = np.shape(s)
    # print(nt)
    frame = np.zeros((fft_len, ts), dtype=np.complex64)
    frame[subcarrier_map, :] = s
    print('frame shape:', np.shape(frame))
    return frame.flatten()


def get_random_qpsk(nsamples):
    d = np.random.randint(0, 2, 2 * nsamples) * -2. + 1.
    d = np.reshape(d, (2, -1))
    d = d.astype(dtype=np.complex)
    return d[0] + 1j * d[1]


def plot_waterfall_db(syms, fft_len, overlap):
    img = prepare_waterfall(syms, fft_len, overlap)
    img = np.abs(img)
    # img = np.log10(img)
    # print(np.max(np.abs(img)))
    plot_2d_grid(img)


def plot_2d_grid(syms):
    plt.imshow(np.abs(syms), cmap=plt.get_cmap('gnuplot2'))
    plt.show()


def prepare_waterfall(syms, fft_len, overlap):
    img = np.reshape(syms, (-1, fft_len)).T
    print(np.shape(img))
    img = np.fft.fft(img, n=fft_len, axis=0)
    return img.T


def plot_waterfall(symbols, fft_len, overlap=-1, cmap=plt.get_cmap('gnuplot2')):
    if overlap < 0:
        overlap = fft_len / 4
    Pxx, freqs, bins, im = plt.specgram(symbols, NFFT=fft_len, Fs=fft_len, noverlap=overlap)

    Pxx = np.fft.fftshift(Pxx, axes=0)
    plt.imshow(10 * np.log10(Pxx.T), cmap=cmap)
    plt.xlabel('FFT bins')
    plt.ylabel('time slots')
    # plt.show()


def main():
    my_set = {'active_subcarriers': 110, 'frame_symbols': 4096, 'fft_len': 128, 'frame_time_symbols': 38, 'cp_len': 6}
    n_subcarriers = my_set['fft_len']
    n_timeslots = my_set['frame_time_symbols']
    overlap = 2
    taps = np.array([1, ] * n_timeslots * overlap, dtype=np.complex)
    kernel = gfdm.py_modulator_kernel_cc(n_timeslots, n_subcarriers, overlap, taps)
    print(kernel.block_size())
    syms = map_to_waveform_resources(get_random_qpsk(my_set['frame_symbols']), np.arange(my_set['active_subcarriers']) + 9, my_set['active_subcarriers'], my_set['fft_len'])

    frame = kernel.modulate(syms)

    plot_waterfall(frame, 128)
    plt.show()





if __name__ == '__main__':
    main()
