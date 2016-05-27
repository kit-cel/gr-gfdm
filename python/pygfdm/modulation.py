#!/usr/bin/env python2.7

import numpy as np
import commpy as cp
import matplotlib.pyplot as plt
# from synchronization import *
from gfdm_plot_utils import plot_gfdm_matrix
from filters import gfdm_filter_taps, gfdm_freq_taps, gfdm_freq_taps_sparse


__all__ = [
    'transmitMatrix',
    'fourierMatrix',
    'samplingMatrix',
    'randomQAMSymbols',
    'gfdm_tx',
    'gfdm_rx']


def gfdm_modulation_matrix(filter_taps, M, K, oversampling_factor=1, group_by_subcarrier=False):
    '''
    This function returns a GFDM modulation matrix
    :param filter_taps: M*K*N length filter tap array
    :param M: number of time slots
    :param K: number of subcarriers
    :param oversampling_factor: factor for oversampling
    :param group_by_subcarrier: if True group by time symbol not subcarrier
    :return: modulation matrix
    [0] Generalized Frequency Division Multiplexing for 5th Generation Cellular Networks
    [1] Generalized frequency division multiplexing: Analysis of an alternative multi-carrier technique for next generation cellular systems

    CAREFUL: It is tricky for a number of reasons.
    First off, only in [1] oversampling is considered. Mostly it's undesirable.
    Secondly, the definitions differ slightly. In [0] frequency modulation is defined by -1j * 2 * np.pi,
    whereas in [1] it is 1j * 2 * np.pi. Notice the sign change!
    '''
    N = M * K

    filter_taps = np.roll(filter_taps, (N * oversampling_factor) // 2)
    A = np.zeros((N * oversampling_factor, N), dtype=np.complex)

    n = np.arange(N * oversampling_factor, dtype=np.complex)
    for m in range(M):
        for k in range(K):
            f_mod = np.exp(1j * 2 * np.pi * (float(k) / (K * oversampling_factor)) * n)
            g = filter_taps * f_mod
            g = np.roll(g, m * K * oversampling_factor)
            A[:, m * K + k] = g

    if group_by_subcarrier:
        indices = np.arange(M * K)
        indices = np.reshape(indices, (-1, K)).T.flatten()
        A = A[:, indices]
    return A


def transmitMatrix(filtertype, alpha, M, K, N):
    # replaces old definition.
    taps = gfdm_filter_taps(filtertype, alpha, M, K, N)
    return gfdm_modulation_matrix(taps, M, K, N, False)


def randomQAMSymbols(length, M):
    '''
     length: number of symbols to generate
     M: M-QAM - Order (4,16,64,...)
    '''
    n = np.sqrt(M / 4)
    if np.around(n) - n > 1e-10:
        raise Exception('M must be power of 4')
    n = int(n)
    n_M_pos = np.array([1 + 2 * i for i in xrange(n)])
    n_M_neg = np.array([-1 - 2 * i for i in xrange(n)])
    choices = np.concatenate((n_M_pos, n_M_neg))
    return np.array(
        [np.random.choice(choices) + 1j * np.random.choice
        (choices) for i in xrange(length)])


def gfdm_tx(x, filtertype, alpha, M, K, L, N):
    '''
    x: Input-Symbols (length M*K)
    filtertype: ['rrc','rc']
    alpha: rolloff-factor
    M: number of timeslots
    K: number of subcarrier
    oversampling_factor: sometimes referred to as N
    '''
    A = transmitMatrix(filtertype, alpha, M, K, N)
    A *= M
    tx = A.dot(x)
    return tx


def gfdm_rx(y, filtertype, alpha, M, K, L, N, QAM, J):
    '''
    y: Transmit-Symbols (length M*K*N)
    filtertype: ['rrc','rc']
    alpha: rolloff-factor
    rx_strat: ['zf','mf']
    M: number of timeslots
    K: numbor of subcarrier
    N: oversampling-factor
    '''
    A = transmitMatrix(filtertype, alpha, M, K, N)
    A_rx = A.conj().transpose()
    rx = A_rx.dot(y)
    return rx


def gfdm_tx_fft2(x, filtertype, alpha, M, K, L, N):
    '''
    x: Input-Array (length: M*K symbols)
    filtertype: ('rrc'|'rc')
    alpha: (0,1) float
    M: number of slots
    K: number of subcarriers
    L: freq-domain length of filter

    Low-complexity transmitter implementation as proposed by G. Fettweis
    '''
    h = gfdm_filter_taps(filtertype, alpha, M, K, 1)
    H = gfdm_freq_taps(h)
    H_sparse = gfdm_freq_taps_sparse(H, M, L)

    # Sort Input subcarrierwise
    x = reshape_input(x, M, K)
    x_out = np.zeros((M * K) + (L - 1) * M, dtype='complex')
    for k in xrange(K):
        # M rows and L columns with respective FFT output
        # pick symbols per subcarrier
        x_k = x[k * M:((k + 1) * M)]
        # perform fft and switch to frequency domain
        x_f = np.fft.fft(x_k)
        # copy values of M-point DFT to obtain MK-point DFT
        x_f_L = np.tile(x_f, L)
        # make data-vector 'sparse'
        # x_f_L = np.concatenate((x_f_K[0:(M*L)/2], x_f_K[-(M*L)/2:]))
        # filter with sparse filter taps in frequency domain
        x_fil = np.multiply(x_f_L, H_sparse)
        # Add data-vector to correct position -max neg frequency : 0 :
        # max_pos_frequency
        x_out[k * M:(k + L) * M] = x_out[k * M:(L + k) * M] + np.fft.fftshift(x_fil)
    # Add 'oversampled' parts of first subcarrier to end and 'oversampled' parts
    # of last subcarrier to start
    x_first = x_out[0:(L - 1) * M / 2]
    x_last = x_out[-(L - 1) * M / 2:]
    x_out = x_out[(L - 1) * M / 2:-(L - 1) * M / 2]
    x_out[0:(L - 1) * M / 2] = x_out[0:(L - 1) * M / 2] + x_last
    x_out[-(L - 1) * M / 2:] = x_out[-(L - 1) * M / 2:] + x_first

    x_t = np.fft.ifft(np.fft.ifftshift(x_out))
    x_t *= 1.0 / K
    return x_t


def gfdm_rx_fft2(y, filtertype, alpha, M, K, L, N, QAM, J):
    '''
    y: transmitted gfdm-block (length: M*K samples)
    filtertype: ('rrc'|'rc')
    alpha: (0,1) float
    M: number of slots
    K: number of subcarriers
    L: freq-domain length of filter
    Low-complexity receiver implementation as proposed by G.Fettweis
    (based on sparse frequency Domain Processing)

    output: demodulated samples in original order (first K samples in timeslot 1, second K ...)
    '''
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M * K, alpha, K, 1)
    h = np.roll(h, h.shape[-1] / 2)
    H_rx = np.fft.fft(h)
    H_sparse = np.concatenate((H_rx[0:M * L / 2], H_rx[-M * L / 2:]))
    y_ifft = np.array([])
    y = (1.0 / K) * y
    # Transfer input to frequency domain and center around 0th frequency bin
    y_f = np.fft.fftshift(np.fft.fft(y))
    # Filter and superposition in frequency domain
    Y_fs = gfdm_rx_filsup(y_f, H_sparse, M, K, L)
    # Demodulate per subcarrier
    y_ifft = gfdm_rx_demod(Y_fs, K)
    if J > 0:
        y_ifft = gfdm_rx_sic(K, M, J, H_sparse, y_ifft, Y_fs, QAM)
        y_ifft = np.reshape(y_ifft, (K * M))
    # Sort output in timeslot,subcarrier order
    y_ifft = reshape_input(y_ifft, K, M)
    return y_ifft


def gfdm_rx_demod(Y_fs, K):
    '''
    Y_fs: received samples filtered and superpositioned in frequency domain (not centered) KxM-array
    K: Number of subcarriers

    output: demodulated samples in subcarrier order (first M samples are on subcarrier 1, second M....)
    '''
    y_ifft = np.array([])
    for k in xrange(K):
        y_ifft = np.concatenate((y_ifft, np.fft.ifft(Y_fs[k])))
    return y_ifft


def gfdm_rx_filsup(y_f, H_sparse, M, K, L):
    '''
    y_f: input samples centered in frequency domain 1xK*M-array
    H_sparse: Rx-filter per subcarrier - length (M*L)
    M: number of time slots
    K: number of subcarrier
    L: width of sparse Rx-filter in number of subcarrier

    output: (K,M) - array
    '''
    y_out = np.empty((K, M), dtype='complex')
    y_f = np.concatenate((y_f[-(L - 1) * M / 2:], y_f, y_f[0:(L - 1) * M / 2]))
    for k in xrange(K):
        # select kth subcarrier
        y_down = y_f[k * M:(k + L) * M]
        # 'uncenter' in Frequency domain
        y_down = np.fft.ifftshift(y_down)
        # apply filter in frequency domain (not centered)
        y_filter = np.multiply(y_down, H_sparse)
        # Superposition L samples in frequency domain
        y_out[k] = np.sum(y_filter.reshape(L, M), axis=0)
    return y_out


def gfdm_rx_sic(K, M, J, H_sparse, d_rx, Y_fs, QAM):
    '''
    K: Number of subcarriers
    M: Number of slots
    J: Number of Iterations for Interference Cancellation
    H_sparse: Rx-filter of length M*L in frequency domain
    d_rx: mapped symbols before Interference cancellation (sorted by subcarrier)
    Y_fs: filtered, superpositioned input samples in frequency domain (not centered) KxM-array
    QAM: QAM order
    '''
    # Receive all subcarriers in F-Domain
    # map each symbol to closest QAM - Point
    # d_rx s
    qam_mod = cp.QAMModem(QAM)
    # Calculate rising/falling flank interference coefficients
    H_rf = np.multiply((H_sparse[0:M] / K), (H_sparse[M:] / K))
    # Reshape mapped symbols into per-subcarrier array
    d_p = np.empty((K, M), dtype='complex')
    # d_p (K,M)
    for k in xrange(K):
        d_p[k] = qam_mod.mapping(d_rx[k * M:(k + 1) * M], 'hard')
    for j in xrange(J):
        y = np.empty((K, M), dtype='complex')
        for k in xrange(K):
            y[k] = Y_fs[k] - H_rf * np.fft.fft(d_p[(k - 1) % K] + d_p[(k + 1) % K])
            # Recalculate d_rx
        d_rx = gfdm_rx_demod(y, K)
        for k in xrange(K):
            d_p[k] = d_rx[k * M:(k + 1) * M]
            d_p[k] = qam_mod.mapping(d_p[k], 'hard')
    return d_rx


def add_cp(x, n):
    return np.append(x[-n:], x)


# [2] "Bit Error Rate Performance of Generalized Frequency Division Multiplexing"
def get_data_matrix(data, K, group_by_subcarrier=False):
    # function yields data matrix according to [2]
    if group_by_subcarrier:
        # alternative grouping. Used in other papers.
        return np.reshape(data, (-1, K))
    else:
        # data grouped as described in [2]
        return np.reshape(data, (K, -1)).T


def reshape_input(x, M, K, group_by_subcarrier=True):
    '''
    1. pick every m*Kth symbol and append to output.
    2. Increase counter one time
    3. perform step 1.+2. M times
    '''
    D = get_data_matrix(x, K, group_by_subcarrier=group_by_subcarrier)
    return D.T.flatten()


def main():
    print 'This is main: currently nothing to do here.'


if __name__ == '__main__':
    main()
