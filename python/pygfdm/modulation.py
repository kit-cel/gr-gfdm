#!/usr/bin/env python2.7

import numpy as np
import commpy as cp
import matplotlib.pyplot as plt
# from synchronization import *
from gfdm_plot_utils import plot_gfdm_matrix


__all__ = [
    'transmitMatrix',
    'fourierMatrix',
    'samplingMatrix',
    'randomQAMSymbols',
    'gfdm_tx',
    'gfdm_rx']


# FIXME TransmitMatrix should group different subcarriers on timeslot-basis


def gfdm_filter_taps(filtertype, alpha, M, K, oversampling_factor):
    N = oversampling_factor
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M * K * N, alpha, N * K, 1)
    elif filtertype == "rc":
        time_h, h = cp.rcosfilter(M * K * N, alpha, N * K, 1)
    return h


def gfdm_freq_taps(h):
    h = np.roll(h, h.shape[-1] / 2)
    # print h == np.fft.fftshift(h)
    H = np.fft.fft(h)
    return H


def gfdm_freq_taps_sparse(H, M, L):
    H_sparse = np.concatenate((H[0:(M * L) / 2], H[-(M * L) / 2:]))
    return H_sparse


def transmitMatrix_core(filter_taps, M, K, N):
    '''
        FIXME: Apparently a BUG is killing this function.
        Create GFDM modulation matrix

        filter_taps: prototype filter taps. centered
        M : number of symbol time slots
        K : number of subcarriers
        N: oversampling factor
    '''
    print M, K, N
    # Move filter cyclic
    G_tx = np.array([np.roll(filter_taps, m - (M * K * N / 2)) for m in xrange(M * K * N)])
    S_mn = samplingMatrix(M, K * N)
    S_nm = samplingMatrix(K, M)
    if N > 1:
        # if oversampling is specified add zeros to samplingMatrix
        S_nm = np.insert(
            S_nm, M * K / 2, np.zeros((M * K * (N - 1), K), dtype='complex'),
            axis=0)
    W_H = fourierMatrix(M * K * N).conj().transpose()
    # Resample Filter
    G_tx_s = np.dot(G_tx, S_mn)
    # Resample FourierMatrix
    W_s = np.dot(S_nm.transpose(), W_H)
    # compute and use all elements of the main diagonal W_s.dot(G_tx_s)
    A = np.array([(np.kron(W_s.transpose()[n], G_tx_s[n]))
                  for n in xrange(K * M * N)])

    return A


def gfdm_modulation_matrix(filter_taps, M, K, oversampling_factor=1, rearrange_indices=True):
    '''
    This function returns a GFDM modulation matrix
    :param filter_taps: M*K*N length filter tap array
    :param M: number of time slots
    :param K: number of subcarriers
    :param oversampling_factor: factor for oversampling
    :param rearrange_indices: if True group by time symbol not subcarrier
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
        g = np.roll(filter_taps, m * K * oversampling_factor)
        for k in range(K):
            f_mod = np.exp(1j * 2 * np.pi * (float(k) / (K * oversampling_factor)) * n)
            g = g * f_mod
            A[:, m * K + k] = g

    if rearrange_indices:
        indices = np.arange(M * K)
        indices = np.reshape(indices, (-1, K)).T.flatten()
        A = A[:, indices]
    return A


def transmitMatrix(filtertype, alpha, M, K, N):
    '''
        Create Convolution Matrix for pulse shaping

        filtertype : (rrc,rc)
        alpha : roll-off-factor
        sampling_rate : sampling rate (in Hz)
        symbol_period : symbol period (in s)
        M : number of symbol time slots
        K : number of subcarriers

        h_matrix: array of impulse responses for time slot (0...M-1)
    '''
    h = gfdm_filter_taps(filtertype, alpha, M, K, N)
    return transmitMatrix_core(h, M, K, N)


def fourierMatrix(N):
    i, j = np.meshgrid(np.arange(N), np.arange(N))
    omega = np.exp(- 2 * np.pi * 1j / N)
    W = np.power(omega, i * j) / np.sqrt(N)
    return W


def samplingMatrix(M, K):
    output = np.zeros((M * K, M), dtype='complex')
    for n in xrange(M * K):
        for m in xrange(M):
            if n == (m * K):
                output[n][m] = 1
    return output


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
    N: oversampling-factor
    '''
    A = transmitMatrix(filtertype, alpha, M, K, N)
    A = A * M
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
    # if rx_strat == "zf":
    #    A_rx = np.linalg.pinv(A)
    # else:
    # A_rx = np.linalg.pinv(A)/M
    A_rx = A.conj().transpose()
    rx = np.array([])
    rx = A_rx.dot(y)
    return rx


def gfdm_tx_fft(x, filtertype, alpha, M, K):
    '''
        Realization of GFDM-Transmitter in FFT:
            Required input: x a np.array of length M*K
            FFT is applied in shifted version (zero-frequency term is centered)
            First symbol is on -freq_max and last symbol ist on freq_max
            h: Prototype-filter impulse response
            s_e[n]:[s_0[n] 0{M-1} s_1[n] .... s_N-1[n] 0{M-1}]
            x[n] = h (*) (IFFT(s_e[n]))
            x_gfdm = sum_M(circshift(x[n],nN))
    '''

    h = gfdm_filter_taps(filtertype, alpha, M, K, 1)
    # Initialization of output vector
    x_out = np.zeros(M * K, dtype='complex')
    # circulary move filter window to symbol 0
    h = np.roll(h, -(M * K / 2))
    # for each gfdm-block
    # for each timeslot
    for m in xrange(M):
        # select the K next symbols
        # symbols = np.fft.ifftshift(x[(m*K):(m+1)*K])
        symbols = np.fft.ifftshift(np.array([x[k * M + m] for k in xrange(K)]))
        # transform K symbols to K carriertones in time-domain
        sym_t = np.fft.ifft(symbols)
        sym_te = np.array([])
        # Repeat result M-times in a vector
        for m2 in xrange(M):
            sym_te = np.concatenate((sym_te, sym_t))
        # multipy with transmit filter -> better convolve?
        sym_te = np.convolve(sym_te, h, mode='same')
        # sym_te = np.multiply(sym_te,h)
        # shift result m*K samples to the right and add it up to the result
        # vector
        x_out = np.add(x_out, np.roll(sym_te, m * K))

    return x_out


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
    # h = np.roll(h, h.shape[-1] / 2)
    # H = np.fft.fft(h)
    H = gfdm_freq_taps(h)
    # H_sparse = np.concatenate((H[0:(M * L) / 2], H[-(M * L) / 2:]))
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
    return x_t #, H, H_sparse


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


def reshape_input(x, M, K):
    '''
    1. pick every m*Kth symbol and append to output.
    2. Increase counter one time
    3. perform step 1.+2. M times
    '''
    # would be equivalent: x_out = np.reshape(x, (-1, K)).T.flatten()
    # Yields the correct data matrix: D = np.reshape(x, (-1, K))
    x_out = np.array([])
    for k in xrange(K):
        for m in xrange(M):
            x_out = np.append(x_out, x[(m * K) + k])

    return x_out



def main():

    M = 8
    K = 2
    alpha = 1.0
    oversampling_factor = 1
    t_extract = 2

    taps = gfdm_filter_taps('rrc', alpha, M, K, oversampling_factor)
    print np.shape(taps)
    A0 = gfdm_modulation_matrix(taps, M, K, oversampling_factor, rearrange_indices=True)
    print 'GFDM shape: ', np.shape(A0)
    # plot_gfdm_matrix(A0)

    A1 = transmitMatrix_core(taps, M, K, oversampling_factor)
    scaling_factor = abs(A0[0, 0]) / abs(A1[0, 0])
    A1 *= scaling_factor
    print 'GR shape: ', np.shape(A1)
    print 'scaling factor:', scaling_factor
    plot_gfdm_matrix(A1)

    # plot_gfdm_matrix(abs(A0) - abs(A1))
    plot_gfdm_matrix(A0 - A1)

    fig = plt.figure()
    data = np.arange(M * K)
    x0 = A0.dot(data)
    plt.plot(x0)
    x1 = gfdm_tx_fft2(data, 'rrc', alpha, M, K, 2, 1)
    # x1 = A1.dot(data)
    plt.plot(x1)
    #
    # t0 = A0[:, t_extract]
    # t1 = A1[:, t_extract]
    # for i in range(4, 8):
    #     plt.plot(A0[:, i].real)
    #     plt.plot(A0[:, i].imag, '--')
    #     plt.plot(A1[:, i].real, '-.')

    plt.show()

if __name__ == '__main__':
    main()
