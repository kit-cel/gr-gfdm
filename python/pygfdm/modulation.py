#!/usr/bin/env python2.7

import numpy as np
import commpy as cp

__all__ = [
    'transmitMatrix',
    'fourierMatrix',
    'samplingMatrix',
    'randomQAMSymbols',
    'gfdm_tx',
    'gfdm_rx']

# FIXME TransmitMatrix should group different subcarriers on timeslot-basis


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
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M*K*N, alpha, N*K, 1)
    elif filtertype == "rc":
        time_h, h = cp.rcosfilter(M*K*N, alpha, N*K, 1)
    # Move filter cyclic
    G_tx = np.array([np.roll(h, m - (M*K*N/2)) for m in xrange(M*K*N)])
    S_mn = samplingMatrix(M, K*N)
    S_nm = samplingMatrix(K, M)
    if N > 1:
        # if oversampling is specified add zeros to samplingMatrix
        S_nm = np.insert(
            S_nm, M * K / 2, np.zeros((M * K * (N - 1), K), dtype='complex'),
            axis=0)
    W_H = fourierMatrix(M*K*N).conj().transpose()
    # Resample Filter
    G_tx_s = np.dot(G_tx, S_mn)
    # Resample FourierMatrix
    W_s = np.dot(S_nm.transpose(), W_H)
    # compute and use all elements of the main diagonal W_s.dot(G_tx_s)
    A = np.array([(np.kron(W_s.transpose()[n], G_tx_s[n]))
                 for n in xrange(K*M*N)])

    return A


def fourierMatrix(N):
    i, j = np.meshgrid(np.arange(N), np.arange(N))
    omega = np.exp(- 2 * np.pi * 1j / N)
    W = np.power(omega, i * j) / np.sqrt(N)
    return W


def samplingMatrix(M, K):
    output = np.zeros((M*K, M), dtype='complex')
    for n in xrange(M*K):
        for m in xrange(M):
            if n == ((m)*K):
                output[n][m] = 1
    return output


def randomQAMSymbols(length, M):
    '''
     length: number of symbols to generate
     M: M-QAM - Order (4,16,64,...)
    '''
    n = np.sqrt(M/4)
    if np.around(n) - n > 1e-10:
        raise Exception('M must be power of 4')
    n = int(n)
    n_M_pos = np.array([1+2*i for i in xrange(n)])
    n_M_neg = np.array([-1-2*i for i in xrange(n)])
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
    A = A*M
    tx =A.dot(x)
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
    #if rx_strat == "zf":
    #    A_rx = np.linalg.pinv(A)
    #else:
    #A_rx = np.linalg.pinv(A)/M
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
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M*K, alpha, K, 1)
    elif filtertype == "rc":
        time_h, h = cp.rcosfilter(M*K, alpha, K, 1)
    # Initialization of output vector
    x_out = np.zeros(M*K, dtype='complex')
    # circulary move filter window to symbol 0
    h = np.roll(h, -(M*K/2))
    # for each gfdm-block
        # for each timeslot
    for m in xrange(M):
        # select the K next symbols
        #symbols = np.fft.ifftshift(x[(m*K):(m+1)*K])
        symbols = np.fft.ifftshift(np.array([x[k*M+m] for k in xrange(K)]))
        # transform K symbols to K carriertones in time-domain
        sym_t = np.fft.ifft(symbols)
        sym_te = np.array([])
        # Repeat result M-times in a vector
        for m2 in xrange(M):
            sym_te = np.concatenate((sym_te, sym_t))
        # multipy with transmit filter -> better convolve?
        sym_te = np.convolve(sym_te, h, mode='same')
        #sym_te = np.multiply(sym_te,h)
        # shift result m*K samples to the right and add it up to the result
        # vector
        x_out = np.add(x_out, np.roll(sym_te, m*K))

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
    if filtertype == "rrc":
        time_h, h = cp.rrcosfilter(M*K, alpha, K, 1)
    elif filtertype == "rc":
        time_h, h = cp.rcosfilter(M*K, alpha, K, 1)
    h = np.roll(h, h.shape[-1]/2)
    H = np.fft.fft(h)
    H_sparse = np.concatenate((H[0:(M*L)/2], H[-(M*L)/2:]))
    # Sort Input subcarrierwise
    x = reshape_input(x, M, K)
    x_out = np.zeros((M*K)+(L-1)*M, dtype='complex')
    for k in xrange(K):
        # M rows and L columns with respective FFT output
        # pick symbols per subcarrier
        x_k = x[k*M:((k+1)*M)]
        # perform fft and switch to frequency domain
        x_f = np.fft.fft(x_k)
        # copy values of M-point DFT to obtain MK-point DFT
        x_f_L = np.tile(x_f, L)
        # make data-vector 'sparse'
        #x_f_L = np.concatenate((x_f_K[0:(M*L)/2], x_f_K[-(M*L)/2:]))
        # filter with sparse filter taps in frequency domain
        x_fil = np.multiply(x_f_L, H_sparse)
        # Add data-vector to correct position -max neg frequency : 0 :
        # max_pos_frequency
        x_out[k*M:(k+L)*M] = x_out[k*M:(L+k)*M] + np.fft.fftshift(x_fil)
    # Add 'oversampled' parts of first subcarrier to end and 'oversampled' parts
    # of last subcarrier to start
    x_first = x_out[0:(L-1)*M/2]
    x_last = x_out[-(L-1)*M/2:]
    x_out = x_out[(L-1)*M/2:-(L-1)*M/2]
    x_out[0:(L-1)*M/2] = x_out[0:(L-1)*M/2] + x_last
    x_out[-(L-1)*M/2:] = x_out[-(L-1)*M/2:] + x_first
    x_t = np.fft.ifft(np.fft.ifftshift(x_out))
    x_t = (1.0/K)*x_t
    return x_t


def gfdm_rx_fft2(y, filtertype, alpha, M, K, L, N, QAM,J):
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
        time_h, h = cp.rrcosfilter(M*K, alpha, K, 1)
    h = np.roll(h, h.shape[-1]/2)
    H_rx = np.fft.fft(h)
    H_sparse = np.concatenate((H_rx[0:M*L/2], H_rx[-M*L/2:]))
    y_ifft = np.array([])
    y = (1.0/K)*y
    # Transfer input to frequency domain and center around 0th frequency bin
    y_f = np.fft.fftshift(np.fft.fft(y))
    # Filter and superposition in frequency domain
    Y_fs = gfdm_rx_filsup(y_f, H_sparse, M, K, L)
    # Demodulate per subcarrier
    y_ifft = gfdm_rx_demod(Y_fs, K)
    if J>0:
        y_ifft = gfdm_rx_sic(K,M,J,H_sparse,y_ifft,Y_fs,QAM)
        y_ifft = np.reshape(y_ifft,(K*M))
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
    y_out = np.empty((K,M),dtype='complex')
    y_f = np.concatenate((y_f[-(L-1)*M/2:],y_f,y_f[0:(L-1)*M/2]))
    for k in xrange(K):
        # select kth subcarrier
        y_down = y_f[k*M:(k+L)*M]
        # 'uncenter' in Frequency domain
        y_down = np.fft.ifftshift(y_down)
        # apply filter in frequency domain (not centered)
        y_filter = np.multiply(y_down, H_sparse)
        # Superposition L samples in frequency domain
        y_out[k] = np.sum(y_filter.reshape(L, M), axis=0)
    return y_out



def gfdm_rx_sic(K,M,J,H_sparse,d_rx,Y_fs,QAM):
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
    H_rf = np.multiply((H_sparse[0:M]/K),(H_sparse[M:]/K))
    # Reshape mapped symbols into per-subcarrier array
    d_p = np.empty((K,M),dtype='complex')
    # d_p (K,M)
    for k in xrange(K):
        d_p[k] = qam_mod.mapping(d_rx[k*M:(k+1)*M],'hard')
    for j in xrange(J):
        y = np.empty((K,M),dtype='complex')
        for k in xrange(K):
            y[k] = Y_fs[k] - H_rf*np.fft.fft(d_p[(k-1) % K] + d_p[(k+1) % K])
            # Recalculate d_rx
        d_rx = gfdm_rx_demod(y,K)
        for k in xrange(K):
            d_p[k] = d_rx[k*M:(k+1)*M]
            d_p[k] = qam_mod.mapping(d_p[k], 'hard')
    return d_rx


def sync_symbol(filtertype, alpha, K, n_mod, N):
    '''
        Generate Schmidls Training Symbols to achieve Receiver Synchronisation
        K: should be an odd number
        Process:
            * First Symbol: Transmit PN-Sequence on all even frequencies while zeros on the odd
            frequencies. Constant signal energy -> Multiply every Symbol with sqrt(2)

            * Second Symbol: Transmit PN-Sequence on all odd frequencies and another PN-Sequence
            on the even frequencies


    '''
    pn_order = 14
    pn_seed = '00101101010010'
    pn_mask = '01001110100111'
    if int(np.floor(K/2.0)) % 2:
        n_even_freq = int(np.floor(K/2.0))
    else:
        n_even_freq = int(np.ceil(K/2.0))
    seq_length = n_even_freq*n_mod
    sym_sequence = np.zeros(K, dtype='complex')
    pn_sequence = cp.pnsequence(pn_order, pn_seed, pn_mask, seq_length)
    qam_mod = cp.modulation.QAMModem(2**n_mod)
    qam_sequence = qam_mod.modulate(pn_sequence)
    for i in xrange(len(sym_sequence)):
        if not i % 2:
            sym_sequence[i] = qam_sequence[i/2]
    ifft_sequence = np.fft.ifftshift(sym_sequence)
    output = gfdm_tx(ifft_sequence, filtertype, alpha, 1, K, N)
    # output = np.fft.ifft(np.sqrt(2)*ifft_sequence)

    return output


def sync_symbol2(filtertype, alpha, K, L, n_mod):
    pn_order = 14
    pn_seed = '01001000111011'
    pn_mask = '01001101001110'
    seq_length = K*n_mod
    pn_sequence = cp.pnsequence(pn_order, pn_seed, pn_mask, seq_length)
    qam_mod = cp.modulation.QAMModem(2**n_mod)
    qam_sequence = qam_mod.modulate(pn_sequence)
    output = gfdm_tx_fft2(np.tile(qam_sequence, 2), filtertype, alpha, 2, K, 2, 1)
    return output


def sync_product(x, L):
    '''
    Auto-Korrelation der ersten L Samples mit den naechsten L Samples
    '''
    return np.sum([x[i].conj()*x[i+L] for i in xrange(L)])


def sync_iter(x, L, cp):
    '''
    Schrittweise Iteration ueber alle Samples (zunaechst die ersten 2*L Samples)
    Danach ueber die restlichen len(x)-2*L Samples
    '''
    P_d = np.array([])
    P_d = np.append(P_d, sync_product(x, L))
    for i in xrange(len(x)-2*L):
        P_d = np.append(
            P_d, P_d[i] + (x[L + i].conj() * x[2 * L + i]) -
            (x[i].conj() * x[L + i]))
    P_d_out = P_d
    P_d = np.append(np.zeros(cp, dtype='complex'), P_d)
    # Integrate cp-samples to eliminate cp-plateau
    P_di = np.array([])
    for i in xrange(cp, len(x)-2*L):
        P_di = np.append(
            P_di, (1.0/(cp+1)*np.sum(np.abs(P_d[i-cp:i])**2)))
    return (P_di, P_d_out)


def sync_energy(x, L):
    '''
    Berechnung der Energie der zweiten Haelfte der Sync-Samples -> normieren
    '''
    R_d = np.array([])
    R_d = np.append(R_d, np.sum([np.abs(x[i+L])**2 for i in xrange(L)]))
    for i in xrange(len(x)-2*L):
        R_d = np.append(R_d, R_d[-1]+np.abs(x[2*L+i])**2-np.abs(x[L+i])**2)
    return R_d


def sync_perform(x, L, cp):
    (P_di, P_d) = sync_iter(x, L, cp)
    #R_d = sync_energy(x, L)
    #M_d = (np.abs(P_d)**2)/(R_d**2)
    return (P_di, P_d)


def sync_CFO(P_d, P_di):
    '''
    Gewinn von d (Beginn der Pilotsequenz) und d_f Frequenzoffset genormt auf 1/T.
    Kann bei nur einem Pilotsymbol nur +/- 1/T betragen.
    '''
    d = np.argmax(P_di)
    dphi = np.angle(P_d[d])
    d_f = dphi/(np.pi)
    print("P_d:{},df:{})".format(P_d[d], d_f))

    return (d, d_f)


def add_cp(x, n):
    return np.append(x[-n:], x)


def reshape_input(x, M, K):
    '''
    1. pick every m*Kth symbol and append to output.
    2. Increase counter one time
    3. perform step 1.+2. M times
    '''
    x_out = np.array([])
    for k in xrange(K):
        for m in xrange(M):
            x_out = np.append(x_out, x[(m*K)+k])
    return x_out
