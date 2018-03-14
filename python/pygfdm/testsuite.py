#!/usr/bin/env python2.7

import numpy as np
import commpy as cp
from . import modulation as mod
from . import utils as ut
import multiprocessing
import functools


class Tester():

    def __init__(self, **kwargs):
        self.env = dict([])
        self.env['m'] = 4
        self.env['k'] = 64
        self.env['l'] = 2
        self.env['qam'] = 4
        self.env['j'] = 0
        self.env['cp'] = 10
        self.env['blocks'] = 10
        self.env['filter'] = 'rrc'
        self.env['rolloff'] = 0.35
        self.env['snr'] = 10
        self.env['os'] = 1
        self.env['fo'] = 0.25
        if kwargs is not None:
            # Non default initialization
            for key, value in kwargs.iteritems():
                self.env[key] = value
        self.funcs = dict([])
        self.funcs['tx'] = mod.gfdm_tx_fft2
        self.funcs['rx'] = mod.gfdm_rx_fft2
        self.funcs['sync_symbol'] = mod.sync_symbol2

    def Transceiver(self):
        import matplotlib.pyplot as mp
        mk = self.env['m']*self.env['k']
        self.sym = mod.randomQAMSymbols(mk*self.env['blocks'], self.env['qam'])
        tx_array = np.array([])
        rx_array = np.array([])
        for b in xrange(self.env['blocks']):
            x = self.sym[b*mk:(b+1)*mk]
            tx = self.funcs['tx'](x,
                                  self.env['filter'],
                                  self.env['rolloff'],
                                  self.env['m'],
                                  self.env['k'],
                                  self.env['l'],
                                  self.env['os'])
            tx_array = np.concatenate((tx_array, tx))
            rx = self.funcs['rx'](tx,
                                  self.env['filter'],
                                  self.env['rolloff'],
                                  self.env['m'],
                                  self.env['k'],
                                  self.env['l'],
                                  self.env['os'],
                                  self.env['qam'],
                                  self.env['j'])
            rx_array = np.concatenate((rx_array, rx))
        (fig, ax) = mp.subplots(1, 2)
        fig.show()
        fig.set_size_inches((10, 5), forward=True)
        ax[0].scatter(self.sym.real, self.sym.imag)
        ax[0].scatter(rx_array.real, rx_array.imag, c=u'r')
        ax[0].set_title('Konstellation')
        ax[0].xaxis.set_label_text('Inphase')
        ax[0].yaxis.set_label_text('Quadratur')
        # ax[0][1].plot(np.abs(tx))
        ax[1].plot(self.sym.real)
        ax[1].plot(rx_array.real)
        ax[1].set_xlim([0, 50])
        ax[1].set_title('Symbolfolge')
        ax[1].xaxis.set_label_text('Zeit')
        ax[1].yaxis.set_label_text('Re{d}')
        # ax[1][1].plot(np.abs(np.fft.fftshift(np.fft.fft(tx_array[0:mk])))/self.env['k'])
        fig.canvas.draw()

    def Tranceiver_nogui(self):
        mk = self.env['m']*self.env['k']
        self.sym = mod.randomQAMSymbols(mk*self.env['blocks'], self.env['qam'])
        tx_array = np.array([])
        rx_array = np.array([])
        for b in xrange(self.env['blocks']):
            x = self.sym[b*mk:(b+1)*mk]
            tx = self.funcs['tx'](x,
                                  self.env['filter'],
                                  self.env['rolloff'],
                                  self.env['m'],
                                  self.env['k'],
                                  self.env['l'],
                                  self.env['os'])
            tx_array = np.concatenate((tx_array, tx))
            rx = self.funcs['rx'](tx,
                                  self.env['filter'],
                                  self.env['rolloff'],
                                  self.env['m'],
                                  self.env['k'],
                                  self.env['l'],
                                  self.env['os'],
                                  self.env['qam'],
                                  self.env['j'])
            rx_array = np.concatenate((rx_array, rx))

    def Transceiver_sync(self):
        import matplotlib.pyplot as mp
        fo = self.env['fo']/(2*self.env['k']*self.env['os'])
        mk = self.env['m']*self.env['k']
        self.sym = mod.randomQAMSymbols(mk*self.env['blocks'], self.env['qam'])
        tx_array = np.array([])
        tx_sync = mod.add_cp(
            self.funcs['sync_symbol'](
                self.env['filter'],
                self.env['rolloff'],
                self.env['k'],
                int(np.log2(self.env['qam'])),
                self.env['os']),
            self.env['cp']*self.env['os'])
        for b in xrange(self.env['blocks']):
            x = self.sym[b*mk:(b+1)*mk]
            tx = self.funcs['tx'](x,
                                  self.env['filter'],
                                  self.env['rolloff'],
                                  self.env['m'],
                                  self.env['k'],
                                  self.env['l'],
                                  self.env['os'])
            tx = mod.add_cp(tx, self.env['cp']*self.env['os'])
            tx_array = np.concatenate((tx_array, tx_sync, tx))

        tx_array = cp.add_frequency_offset(tx_array, 1, fo)
        tx_array = np.concatenate(
            (tx_array, np.zeros(
                 int(3 * mk / 4) * self.env['os'],
                 dtype='complex')))
        (P_di, P_d) = mod.sync_perform(tx_array,
                                       self.env['k']*self.env['os'],
                                       self.env['cp']*self.env['os'])
        d = np.array([])
        d_f = np.array([])
        bs = (mk+2*self.env['k']+self.env['cp'])*self.env['os']
        rx_array = np.array([])
        for b in xrange(int(np.floor(tx_array.shape[-1]/bs))):
            (db, d_fb) = mod.sync_CFO(P_d[b*bs:(b+1)*bs],
                                      P_di[b*bs:(b+1)*bs])
            # Block start minus Circular Prefix
            db = db+b*bs
            d = np.concatenate((d, np.array([db])))
            d_f = np.concatenate((d_f, np.array([d_fb])))
            # Start of the 'real' gfdm block
            gfdm_bs = db+self.env['os']*(2*self.env['k']+self.env['cp'])
            rx = self.funcs['rx'](
                tx_array[gfdm_bs:gfdm_bs+mk*self.env['os']],
                self.env['filter'],
                self.env['rolloff'],
                self.env['m'],
                self.env['k'],
                self.env['l'],
                self.env['os'],
                self.env['qam'],
                self.env['j'])
            rx_array = np.concatenate((rx_array, rx))
        (fig, ax) = mp.subplots(2, 2)
        fig.show()
        #ax[0][0].stem(d, np.ones(d.shape[-1]))
        ax[0][0].plot(self.sym.real)
        ax[0][0].plot(rx_array.real)
        # ax[0][1].plot(d_f)
        ax[0][1].plot(tx_array.real)
        ax[1][0].stem(d, np.ones(d.shape[-1]))
        ax[1][0].plot(P_di)
        ax[1][0].plot(np.abs(P_d)**2)
        ax[1][1].scatter(rx_array.real, rx_array.imag)

        fig.canvas.draw()

    def Transceiver_iter(self):
        mk = self.env['m']*self.env['k']
        self.sym = mod.randomQAMSymbols(
            mk*self.env['blocks'], self.env['qam'])
        eb_n0 = dict()
        if self.env['j'] > 0:
            for j_it in xrange(self.env['j']):
                tx_array = np.array([])
                rx_array = np.array([])

                for b in xrange(self.env['blocks']):
                    x = self.sym[b*mk:(b+1)*mk]
                    tx = self.funcs['tx'](x,
                                          self.env['filter'],
                                          self.env['rolloff'],
                                          self.env['m'],
                                          self.env['k'],
                                          self.env['l'],
                                          self.env['os'])
                    tx_array = np.concatenate((tx_array, tx))
                    rx = self.funcs['rx'](tx,
                                          self.env['filter'],
                                          self.env['rolloff'],
                                          self.env['m'],
                                          self.env['k'],
                                          self.env['l'],
                                          self.env['os'],
                                          self.env['qam'],
                                          j_it)
                    rx_array = np.concatenate((rx_array, rx))
                sc_interference = self.sym - rx_array
                mw_real = (1.0/len(sc_interference))*np.sum(sc_interference.real)
                mw_imag = (1.0/len(sc_interference))*np.sum(sc_interference.imag)
                var_real = (1.0/len(sc_interference))*np.sum(sc_interference.real**2)
                var_imag = (1.0/len(sc_interference))*np.sum(sc_interference.imag**2)
                mw_sym = (1.0/len(self.sym))*(np.sum(np.abs(self.sym)))
                var_sum = var_real+var_imag
                eb_n0[j_it] = 10*np.log10(mw_sym/(2*np.log2(self.env['qam'])*var_sum))

        return eb_n0

    def sc_interference(self, k_range, m_range, qam_range, j):
        self.env['j'] = j
        sc_int = dict()
        for k in k_range:
            self.env['k'] = k
            m_int = dict()
            for m in m_range:
                self.env['m'] = m
                qam_int = dict()
                for qam in qam_range:
                    self.env['qam'] = qam
                    print('Calculating k: {},m: {},qam: {}'.format(k,m,qam))
                    qam_int[qam] = self.Transceiver_iter()
                m_int[m] = qam_int
            sc_int[k] = m_int
        return sc_int




def BER_iter(m, k, qam):
    '''
    Calculate BER for different M,K,rolloff and qam
    '''
    rolloff = 0.35
    ber_m = dict([])
    ber_k = dict([])
    for k_i in k:
        print("M:{} ,K:{}".format(m, k_i))
        ber_qam = dict([])
        for qam_i in qam:
            ber_qam[qam_i] = ut.BER(m, k_i, rolloff, qam_i, -10, 30, 0.5)
        ber_k[k_i] = ber_qam
    ber_m[m] = ber_k
    return ber_m


def BER_mp(M_steps, K_steps, QAM_steps):
    M = np.logspace(3, 3+M_steps, num=1+M_steps, base=2, dtype='int')
    K = np.logspace(4, 4+K_steps, num=1+K_steps, base=2, dtype='int')
    QAM = np.logspace(1, QAM_steps, num=QAM_steps, base=4, dtype='int')

    BER_iter_part = functools.partial(BER_iter, k=K, qam=QAM)
    p = multiprocessing.Pool(processes=4)
    ans = p.map(BER_iter_part, M)
    return ans
