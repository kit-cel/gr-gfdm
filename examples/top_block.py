#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Top Block
# Generated: Tue May 24 12:34:05 2016
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

def struct(data): return type('Struct', (object,), data)()
from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.qtgui import Range, RangeWidget
from optparse import OptionParser
import gfdm
import numpy
import sys


class top_block(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Top Block")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Top Block")
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "top_block")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())

        ##################################################
        # Variables
        ##################################################
        self.samp_rate_range = samp_rate_range = 32e3
        self.samp_rate = samp_rate = 32000
        self.mod_length_key = mod_length_key = "frame_len"
        self.gfdm_var = gfdm_var = struct({"nsubcarrier": 16, "ntimeslots": 16, "fftlen": 16*16, "syncfftlen": 16*2, "cplen": 10, "filteralpha": 0.3, })
        self.gfdm_pregen = gfdm_pregen = gfdm.preamble_generator(gfdm_var.nsubcarrier, gfdm_var.filteralpha, gfdm_var.syncfftlen).base()
        self.gfdm_noise = gfdm_noise = 0
        self.gfdm_ic_iter = gfdm_ic_iter = 0
        
        self.gfdm_constellation = gfdm_constellation = digital.constellation_16qam().base()
        

        ##################################################
        # Blocks
        ##################################################
        self._samp_rate_range_range = Range(0, 2e6, 1000, 32e3, 50)
        self._samp_rate_range_win = RangeWidget(self._samp_rate_range_range, self.set_samp_rate_range, "Sample Rate Range", "counter_slider", float)
        self.top_grid_layout.addWidget(self._samp_rate_range_win, 0,1,1,1)
        self._gfdm_noise_range = Range(0, float(1.0/gfdm_var.nsubcarrier), 0.001, 0, 200)
        self._gfdm_noise_win = RangeWidget(self._gfdm_noise_range, self.set_gfdm_noise, "Noise voltage", "counter_slider", float)
        self.top_grid_layout.addWidget(self._gfdm_noise_win, 0,3,1,1)
        self.gfdm_modulator_cc_0 = gfdm.modulator_cc(gfdm_var.nsubcarrier, gfdm_var.ntimeslots, gfdm_var.filteralpha, gfdm_var.fftlen, gfdm_var.syncfftlen, mod_length_key)
        self._gfdm_ic_iter_range = Range(0, 50, 1, 0, 50)
        self._gfdm_ic_iter_win = RangeWidget(self._gfdm_ic_iter_range, self.set_gfdm_ic_iter, "IC Iteration Count", "counter_slider", int)
        self.top_grid_layout.addWidget(self._gfdm_ic_iter_win, 0,0,1,1)
        self.gfdm_cyclic_prefixer_cc_0 = gfdm.cyclic_prefixer_cc(gfdm_var.cplen, "frame_len")
        self.digital_chunks_to_symbols_xx_0 = digital.chunks_to_symbols_bc((gfdm_constellation.points()), 1)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_tag_debug_0_1 = blocks.tag_debug(gr.sizeof_gr_complex*1, "", mod_length_key); self.blocks_tag_debug_0_1.set_display(True)
        self.blocks_stream_to_tagged_stream_0 = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 1, gfdm_var.nsubcarrier*gfdm_var.ntimeslots, mod_length_key)
        self.analog_random_source_x_0 = blocks.vector_source_b(map(int, numpy.random.randint(0, len(gfdm_constellation.points()), 8192)), True)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_random_source_x_0, 0), (self.digital_chunks_to_symbols_xx_0, 0))    
        self.connect((self.blocks_stream_to_tagged_stream_0, 0), (self.gfdm_modulator_cc_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.blocks_stream_to_tagged_stream_0, 0))    
        self.connect((self.digital_chunks_to_symbols_xx_0, 0), (self.blocks_throttle_0, 0))    
        self.connect((self.gfdm_cyclic_prefixer_cc_0, 0), (self.blocks_tag_debug_0_1, 0))    
        self.connect((self.gfdm_modulator_cc_0, 0), (self.gfdm_cyclic_prefixer_cc_0, 0))    

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "top_block")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_samp_rate_range(self):
        return self.samp_rate_range

    def set_samp_rate_range(self, samp_rate_range):
        self.samp_rate_range = samp_rate_range

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)

    def get_mod_length_key(self):
        return self.mod_length_key

    def set_mod_length_key(self, mod_length_key):
        self.mod_length_key = mod_length_key

    def get_gfdm_var(self):
        return self.gfdm_var

    def set_gfdm_var(self, gfdm_var):
        self.gfdm_var = gfdm_var

    def get_gfdm_pregen(self):
        return self.gfdm_pregen

    def set_gfdm_pregen(self, gfdm_pregen):
        self.gfdm_pregen = gfdm_pregen

    def get_gfdm_noise(self):
        return self.gfdm_noise

    def set_gfdm_noise(self, gfdm_noise):
        self.gfdm_noise = gfdm_noise

    def get_gfdm_ic_iter(self):
        return self.gfdm_ic_iter

    def set_gfdm_ic_iter(self, gfdm_ic_iter):
        self.gfdm_ic_iter = gfdm_ic_iter

    def get_gfdm_constellation(self):
        return self.gfdm_constellation

    def set_gfdm_constellation(self, gfdm_constellation):
        self.gfdm_constellation = gfdm_constellation


def main(top_block_cls=top_block, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
