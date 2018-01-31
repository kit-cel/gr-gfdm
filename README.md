GNU Radio GFDM Modulator/Demodulator
================

The *gr-gfdm* project is a Free Software Package which aims to provide an implementation of Generalized Frequency Division Multiplexing (GFDM) in the GNU Radio framework. GFDM is a proposed waveform for use in 5G.

This project was initiated as a Bachelor thesis at the *Communication Engineering Lab (CEL)* at the *Karlsruhe Institute of Technology (KIT)*, Germany, <http://www.cel.kit.edu>.

Concept
-------------
Tackling current issues with *OFDM* several new waveforms are proposed for consideration in 5G. *FBMC*, *UFMC*, *BFDM* and *GFDM* are the names of the proposed waveforms and they are filtered multicarrier systems. Orthogonality of neighboring carriers is no constraint.

*GFDM* was proposed by the Vodafone Chair of TU Dresden and accurately described in [1]. Due to the high complexity of a Transmit-Matrix-based approach a low coplexity receiver [2] and transmitter [3] are proposed.

*GFDM* can be described by parallelizing several *SC-FDE* streams on subcarrier. The transmit symbols are localized in a time/frequency-grid and pulshaping is applied subcarrier-wise. After pulshaping and localizing on the correct subcarrier-frequency the symbolstreams are superpositioned and can be transmitted. On receiver side the symbols on the subcarrier can be extracted by applying a *MF*, *ZF* or *MMSE*-filter of the previous pulseshaping filter. Non-orthogonality of neighboring subcarrier introduces *ICI* if demodulating with *MF*. A successive interference cancellation algorithm is proposed to remove interference.

Due to its block-nature a block synchronisation with improved Schmidl & Cox - Symbols can be achievd.

Capabilities
-------------

For generating a transmit symbol stream *gr-gfdm* provides following blocks:

- resource mapper
- modulator
- cyclic prefixer

For reception of a GFDM symbol stream *gr-gfdm* provides following blocks:

- synchronization
- channel estimator
- prefix remover
- simple demodulator
- advanced demodulator with IC
- demapper 

Signal processing in the blocks is performed in *kernels* which have no dependency on *GNU Radio*. Thus inclusion in interfaces to other programming languages is easier.

The modulator and demodulator are implemented using the low complexity approach with Sparse Frequency Domain processing and heavy use of *FFTW* and *VOLK* to accelerate signal processing on modern *GPP*s.

Tests of optimized C++ code are performed with a reference implementation in Python.
Additionally the C++ code is wrapped with [cython](python/cython/README.md).
 


Requirements
------------
- GNU Radio OOT
  - GNU Radio 3.7.9 (verified)
  - GR 3.7 API
  - GR-FFT
  - GR-FILTER
  - VOLK
  - UHD (examples)
    
- PyGFDM
  - Numpy
  - Scipy
  - Commpy
  - Matplotlib
  


Build/Install instructions
------------------------------------

1. Install/Build GNU Radio with at least support for FFT,FILTER and VOLK

2. Get *gr-gfdm* from github - `git clone https://github.com/kit-cel/gr-gfdm.git`

3. Configure *gr-gfdm* - `mkdir build && cd build && cmake ../`

4. Build and install *gr-gfdm* - `make && sudo make install` (default install target: /usr/local)

5. Configure custom blocks path in GNU Radio Companion to use `/usr/local/share/gnuradio/grc/blocks`

Troubleshooting/Bugs
------------------------------------

In case you encounter bugs not related to *GNU Radio* core functions in *gr-gfdm* open an issue at <https://github.com/kit-cel/gr-gfdm/issues>.
Otherwise consider reporting the issue to the GNU Radio project <https://www.gnuradio.org>.

References
-------------
1. N. Michailow et al. “Generalized Frequency Division Multiplexing for 5th Generation Cellular Networks”. In: Communications, IEEE Transactions on 62.9 (2014), S. 3045–3061. doi: 10.1109/TCOMM.2014.2345566.

2. I.S. Gaspar et al. “Low Complexity GFDM Receiver Based on Sparse Frequency Domain Processing”. In: Vehicular Technology Conference (VTC Spring), 2013 IEEE 77th. IEEE, 2013, S. 1–6. doi: 10.1109/VTCSpring.2013.6692619.

3. N. Michailow et al. “Generalized frequency division multiplexing: Analysis of an alternative multi-carrier technique for next generation cellular systems”. In: Wireless Communication Systems (ISWCS), 2012 International Symposium on. IEEE, 2012, S. 171–175. doi: 10.1109/ISWCS.2012.6328352.

