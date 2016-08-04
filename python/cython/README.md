Cython module for gr-gfdm
=======

This module is intended to provide direct access to the C++ objects via Python. Unlike GNU Radio's SWIG system it provides an interface for the actual work method. There might be situations where you want to have direct access to the optimized GFDM functions from Python.

Dependencies
=======
It is assumed that the classes wrapped via Cython are NOT depending on GNU Radio. Needed components are:

* cython (obviously, since this is the language of choice, version >= 0.23)
* VOLK (libvolk-dev or version from github)
* FFTW3F (version also used by GNU Radio, package: libfftw3-dev)
* pkg-config (makes Cython setup for compilation easier)

Compilation
=====
just run 'python setup.py build' from 'gr-gfdm/python'. If you want all build files to be located in the source directory, append '--inplace'.

Functionality
=====
after compilation a *cgfdm* module is available in Python. Just load it like any other Python module. It provides the interface for several blocks. See *cgfdm.pyx* for exposed objects and functions.

Installation
=====
'python setup.py install' does install *cgfdm* along with *pygfdm*. The '--prefix >my_install_prefix<' option will install it into a custom directory. Afterwards 'import pygfdm' and 'import 'cgfdm' are available.

