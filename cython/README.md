Cython module for gr-gfdm
=======

This module is intended to provide direct access to the C++ objects via Python. Unlike GNU Radio's SWIG system it provides an interface for the actual work method. There might be situations where you want to have direct access to the optimized GFDM functions from Python.

Dependencies
=======
It is assumed that the classes wrapped via Cython are NOT depending on GNU Radio. Needed components are:

* VOLK (libvolk-dev or similar)
* FFTW3F (version also used by GNU Radio)
* pkg-config (makes Cython setup for compilation easier)
* 

Functionality
=====
after compilation a *gfdm_wrapper* module is available in Python. Just load it like any other Python module. It provides the interface for several blocks. See *gfdm_wrapper.pyx* for exposed objects and functions.