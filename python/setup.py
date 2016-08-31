#!/usr/bin/env python

'''
GFDM Cython kernels
~~~~~~~

This setup.py file exists to assist with Cython install.
It may be fragile in case assumptions and dependencies are not met.
It depends on various tools!
Python: setuptools, Cython, subprocess, shlex, inspect, os
C/C++: libfftw3-dev, libvolk-dev (or VOLK from github)
Cython files must be located in ./python/cython!
C/C++ files must be located in ./lib and ./include/gfdm
Additional Python modules must be in python/pygfdm!

This installs 2 additional Python modules!
pygfdm: Python GFDM tools
cgfdm: Cython wrapper to fast C++ GFDM implementation.
'''

from __future__ import print_function, division
from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import subprocess
import shlex
import inspect
import os


# the actual classes to be compiled!
cython_targets = ["modulator_kernel_cc", "add_cyclic_prefix_cc", "resource_mapper_kernel_cc", "receiver_kernel_cc"]
# assume those are the only additional libraries to link against.
libraries = ['fftw3f', 'volk']

deps = [
    'numpy',
    'matplotlib',
    'scipy',
    'scikit-commpy',
]


def get_pkg_option(option, pkg_name):
    return shlex.split(subprocess.check_output(("pkg-config", "--" + option, pkg_name)))


def get_library_shared_objects(lib):
    lib_dirs = []
    lib_arg = get_pkg_option('libs', lib)
    for la in lib_arg:
        if '-L' in la:
            lib_dirs += [la.replace('-L', ''), ]
    return lib_dirs


def find_shared_objects(libraries):
    lib_dirs = []
    for l in libraries:
        lib_dirs += get_library_shared_objects(l)
    return lib_dirs


def get_library_headers(lib):
    shared_dirs = []
    comp_arg = get_pkg_option('cflags', lib)
    for c in comp_arg:
        shared_dirs += [c.replace('-I', ''), ]
    return shared_dirs


def find_headers(libraries):
    shared_dirs = []
    for l in libraries:
        shared_dirs += get_library_headers(l)
    return shared_dirs


def find_source_files(project_top_level_dir, targets):
    source_files = []
    for target in targets:
        header_file = os.path.join(project_top_level_dir, 'include/gfdm', target + '.h')
        source_file = os.path.join(project_top_level_dir, 'lib', target + '.cc')
        if not os.path.exists(header_file) or not os.path.exists(source_file):
            raise ValueError('ERROR: Could not find header(' + header_file + ') and source(' + source_file + ') files!')
        source_files += [source_file, ]
    return source_files


def get_project_top_level_dir():
    own_file_path = inspect.getfile(inspect.currentframe())
    own_file_path = os.path.abspath(own_file_path)  # make sure full path is available!
    # print 'own_file_path: ', own_file_path
    # print os.path.abspath(own_file_path)
    path, filename = os.path.split(own_file_path)
    if not filename == 'setup.py':
        raise ValueError('Assumption FAILED: assumed to run setup.py in order to compile Cython module.')
    if not path.endswith('gr-gfdm/python'):
        print(path)
        raise ValueError("Assumption FAILED: assumed project structure is 'gr-gfdm/cython'")
    return path.strip('python')


def get_available_dependencies(deps):
    avail = []
    for d in deps:
        try:
            __import__(d)
            print(d)
            avail += [d, ]
        except ImportError:
            print('Python module:', d, 'not available!')
    return avail


project_top_level_dir = get_project_top_level_dir()
print('project TOP level directory:', project_top_level_dir)


source_files = [os.path.join(project_top_level_dir, 'python/cython/cgfdm.pyx'), ]
source_files += find_source_files(project_top_level_dir, cython_targets)
print(source_files)

include_dirs = [os.path.join(project_top_level_dir, 'include/gfdm'), os.path.join(project_top_level_dir, 'include'), ]
include_dirs += find_headers(libraries)
library_dirs = find_shared_objects(libraries)

print('include_dirs: ', include_dirs)
print('library_dirs: ', library_dirs)


ext_modules = [
    Extension("cgfdm",
              include_dirs=include_dirs,
              library_dirs=library_dirs,
              sources=source_files,
              libraries=libraries,
              language="c++",
              )
]


# setuptools Doc: https://pythonhosted.org/an_example_pypi_project/setuptools.html
setup(
    name="python-gfdm",
    version="0.0.1",
    author="Johannes Demel, Andrej Rode",
    author_email="demel@ant.uni-bremen.de, andrej.rode@student.kit.edu",
    description='Python GFDM utils',
    keywords="GFDM Python Cython",
    packages=['pygfdm'],
    install_requires=deps,
    ext_modules=cythonize(ext_modules),
    # ext_modules=cythonize(ext_modules, annotate=True)  # this produces an .HTML which helps identify function overhead.
)
