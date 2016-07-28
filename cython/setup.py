from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import subprocess
import shlex
import inspect
import os


# the actual classes to be compiled!
cython_targets = ["modulator_kernel_cc", "add_cyclic_prefix_cc", ]


def get_pkg_option(option, pkg_name):
    return shlex.split(subprocess.check_output(("pkg-config", "--" + option, pkg_name)))

# a few assumptions are made here!
# libfftw3-dev, cython, pkg-config, VOLK are installed.
# cython is a top-level dir of gr-gfdm!
# header files are in 'include/gfdm'
# source files are in 'lib'

own_file_path = inspect.getfile(inspect.currentframe())
own_file_path = os.path.abspath(own_file_path)  # make sure full path is available!
print 'own_file_path: ', own_file_path
print os.path.abspath(own_file_path)
path, filename = os.path.split(own_file_path)
if not filename == 'setup.py':
    raise ValueError('Assumption FAILED: assumed to run setup.py in order to compile Cython module.')
if not path.endswith('gr-gfdm/cython'):
    print path
    raise ValueError("Assumption FAILED: assumed project structure is 'gr-gfdm/cython'")

project_top_level_dir = path.strip('cython')
source_files = ["gfdm_wrapper.pyx"]
for target in cython_targets:
    header_file = os.path.join(project_top_level_dir, 'include/gfdm', target + '.h')
    source_file = os.path.join(project_top_level_dir, 'lib', target + '.cc')
    if not os.path.exists(header_file) or not os.path.exists(source_file):
        raise ValueError('ERROR: Could not find header(' + header_file + ') and source(' + source_file + ') files!')
    source_files += [source_file, ]
print source_files

# assume those are the only additional libraries to link against.
libraries = ['fftw3f', 'volk']
include_dirs = [os.path.join(project_top_level_dir, 'include/gfdm'), os.path.join(project_top_level_dir, 'include'), ]
library_dirs = []
extra_compile_args = []
extra_link_args = []
for l in libraries:
    lib_arg = get_pkg_option('libs', l)
    comp_arg = get_pkg_option('cflags', l)
    for c in comp_arg:
        r = c.replace('-I', '')
        include_dirs += [r, ]
    for la in lib_arg:
        if la.find('-L') > -1:
            r = la.replace('-L', '')
            library_dirs += [r, ]
    extra_compile_args += comp_arg
    extra_link_args += lib_arg

print 'include_dirs: ', include_dirs
print 'library_dirs: ', library_dirs
print extra_compile_args
print extra_link_args

ext_modules = [
    Extension("gfdm_wrapper",
              include_dirs=include_dirs,
              library_dirs=library_dirs,
              sources=source_files,
              libraries=libraries,
              language="c++",
              )
]

setup(
    ext_modules=cythonize(ext_modules)
    # ext_modules=cythonize(ext_modules, annotate=True)  # this produces an .HTML which helps identify function overhead.
)
