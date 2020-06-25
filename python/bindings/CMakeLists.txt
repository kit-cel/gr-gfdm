# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

GR_PYTHON_CHECK_MODULE_RAW(
    "pygccxml"
    "import pygccxml"
    PYGCCXML_FOUND
    )

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND gfdm_python_files
    ../../lib/modulator_kernel_cc.cc
    ../../lib/add_cyclic_prefix_cc.cc
    ../../lib/resource_mapper_kernel_cc.cc
    ../../lib/receiver_kernel_cc.cc
    ../../lib/gfdm_kernel_utils.cc
    modulator_cc_python.cc
    modulator_python.cc
    cyclic_prefix_python.cc
    cyclic_prefixer_cc_python.cc
    resource_mapper_python.cc
    demodulator_python.cc
    advanced_receiver_sb_cc_python.cc
    channel_estimator_cc_python.cc
    extract_burst_cc_python.cc
    remove_prefix_cc_python.cc
    resource_demapper_cc_python.cc
    resource_mapper_cc_python.cc
    short_burst_shaper_python.cc
    simple_modulator_cc_python.cc
    simple_receiver_cc_python.cc
    transmitter_cc_python.cc
    python_bindings.cc)

GR_PYBIND_MAKE_OOT(gfdm
   ../..
   gr::gfdm
   "${gfdm_python_files}")

install(TARGETS gfdm_python DESTINATION ${GR_PYTHON_DIR}/gfdm COMPONENT pythonapi)
