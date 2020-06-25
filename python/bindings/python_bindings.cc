/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <pybind11/pybind11.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

namespace py = pybind11;

// Headers for binding functions
/**************************************/
/* The following comment block is used for
/* gr_modtool to insert function prototypes
/* Please do not delete
/**************************************/
// BINDING_FUNCTION_PROTOTYPES(
void bind_advanced_receiver_sb_cc(py::module& m);
void bind_channel_estimator_cc(py::module& m);
void bind_extract_burst_cc(py::module& m);
void bind_modulator_cc(py::module& m);
void bind_modulator(py::module& m);
void bind_cyclic_prefixer(py::module& m);
void bind_cyclic_prefixer_cc(py::module& m);
void bind_remove_prefix_cc(py::module& m);
void bind_resource_demapper_cc(py::module& m);
void bind_resource_mapper_cc(py::module& m);
void bind_resource_mapper(py::module& m);
void bind_short_burst_shaper(py::module& m);
void bind_demodulator(py::module& m);
void bind_simple_modulator_cc(py::module& m);
void bind_simple_receiver_cc(py::module& m);
void bind_transmitter_cc(py::module& m);
// ) END BINDING_FUNCTION_PROTOTYPES


// We need this hack because import_array() returns NULL
// for newer Python versions.
// This function is also necessary because it ensures access to the C API
// and removes a warning.
void* init_numpy()
{
    import_array();
    return NULL;
}

PYBIND11_MODULE(gfdm_python, m)
{
    // Initialize the numpy C API
    // (otherwise we will see segmentation faults)
    init_numpy();

    // Allow access to base block methods
    py::module::import("gnuradio.gr");

    /**************************************/
    /* The following comment block is used for
    /* gr_modtool to insert binding function calls
    /* Please do not delete
    /**************************************/
    // BINDING_FUNCTION_CALLS(
    bind_advanced_receiver_sb_cc(m);
    bind_channel_estimator_cc(m);
    bind_extract_burst_cc(m);
    bind_modulator_cc(m);
    bind_modulator(m);
    bind_cyclic_prefixer(m);
    bind_cyclic_prefixer_cc(m);
    bind_remove_prefix_cc(m);
    bind_resource_demapper_cc(m);
    bind_resource_mapper_cc(m);
    bind_resource_mapper(m);
    bind_short_burst_shaper(m);
    bind_demodulator(m);
    bind_simple_modulator_cc(m);
    bind_simple_receiver_cc(m);
    bind_transmitter_cc(m);
    // ) END BINDING_FUNCTION_CALLS
}
