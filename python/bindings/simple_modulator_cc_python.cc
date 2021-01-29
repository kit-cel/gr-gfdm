/*
 * Copyright 2020 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(simple_modulator_cc.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(34c2be9ad92bac9e8e9eaa7963675e70)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gfdm/simple_modulator_cc.h>
// pydoc.h is automatically generated in the build directory
#include <simple_modulator_cc_pydoc.h>

void bind_simple_modulator_cc(py::module& m)
{

    using simple_modulator_cc = ::gr::gfdm::simple_modulator_cc;


    py::class_<simple_modulator_cc,
               gr::sync_block,
               gr::block,
               gr::basic_block,
               std::shared_ptr<simple_modulator_cc>>(
        m, "simple_modulator_cc", D(simple_modulator_cc))

        .def(py::init(&simple_modulator_cc::make),
             py::arg("n_timeslots"),
             py::arg("n_subcarriers"),
             py::arg("overlap"),
             py::arg("frequency_taps"),
             D(simple_modulator_cc, make))


        ;
}