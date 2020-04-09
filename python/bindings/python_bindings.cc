
/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


void bind_modulator(py::module& m);
void bind_cyclic_prefixer(py::module& m);
void bind_resource_mapper(py::module& m);
void bind_demodulator(py::module& m);


PYBIND11_MODULE(gfdm_python, m)
{
    bind_modulator(m);
    bind_cyclic_prefixer(m);
    bind_resource_mapper(m);
    bind_demodulator(m);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
