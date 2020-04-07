
/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <cstdint>

#include <gfdm/modulator_kernel_cc.h>


namespace py = pybind11;


void bind_modulator(py::module& m)
{
    using namespace gr::gfdm;
    py::class_<modulator_kernel_cc>(m, "Modulator")

        .def(py::init<int, int, int, std::vector<modulator_kernel_cc::gfdm_complex> >())
        .def("block_size", &modulator_kernel_cc::block_size)
        .def("filter_taps", &modulator_kernel_cc::filter_taps)
        .def("modulate", [](modulator_kernel_cc& self,
                            const py::array_t<modulator_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array){
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            if(inb.size != self.block_size()){
                throw std::runtime_error("Input vector size(" +
                                         std::to_string(inb.size) +
                                         ") MUST be equal to Modulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            auto result =py::array_t<modulator_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.generic_work((modulator_kernel_cc::gfdm_complex*) resb.ptr,
                              (modulator_kernel_cc::gfdm_complex*) inb.ptr);
            return result;
        })
        ;
}
