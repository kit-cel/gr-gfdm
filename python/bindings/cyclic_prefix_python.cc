/*
 * Copyright 2020 Johannes Demel
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
/* BINDTOOL_HEADER_FILE(add_cyclic_prefix_cc.h)                                    */
/* BINDTOOL_HEADER_FILE_HASH(70316a02a0eb724098920d68538988c2)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>

#include <gfdm/add_cyclic_prefix_cc.h>


namespace py = pybind11;


void bind_cyclic_prefixer(py::module& m)
{
    using namespace gr::gfdm;
    py::class_<add_cyclic_prefix_cc>(m, "Cyclic_prefixer")

        .def(py::init<int,
                      int,
                      int,
                      int,
                      std::vector<add_cyclic_prefix_cc::gfdm_complex>,
                      int>(),
             py::arg("block_len"),
             py::arg("cp_len"),
             py::arg("cs_len"),
             py::arg("ramp_len"),
             py::arg("window_taps"),
             py::arg("cyclic_shift") = 0)
        .def("block_size", &add_cyclic_prefix_cc::block_size)
        .def("frame_size", &add_cyclic_prefix_cc::frame_size)
        .def("cyclic_shift", &add_cyclic_prefix_cc::cyclic_shift)
        .def("add_cyclic_prefix",
             [](add_cyclic_prefix_cc& self,
                const py::array_t<add_cyclic_prefix_cc::gfdm_complex,
                                  py::array::c_style | py::array::forcecast> array) {
                 py::buffer_info inb = array.request();
                 if (inb.ndim != 1) {
                     throw std::runtime_error("Only ONE-dimensional vectors allowed!");
                 }
                 if (inb.size != self.block_size()) {
                     throw std::runtime_error(
                         "Input vector size(" + std::to_string(inb.size) +
                         ") MUST be equal to Cyclic_prefix.block_size(" +
                         std::to_string(self.block_size()) + ")!");
                 }
                 auto result =
                     py::array_t<add_cyclic_prefix_cc::gfdm_complex>(self.frame_size());
                 py::buffer_info resb = result.request();

                 self.generic_work((add_cyclic_prefix_cc::gfdm_complex*)resb.ptr,
                                   (add_cyclic_prefix_cc::gfdm_complex*)inb.ptr);
                 return result;
             })
        .def("remove_cyclic_prefix",
             [](add_cyclic_prefix_cc& self,
                const py::array_t<add_cyclic_prefix_cc::gfdm_complex,
                                  py::array::c_style | py::array::forcecast> array) {
                 py::buffer_info inb = array.request();
                 if (inb.ndim != 1) {
                     throw std::runtime_error("Only ONE-dimensional vectors allowed!");
                 }
                 if (inb.size != self.frame_size()) {
                     throw std::runtime_error(
                         "Input vector size(" + std::to_string(inb.size) +
                         ") MUST be equal to Cyclic_prefix.frame_size(" +
                         std::to_string(self.block_size()) + ")!");
                 }
                 auto result =
                     py::array_t<add_cyclic_prefix_cc::gfdm_complex>(self.block_size());
                 py::buffer_info resb = result.request();

                 self.remove_cyclic_prefix((add_cyclic_prefix_cc::gfdm_complex*)resb.ptr,
                                           (add_cyclic_prefix_cc::gfdm_complex*)inb.ptr);
                 return result;
             });
}
