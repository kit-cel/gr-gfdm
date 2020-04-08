
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

#include <gfdm/resource_mapper_kernel_cc.h>


namespace py = pybind11;


void bind_resource_mapper(py::module& m)
{
    using namespace gr::gfdm;
    py::class_<resource_mapper_kernel_cc>(m, "Resource_mapper")

        .def(py::init<int, int, int, std::vector<int>, bool>())
        .def("block_size", &resource_mapper_kernel_cc::block_size)
        .def("frame_size", &resource_mapper_kernel_cc::frame_size)
        .def("map_to_resources", [](resource_mapper_kernel_cc& self,
                            const py::array_t<resource_mapper_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array){
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
            auto result =py::array_t<resource_mapper_kernel_cc::gfdm_complex>(self.frame_size());
            py::buffer_info resb = result.request();

            self.map_to_resources((resource_mapper_kernel_cc::gfdm_complex*) resb.ptr,
                                  (resource_mapper_kernel_cc::gfdm_complex*) inb.ptr,
                                  self.block_size());
            return result;
        })
        .def("demap_from_resources", [](resource_mapper_kernel_cc& self,
                            const py::array_t<resource_mapper_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array){
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            if(inb.size != self.frame_size()){
                throw std::runtime_error("Input vector size(" +
                                         std::to_string(inb.size) +
                                         ") MUST be equal to Modulator.block_size(" +
                                         std::to_string(self.frame_size()) +
                                         ")!");
            }
            auto result =py::array_t<resource_mapper_kernel_cc::gfdm_complex>(self.block_size());
            py::buffer_info resb = result.request();

            self.demap_from_resources((resource_mapper_kernel_cc::gfdm_complex*) resb.ptr,
                                      (resource_mapper_kernel_cc::gfdm_complex*) inb.ptr,
                                      self.block_size());
            return result;
        })
        ;
}
