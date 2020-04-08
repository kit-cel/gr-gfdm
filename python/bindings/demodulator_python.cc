
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

#include <gfdm/receiver_kernel_cc.h>


namespace py = pybind11;



void bind_demodulator(py::module& m)
{
    using namespace gr::gfdm;

    py::class_<receiver_kernel_cc>(m, "Demodulator")

        .def(py::init<int, int, int, std::vector<receiver_kernel_cc::gfdm_complex> >())
        .def("timeslots", &receiver_kernel_cc::timeslots)
        .def("subcarriers", &receiver_kernel_cc::subcarriers)
        .def("overlap", &receiver_kernel_cc::overlap)
        .def("block_size", &receiver_kernel_cc::block_size)
        .def("filter_taps", &receiver_kernel_cc::filter_taps)
        .def("demodulate", [](receiver_kernel_cc& self,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array){
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
            auto result =py::array_t<receiver_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.generic_work((receiver_kernel_cc::gfdm_complex*) resb.ptr,
                              (receiver_kernel_cc::gfdm_complex*) inb.ptr);
            return result;
        })
        .def("fft_filter_downsample", [](receiver_kernel_cc& self,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array){
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
            auto result =py::array_t<receiver_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.fft_filter_downsample((receiver_kernel_cc::gfdm_complex*) resb.ptr,
                              (receiver_kernel_cc::gfdm_complex*) inb.ptr);
            return result;
        })
        .def("transform_subcarriers_to_td", [](receiver_kernel_cc& self,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array){
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
            auto result =py::array_t<receiver_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.transform_subcarriers_to_td((receiver_kernel_cc::gfdm_complex*) resb.ptr,
                              (receiver_kernel_cc::gfdm_complex*) inb.ptr);
            return result;
        })
        .def("demodulate_equalize", [](receiver_kernel_cc& self,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> eq_arr){
            py::buffer_info inb = array.request();
            py::buffer_info eq_inb = eq_arr.request();
            if(inb.ndim != 1 or eq_inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            if(inb.size != self.block_size()){
                throw std::runtime_error("Input vector size(" +
                                         std::to_string(inb.size) +
                                         ") MUST be equal to Demodulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            if(eq_inb.size != self.block_size()){
                throw std::runtime_error("Channel vector size(" +
                                         std::to_string(eq_inb.size) +
                                         ") MUST be equal to Demodulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            auto result =py::array_t<receiver_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.generic_work_equalize((receiver_kernel_cc::gfdm_complex*) resb.ptr,
                                       (receiver_kernel_cc::gfdm_complex*) inb.ptr,
                                       (receiver_kernel_cc::gfdm_complex*) eq_inb.ptr);
            return result;
        })
        .def("fft_equalize_filter_downsample", [](receiver_kernel_cc& self,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> eq_arr){
            py::buffer_info inb = array.request();
            py::buffer_info eq_inb = eq_arr.request();
            if(inb.ndim != 1 or eq_inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            if(inb.size != self.block_size()){
                throw std::runtime_error("Input vector size(" +
                                         std::to_string(inb.size) +
                                         ") MUST be equal to Demodulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            if(eq_inb.size != self.block_size()){
                throw std::runtime_error("Channel vector size(" +
                                         std::to_string(eq_inb.size) +
                                         ") MUST be equal to Demodulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            auto result =py::array_t<receiver_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.fft_equalize_filter_downsample((receiver_kernel_cc::gfdm_complex*) resb.ptr,
                                                (receiver_kernel_cc::gfdm_complex*) inb.ptr,
                                                (receiver_kernel_cc::gfdm_complex*) eq_inb.ptr);
            return result;
        })
        .def("cancel_sc_interference", [](receiver_kernel_cc& self,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> array,
                            const py::array_t<receiver_kernel_cc::gfdm_complex, py::array::c_style | py::array::forcecast> eq_arr){
            py::buffer_info inb = array.request();
            py::buffer_info eq_inb = eq_arr.request();
            if(inb.ndim != 1 or eq_inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            if(inb.size != self.block_size()){
                throw std::runtime_error("Input vector size(" +
                                         std::to_string(inb.size) +
                                         ") MUST be equal to Demodulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            if(eq_inb.size != self.block_size()){
                throw std::runtime_error("Channel vector size(" +
                                         std::to_string(eq_inb.size) +
                                         ") MUST be equal to Demodulator.block_size(" +
                                         std::to_string(self.block_size()) +
                                         ")!");
            }
            auto result =py::array_t<receiver_kernel_cc::gfdm_complex>(inb.size);
            py::buffer_info resb = result.request();

            self.cancel_sc_interference((receiver_kernel_cc::gfdm_complex*) resb.ptr,
                                                (receiver_kernel_cc::gfdm_complex*) inb.ptr,
                                                (receiver_kernel_cc::gfdm_complex*) eq_inb.ptr);
            return result;
        })
        ;
}
