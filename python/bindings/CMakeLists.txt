########################################################################
# Python Bindings
########################################################################

pybind11_add_module(gfdm_python
    ${CMAKE_SOURCE_DIR}/lib/modulator_kernel_cc.cc
    ${CMAKE_SOURCE_DIR}/lib/add_cyclic_prefix_cc.cc
    modulator_python.cc
    cyclic_prefix_python.cc
    python_bindings.cc
)

target_include_directories(gfdm_python
    PUBLIC ${PYTHON_NUMPY_INCLUDE_DIR}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
    PUBLIC $<INSTALL_INTERFACE:include>
    ${PYBIND11_INCLUDE_DIR}
)


target_link_libraries(gfdm_python PUBLIC
    ${PYTHON_LIBRARIES}
    gnuradio-gfdm
    gnuradio::gnuradio-runtime
    gnuradio-digital
    ${FFTW3f_LIBRARIES}
    ${VOLK_LIBARIES}
)

if(CMAKE_COMPILER_IS_GNUCC AND NOT APPLE)
    SET_TARGET_PROPERTIES(gfdm_python
    PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
endif()

install(TARGETS gfdm_python DESTINATION ${GR_PYTHON_DIR}/gfdm COMPONENT pythonapi)
