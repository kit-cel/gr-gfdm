# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

list(APPEND gfdm_sources
    transmitter_cvc_impl.cc
    framer_cc_impl.cc
    gfdm_utils.cc
    gfdm_receiver.cc
    modulator_cc_impl.cc
    #receiver_cc_impl.cc
    receiver2_cc_impl.cc
    advanced_receiver_cc_impl.cc
    sync_cc_impl.cc
    cyclic_prefixer_cc_impl.cc)

set(gfdm_sources "${gfdm_sources}" PARENT_SCOPE)
if(NOT gfdm_sources)
	MESSAGE(STATUS "No C++ sources... skipping lib/")
	return()
endif(NOT gfdm_sources)

add_library(gnuradio-gfdm SHARED ${gfdm_sources})
target_link_libraries(gnuradio-gfdm ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES})
set_target_properties(gnuradio-gfdm PROPERTIES DEFINE_SYMBOL "gnuradio_gfdm_EXPORTS")

if(APPLE)
    set_target_properties(gnuradio-gfdm PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-gfdm
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_gfdm_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_gfdm.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_gfdm.cc
)

add_executable(test-gfdm ${test_gfdm_sources})

target_link_libraries(
  test-gfdm
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  gnuradio-gfdm
  ${GNURADIO_FFT_LIBRARIES}
  ${GNURADIO_FILTER_LIBRARIES}
)

GR_ADD_TEST(test_gfdm test-gfdm)
