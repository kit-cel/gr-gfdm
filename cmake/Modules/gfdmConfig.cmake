INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_GFDM gfdm)

FIND_PATH(
    GFDM_INCLUDE_DIRS
    NAMES gfdm/api.h
    HINTS $ENV{GFDM_DIR}/include
        ${PC_GFDM_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GFDM_LIBRARIES
    NAMES gnuradio-gfdm
    HINTS $ENV{GFDM_DIR}/lib
        ${PC_GFDM_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GFDM DEFAULT_MSG GFDM_LIBRARIES GFDM_INCLUDE_DIRS)
MARK_AS_ADVANCED(GFDM_LIBRARIES GFDM_INCLUDE_DIRS)

