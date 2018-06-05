#
# PairInteraction - Main CMakeLists.txt
#
# Mark J. Olah [mjo@cs.unm.edu] (2016)
# copyright Mark J. Olah and The Regents of the University of New Mexico (2014).
#            This code is free for non-commercial use and modification, provided
#            this copyright notice remains unmodified and attached to the code

cmake_minimum_required( VERSION 3.0 )
project(RDCapture CXX)

set(LOCAL_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake CACHE STRING "Local cmake library")

include(${LOCAL_CMAKE_DIR}/SMAToolbox-configure.cmake)
include(${LOCAL_CMAKE_DIR}/SMAToolbox-configure-gsl.cmake)
include(${LOCAL_CMAKE_DIR}/SMAToolbox-configure-blas.cmake)
include(${LOCAL_CMAKE_DIR}/SMAToolbox-configure-lapack.cmake)
include(${LOCAL_CMAKE_DIR}/MexIface.cmake)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../RDCapture/src)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../RDCapture/src)

add_subdirectory(src)
