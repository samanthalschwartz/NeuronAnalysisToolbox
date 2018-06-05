#
# DEstimator - Main CMakeLists.txt
#
# Mark J. Olah [mjo@cs.unm.edu] (2014)

cmake_minimum_required( VERSION 2.8 )
project(DEstimator CXX)

set(LOCAL_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake CACHE STRING "Local cmake library")

include(${LOCAL_CMAKE_DIR}/SMAToolbox-configure.cmake)
include(${LOCAL_CMAKE_DIR}/SMAToolbox-configure-blas.cmake)
include(${LOCAL_CMAKE_DIR}/MexIface.cmake)

add_subdirectory(src)
