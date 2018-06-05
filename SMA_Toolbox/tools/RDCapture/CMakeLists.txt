#
# SRRender - Main CMakeLists.txt
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
include(${LOCAL_CMAKE_DIR}/MexIface.cmake)

add_subdirectory(src)
