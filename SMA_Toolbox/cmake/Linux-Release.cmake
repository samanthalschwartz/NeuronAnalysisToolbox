# SMAToolbox CMake build system
# Mark J. Olah (mjo@cs.unm.edu)
# 03-2014
#
# Linux release configuration
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug|Release)" FORCE)
set(CMAKE_CXX_COMPILER /usr/bin/g++-4.9.4 CACHE STRING "C++ compiler" FORCE) #This is the latest version used by MATLAB 2015a
set(DEBUG_FILE_EXT "" CACHE STRING "Directory extension for debug or release" FORCE)
set(CMAKE_INSTALL_PREFIX ../../.. CACHE FILEPATH "Install location" FORCE)
