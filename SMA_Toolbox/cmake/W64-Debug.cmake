# SMAToolbox CMake build system
# Mark J. Olah (mjo@cs.unm.edu)
# 03-2014
#
# W64 debug configuration
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type (Debug|Release)" FORCE)
set(DEBUG_FILE_EXT ".debug" CACHE STRING "Directory extension for debug or release" FORCE)
set(CMAKE_TOOLCHAIN_FILE ../../../cmake/Toolchain-w64.cmake CACHE FILEPATH "CrossCompile Toolchain" FORCE)
set(CMAKE_INSTALL_PREFIX ../../../mex/mex.w64.debug CACHE FILEPATH "Install location" FORCE)
