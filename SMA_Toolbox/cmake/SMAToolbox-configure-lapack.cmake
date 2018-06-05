# SMAToolbox CMake build system
# Mark J. Olah (mjo@cs.unm.edu)
# 03-2014
# Include this file to configure LAPACK integration with armadillo
if(WIN32)
    #These compiled libraries are provided by armadillo
    find_library(LAPACK_LIBRARIES lapack_win64_MT.dll)
    set(W64_DLLS ${W64_DLLS} lapack_win64_MT.dll)
else()
    find_package(LAPACK)
endif()
