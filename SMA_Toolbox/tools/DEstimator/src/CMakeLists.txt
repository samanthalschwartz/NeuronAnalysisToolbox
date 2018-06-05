#
# DEstimator - src CMakeLists.txt
#
# Mark J. Olah [mjo@cs.unm.edu] (2014)

set(LIB_SRCS destimator.cpp DCompLib.cpp)

set(TEST_SRCS test_destimator.cpp)
set(TEST_TARGET test_destimator)

set(LIB_TARGET destimator${DEBUG_FILE_EXT})
set(IFACE_TARGET DEstimator_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
