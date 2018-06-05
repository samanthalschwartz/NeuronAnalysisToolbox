# Mappel - Boxxer CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (07/2014)
#
set(LIB_SRCS filter_kernels.cpp Boxxer2D.cpp Boxxer3D.cpp GaussFilter.cpp Maxima.cpp)
set(TEST_SRCS test_boxxer.cpp)

set(LIB_TARGET boxxer${DEBUG_FILE_EXT})
set(TEST_TARGET test_boxxer${DEBUG_FILE_EXT})
set(IFACE_TARGET Boxxer2D_Iface Boxxer3D_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
