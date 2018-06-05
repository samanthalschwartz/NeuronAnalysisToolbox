# SRRender -  CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (03/2015)
#
set(LIB_SRCS srrender.cpp)
set(TEST_SRCS test_srrender.cpp)

set(LIB_TARGET srrender${DEBUG_FILE_EXT})
set(TEST_TARGET test_srrender${DEBUG_FILE_EXT})
set(IFACE_TARGET SRRender2DDouble_Iface SRRender2DSingle_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
