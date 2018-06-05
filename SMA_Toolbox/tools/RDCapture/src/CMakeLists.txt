# SRRender -  CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (07/2016)
#
set(LIB_SRCS bessel.cpp rombergint.cpp laplacetform.cpp rdcapture.cpp)
set(TEST_SRCS test_rdcapture.cpp)

set(LIB_TARGET rdcapture${DEBUG_FILE_EXT})
set(TEST_TARGET test_rdcapture${DEBUG_FILE_EXT})
set(IFACE_TARGET RDCapture_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
target_link_libraries(${LIB_TARGET} ${GSL_LIBRARY})
