# Tracker CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (04/2015)

set(LIB_SRCS Tracker.cpp LAPTrack.cpp LAP_JVSparse.cpp)
set(TEST_SRCS lap_test.cpp)

set(LIB_TARGET tracker${DEBUG_FILE_EXT})
set(TEST_TARGET test_tracker${DEBUG_FILE_EXT})

set(IFACE_TARGET LAPTrack_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
