# Mappel - HSData CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (07/2014)
#
set(LIB_SRCS HSData.cpp)
set(TEST_SRCS test_hsdata.cpp)

set(LIB_TARGET hsdata${DEBUG_FILE_EXT})
set(TEST_TARGET test_hsdata${DEBUG_FILE_EXT})
set(IFACE_TARGET HSData_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
target_link_libraries(${LIB_TARGET} ${LIBICS} ${ZLIB_LIBRARIES})
