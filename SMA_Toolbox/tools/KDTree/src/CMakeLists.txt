# Mappel -KDTree CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (08/2014)
#

set(LIB_SRCS kdtree.cpp)
set(TEST_SRCS test_kdtree.cpp)

set(LIB_TARGET kdtree${DEBUG_FILE_EXT})
set(TEST_TARGET test_kdtree${DEBUG_FILE_EXT})

set(IFACE_TARGET KDTreeSingle_Iface KDTreeDouble_Iface KDTreeInt_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
