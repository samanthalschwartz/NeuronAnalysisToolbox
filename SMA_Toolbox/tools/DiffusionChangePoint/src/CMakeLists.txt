# DiffusionChangePoint CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (07/2015)

set(LIB_SRCS  interaction_change_point.cpp)
# set(TEST_SRCS lap_test.cpp)

set(LIB_TARGET diffusionChangePoint${DEBUG_FILE_EXT})

set(IFACE_TARGET InteractionChangePoint_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
