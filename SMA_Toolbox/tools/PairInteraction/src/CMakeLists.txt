# PairInteraction -  CMakeLists.txt
# Mark J. Olah [mjo@cs.unm.edu] (08/2016)
#
set(LIB_SRCS pairinteractionsmc.cpp unscentedtransform.cpp computations.cpp rngmanager.cpp vonmises.cpp gaussianangle.cpp)
set(TEST_SRCS test_pairinteraction.cpp)

set(LIB_TARGET pairinteraction${DEBUG_FILE_EXT})
set(TEST_TARGET test_pair${DEBUG_FILE_EXT})
set(IFACE_TARGET PairInteractionSMC_Iface)

include(${LOCAL_CMAKE_DIR}/SMAToolbox-Build.cmake)
target_link_libraries(${LIB_TARGET} ${GSL_LIBRARY})

#Replace this with better version later.
if(WIN32)
    find_library(RDCapture librdcapture${DEBUG_FILE_EXT}.dll PATH ${HOME}/LidkeLab/SMA_Toolbox/mex/mex.w64/)
else()
    find_library(RDCapture librdcapture${DEBUG_FILE_EXT}.so)
endif()
target_link_libraries(${LIB_TARGET} ${RDCapture})
# target_link_libraries(${TEST_TARGET} ${RDCapture})
