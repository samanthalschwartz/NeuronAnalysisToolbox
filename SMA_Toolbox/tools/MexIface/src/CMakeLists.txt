# Mappel - mex wrapper CMakeLists.txt
#
# Mark J. Olah [mjo@cs.unm.edu] (2014)
#
#copyright Mark J. Olah and The Regents of the University of New Mexico (2014).
# This code is free for non-commercial use and modification, provided
# this copyright notice remains unmodified and attached to the code

## Source Files ##
set(MEX_IFACE_SRCS Mex_Iface.cpp MexUtils.cpp explore.cpp)

## Iface Core Object librarary. ##
add_library( iface-core OBJECT ${MEX_IFACE_SRCS} )
