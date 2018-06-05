

#define __NO_STD_VECTOR // Use cl::vector instead of STL version
#include "cl.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <mex.h>

cl::Platform getplatform(cl_device_type devicetype);
inline void checkErr(cl_int err, const char * name);
