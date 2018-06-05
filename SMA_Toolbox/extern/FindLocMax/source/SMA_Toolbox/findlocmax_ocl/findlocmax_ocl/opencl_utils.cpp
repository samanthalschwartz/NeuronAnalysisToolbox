
#include "opencl_utils.hpp"


inline void checkErr(cl_int err, const char * name)
{
	if (err != CL_SUCCESS) {
		mexPrintf("ERROR %d: %s\n",err,name); 
		mexErrMsgTxt("err");
		exit(EXIT_FAILURE);
	}
}

cl::Platform getplatform(cl_device_type devicetype)
{
	cl::vector< cl::Platform > platforms;
	cl::vector< cl::Device > devices;
	int pid;

	cl::Platform::get(&platforms);
	checkErr(platforms.size()!=0 ? CL_SUCCESS : -1, "cl::Platform::get");
	int Nplatforms = platforms.size();
	mexPrintf("Number of Platforms is: %d\n",Nplatforms);
	if (Nplatforms < 1)
		mexErrMsgTxt("An OpenCL platform is required");

	for (int ii=0;ii<(int) platforms.size();ii++){
		std::string platformVendor;
		platforms[ii].getInfo((cl_platform_info)CL_PLATFORM_VENDOR, &platformVendor);
		mexPrintf("Platform %d is by: %s\n",ii,platformVendor.c_str());}

	//show types of devices on each platform
	for (int ii=0;ii<(int) platforms.size();ii++){
		if (CL_SUCCESS==platforms[ii].getDevices(CL_DEVICE_TYPE_GPU,&devices))
			mexPrintf("Platform %d has GPU device type\n",ii);
		if (CL_SUCCESS==platforms[ii].getDevices(CL_DEVICE_TYPE_CPU,&devices))
			mexPrintf("Platform %d has CPU device type\n",ii);}

	//find first platform with GPU
	for (int ii=0;ii<(int) platforms.size();ii++){
		if (CL_SUCCESS==platforms[ii].getDevices(devicetype,&devices)){
			pid=ii;break;}}
		mexPrintf("Using Platform %d \n",pid);

	cl::Platform platform=platforms[pid];
	return platform;
}