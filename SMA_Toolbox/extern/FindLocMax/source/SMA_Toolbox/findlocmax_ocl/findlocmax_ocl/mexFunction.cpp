#include <windows.h>
#pragma comment(lib, "kernel32.lib")

#include "opencl_utils.hpp"

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define pi 3.141592f

void getbox(float *data,int ii,int N, int sz,int sz0,int sz1,int sz2,int *centers, float* dataout,int* boxstart);

void gauss_inplace(cl::Buffer bufferA, cl::Kernel kernel1,cl::Kernel kernel2, cl::CommandQueue queue, float sigma, int sz0, int sz1, int sz2)
{
	cl_int err;
	float b0,b1,b2,b3,B,q,qq,qqq;
	if (sigma>2.5)
		q=0.98711f*sigma-.96330f;
	else
		q=3.97156f-4.14554f*sqrt(1.0f-0.26891f*sigma);

	qq=q*q;qqq=qq*q;

	b0=1.57825f+2.44413f*q+1.4281f*qq+.422205f*qqq;
	b1=2.44413f*q+2.85619f*qq+1.26661f*qqq;
	b2=-1.4281f*qq-1.26661f*qqq;
	b3=0.422205f*qqq;
	B=1-(b1+b2+b3)/b0;

	// Set arguments to kernel
	kernel1.setArg(0, bufferA);
	kernel1.setArg(1, sz0);
	kernel1.setArg(2, sz1);
	kernel1.setArg(3, sz2);
	kernel1.setArg(4, b0);
	kernel1.setArg(5, b1);
	kernel1.setArg(6, b2);
	kernel1.setArg(7, b3);
	kernel1.setArg(8, B);

	// Run the kernel on specific ND range
	cl::NDRange global(sz1,sz2);
	cl::NDRange local(sz1,1); // thread per column

	err=queue.enqueueNDRangeKernel(kernel1, cl::NullRange, global, local);
	checkErr(err, "enqueueNDRangeKernel");

	// Set arguments to kernel
	kernel2.setArg(0, bufferA);
	kernel2.setArg(1, sz0);
	kernel2.setArg(2, sz1);
	kernel2.setArg(3, sz2);
	kernel2.setArg(4, b0);
	kernel2.setArg(5, b1);
	kernel2.setArg(6, b2);
	kernel2.setArg(7, b3);
	kernel2.setArg(8, B);

	cl::NDRange global2(sz0,sz2);
	cl::NDRange local2(sz0,1);

	err=queue.enqueueNDRangeKernel(kernel2, cl::NullRange, global2, local2);
	checkErr(err, "enqueueNDRangeKernel");

}

void subtract_inplace(cl::Buffer bufferA, cl::Buffer bufferB, cl::Kernel kernel, cl::CommandQueue queue, int sz0, int sz1, int sz2)
{
	cl_int err;
	// Set arguments to kernel
	kernel.setArg(0, bufferA);
	kernel.setArg(1, bufferB);
	kernel.setArg(2, sz0);
	kernel.setArg(3, sz1);
	kernel.setArg(4, sz2);

	// Run the kernel on specific ND range
	cl::NDRange global(sz0,sz1,sz2);
	cl::NDRange local(sz0,1,1);

	err=queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
	checkErr(err, "enqueueNDRangeKernel");	
}

void locmax(cl::Buffer bufferA, cl::Buffer bufferB, cl::Kernel kernel1, cl::Kernel kernel2,cl::CommandQueue queue, int kernelsize, float minval, int sz0, int sz1, int sz2)
{
	cl_int err;
	//first dimension
	kernel1.setArg(0, bufferA);
	kernel1.setArg(1, bufferB);
	kernel1.setArg(2, kernelsize);	
	kernel1.setArg(3, sz0);
	kernel1.setArg(4, sz1);
	kernel1.setArg(5, sz2);
	kernel1.setArg(6, minval);
	cl::NDRange global(sz0,sz1,sz2);
	cl::NDRange local(sz0,1,1);
	err=queue.enqueueNDRangeKernel(kernel1, cl::NullRange, global, local);
	checkErr(err, "enqueueNDRangeKernel");

	//second dimension
	kernel2.setArg(0, bufferB);
	kernel2.setArg(1, bufferA);
	kernel2.setArg(2, kernelsize);	
	kernel2.setArg(3, sz0);
	kernel2.setArg(4, sz1);
	kernel2.setArg(5, sz2);
	kernel2.setArg(6, minval);
	cl::NDRange global2(sz0,sz1,sz2);
	cl::NDRange local2(1,sz1,1);
	err=queue.enqueueNDRangeKernel(kernel2, cl::NullRange, global2, local2);
	checkErr(err, "enqueueNDRangeKernel");
}



//*******************************************************************************************
void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[]) {
	/*!
	*  \brief Entry point in the code for Matlab.  Equivalent to main().
	*  \param nlhs number of left hand mxArrays to return
	*  \param plhs array of pointers to the output mxArrays
	*  \param nrhs number of input mxArrays
	*  \param prhs array of pointers to the input mxArrays.
	*/


	const mwSize *im_size=0;
	float *image, sigma1,sigma2, *out, minval;
	int Ndim,boxsize,cutboxsize;

	if (nrhs != 6)
		mexErrMsgTxt("Input must be image, sigma1, sigma2, locmaxboxsize, minval, cutboxsize");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Image must be single data type");

	Ndim = (int) mxGetNumberOfDimensions(prhs[0]);
	im_size = mxGetDimensions(prhs[0]);

	if (Ndim>3)
		mexErrMsgTxt("Image must 3D or smaller");

	int Nelem=(int) mxGetNumberOfElements(prhs[0]);

	int sz0=(int) im_size[0];
	int sz1=(int) im_size[1];
	int sz2=(int) 1;
	if (Ndim>2) sz2=(int) im_size[2];
	mexPrintf("sz0: %d sz1: %d sz2: %d Ndim %d\n",sz0,sz1,sz2,Ndim);

	//retrieve all inputs
	image = (float*) mxGetData(prhs[0]);
	sigma1 = (float) mxGetScalar(prhs[1]);
	sigma2 = (float) mxGetScalar(prhs[2]);
	boxsize = (int) mxGetScalar(prhs[3]);
	minval = (float) mxGetScalar(prhs[4]);
	cutboxsize = (int) mxGetScalar(prhs[5]);

	//prepare output vector
	plhs[3] = mxCreateNumericArray(Ndim, im_size, mxSINGLE_CLASS, mxREAL);
	out=(float*)mxGetData(plhs[3]);

	//OpenCL Start------------------------------------------------------

	cl_int err;

	cl::vector< cl::Device > devices;
	cl::vector< cl::Kernel > kernels;
	cl::vector< cl::Event  > events;

	cl::Platform platform=getplatform(CL_DEVICE_TYPE_GPU);

	// Make context
	cl_context_properties cprops[3] ={CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 0};
	cl::Context context(CL_DEVICE_TYPE_GPU,cprops,NULL,NULL,&err);
	checkErr(err, "Context::Context()"); 

	// Get a list of devices on this platform
	devices = context.getInfo<CL_CONTEXT_DEVICES>();
	checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");

	// Create a command queue and use the first device
	cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
	checkErr(err, "CommandQueue::CommandQueue()");

	// Read source file
	std::ifstream file("gauss_kernel.cl");
	checkErr(file.is_open() ? CL_SUCCESS:-1, "gauss_kernel.cl");

	std::string prog(std::istreambuf_iterator<char>(file),(std::istreambuf_iterator<char>()));
	cl::Program::Sources source(1,std::make_pair(prog.c_str(), prog.length()+1));
	checkErr(err, "cl::Program::Sources");

	// Make program of the source code in the context
	cl::Program program(context, source);
	err = program.build(devices,"");
	mexPrintf("Build Log:\t %s\n",program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]).c_str());
	checkErr(err, "Program::build()");

	// Build program for these specific devices
	//program.build(devices);

	// get max memory allocation size
	cl_ulong maxMemAllocSize;
	err =devices[0].getInfo<cl_ulong>(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &maxMemAllocSize);
	mexPrintf("CL_DEVICE_MAX_MEM_ALLOC_SIZE:\t %d\n",maxMemAllocSize);
	checkErr(err, "getInfo");

	int maxelem=maxMemAllocSize/4;
	int nloops=(int) ceil((float) Nelem/((float)maxelem));
	//keep as complete images (multiples of sz2)
	int sz2_sub=floor((float)sz2/nloops);
	// Create memory buffers
	cl::Buffer bufferA = cl::Buffer(context, CL_MEM_READ_WRITE, sz0*sz1*sz2_sub * sizeof(float));
	cl::Buffer bufferB = cl::Buffer(context, CL_MEM_READ_WRITE, sz0*sz1*sz2_sub * sizeof(float));
	// Make kernels
	cl::Kernel kernel0(program, "gauss_kernel0");
	cl::Kernel kernel1(program, "gauss_kernel1");
	cl::Kernel kernel2(program, "subtract_kernel");
	cl::Kernel kernel3(program, "maxdim0_kernel");
	cl::Kernel kernel4(program, "maxdim1_kernel");
	checkErr(err, "Kernel::Kernel()");

	for (int ii=0;ii<nloops;ii++)
	{
		int sz2_inst=min(sz2_sub,sz2-ii*sz2_sub);	
		mexPrintf("sz2_inst:\t %d\n",sz2_inst);

		// Copy image to the memory buffer A and B
		
		err=queue.enqueueWriteBuffer(bufferA, CL_TRUE, 0, sz0*sz1*sz2_inst * sizeof(float), &image[ii*sz0*sz1*sz2_sub],NULL,&events[1]);
		const cl::vector<cl::Event> e1=events;
		err=queue.enqueueCopyBuffer(bufferA, bufferB, 0, 0, sz0*sz1*sz2_inst * sizeof(float),&e1);
		checkErr(err, "enqueueWriteBuffer");

		//run kernels
		gauss_inplace(bufferA, kernel0,kernel1, queue,   sigma1, sz0, sz1, sz2_inst);
		gauss_inplace(bufferB, kernel0,kernel1, queue, sigma2, sz0, sz1, sz2_inst);
		subtract_inplace(bufferA, bufferB, kernel2, queue, sz0, sz1, sz2_inst);
		locmax(bufferA, bufferB, kernel3, kernel4, queue, boxsize, minval, sz0, sz1, sz2_inst);

		// Read buffer C into a local list
		err=queue.enqueueReadBuffer(bufferA, CL_TRUE, 0, sz0*sz1*sz2_inst * sizeof(float), &out[ii*sz0*sz1*sz2_sub]);
		checkErr(err, "enqueueReadBuffer");
	}

	//OpenCL End---------------------------------------------------------
	// find coordinates and write back

	int N=0;
	mwSize outsize1[2];

	for (int ii=0;ii<sz0*sz1*sz2;ii++)
		if (out[ii])N++;

	outsize1[0]=N;
	outsize1[1]=3;
	plhs[2]=mxCreateNumericArray(2, outsize1, mxINT32_CLASS, mxREAL);
	int* centers=(int*)mxGetData(plhs[2]);
	int count=0;

	// here it must be jj*sizeX not jj*sizeY in the if statement (BR)
	for (int kk=0;kk<sz2;kk++)for (int jj=0;jj<sz1;jj++)for (int ii=0;ii<sz0;ii++){
		if (out[kk*sz1*sz0 + jj*sz0 + ii]) {
			centers[count]=ii;
			centers[count+N]=jj;
			centers[count+2*N]=kk;
			count++;
		}}

	const mwSize outsize[3]={cutboxsize,cutboxsize,N};
	plhs[0]=mxCreateNumericArray(3,outsize,mxSINGLE_CLASS,mxREAL);
	float *dataout=(float *)mxGetData(plhs[0]);

	const mwSize outsize2[3]={N,3};
	plhs[1]=mxCreateNumericArray(2,outsize2,mxINT32_CLASS,mxREAL);
	int * boxstart=(int *)mxGetData(plhs[1]);

	for (int ii=0;ii<N;ii++)
		getbox(image,ii,N,cutboxsize,sz0,sz1,sz2,centers,dataout,boxstart);

	return;
}

void getbox(float *data,int ii,int N, int cutboxsize,int sz0,int sz1,int sz2,int *centers, float* dataout,int* boxstart)
{
	//This function copies the specified subregion in data to dataout 

	int jj;
	int start0,start1,end0,end1;

	const int szl=(int) floor(cutboxsize/2.0+0.5);

	//get coordinates
	const int x=(int) centers[ii]; 
	const int y=(int) centers[N+ii]; 
	const int z=(int) centers[2*N+ii]; 

	start0=max(x-szl+1,0);
	end0=start0+cutboxsize-1;
	if (end0>(sz0-1)){end0=sz0-1;start0=end0-sz0+1;}

	start1=max(y-szl+1,0);
	end1=start1+cutboxsize-1;
	if (end1>(sz1-1)){end1=sz1-1;start1=end1-sz1+1;}

	for (jj=0;jj<cutboxsize;jj++) 
		memcpy(dataout+(cutboxsize*cutboxsize*ii+cutboxsize*jj),data+(sz1*sz0*z+sz0*(start1+jj)+start0),cutboxsize*sizeof(float));

	boxstart[ii]=start0;
	boxstart[N+ii]=start1;
	boxstart[2*N+ii]=z;
	return;
}