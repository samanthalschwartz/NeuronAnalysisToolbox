/** @file HSData_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-21-2014
 * @brief The class declaration and inline and templated functions for HSData_Iface.
 */
#include "HSData_Iface.h"
#include "Handle.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    HSData_Iface iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}

HSData_Iface::HSData_Iface()
    : Mex_Iface("HSData_Iface")
{
    methodmap["loadFrames"]=boost::bind(&HSData_Iface::objLoadFrames, this);
    methodmap["loadRawFrames"]=boost::bind(&HSData_Iface::objLoadRawFrames, this);
}

void HSData_Iface::objConstruct()
{
    checkNumArgs(1,4);
    auto raw_size=getIVec();
    auto gain=getDouble();
    auto bg=getFMat();
    auto filepath=getString();
    auto *hsd=new HSData(raw_size(0), raw_size(1), raw_size(2), raw_size(3), gain, bg, filepath);
    outputMXArray(Handle<HSData>::makeHandle(hsd));
}

void HSData_Iface::objDestroy()
{
    checkNumArgs(0,1);
    Handle<HSData>::destroyObject(rhs[0]);
}

void HSData_Iface::getObjectFromHandle(const mxArray *mxhandle)
{
    obj=Handle<HSData>::getObject(mxhandle);
}

void HSData_Iface::objLoadFrames()
{
    checkNumArgs(1,0);
    auto im=makeFHyperStack(obj->sizeL, obj->sizeY, obj->sizeX, obj->nT);
    PixelT *buf=static_cast<PixelT*>(mxGetData(lhs[0]));
    obj->loadAllFrames(buf);
}

void HSData_Iface::objLoadRawFrames()
{
    checkNumArgs(1,0);
    auto im=makeU16HyperStack(obj->raw_sizeL, obj->sizeY, obj->sizeX, obj->nT);
    RawPixelT *buf=static_cast<RawPixelT*>(mxGetData(lhs[0]));
    obj->loadAllRawFrames(buf);
}


