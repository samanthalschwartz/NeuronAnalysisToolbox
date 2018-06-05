/** @file DEstimator_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-21-2014
 * @brief The class declaration and inline and templated functions for DEstimator_Iface.
 */
#include "DEstimator_Iface.h"
#include "Handle.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    DEstimator_Iface DE_iface;
    DE_iface.mexFunction(nlhs, lhs, nrhs, rhs);
}

DEstimator_Iface::DEstimator_Iface() 
    : Mex_Iface("DEstimator_Iface") 
{
    methodmap["LLH"]=boost::bind(&DEstimator_Iface::objLLH, this);
    methodmap["LLHdim"]=boost::bind(&DEstimator_Iface::objLLHdim, this);
    staticmethodmap["LLH_laplace1D"]=boost::bind(&DEstimator_Iface::objStaticLLH, this, "laplace");
    staticmethodmap["LLH_recursive1D"]=boost::bind(&DEstimator_Iface::objStaticLLH, this, "recursive");
    staticmethodmap["LLH_markov1D"]=boost::bind(&DEstimator_Iface::objStaticLLH, this, "markov");
}

void DEstimator_Iface::objConstruct() 
{
    // [in] Obs: N x Ndim double.  Vector of observed positions
    // [in] T: N x 1 double.  Vector of observation times
    // [in] SE: N x Ndim double.  Vector of observation standard errors
    // [in] exposureT: scalar double.  The duration over which an observation is made
    checkNumArgs(1,4);
    auto Obs=getMat<FloatT>();
    auto T=getVec<FloatT>();
    auto SE=getMat<FloatT>();
    auto exposureT = static_cast<FloatT>(getDouble());
    auto *boxxer = new DEstimator<FloatT>(Obs, T, SE, exposureT);
    outputMXArray(Handle<DEstimator<FloatT>>::makeHandle(boxxer));
}

void DEstimator_Iface::objDestroy()
{
    checkNumArgs(0,1);
    Handle<DEstimator<FloatT>>::destroyObject(rhs[0]);
}

void DEstimator_Iface::getObjectFromHandle(const mxArray *mxhandle)
{
    obj=Handle<DEstimator<FloatT>>::getObject(mxhandle);
}

void DEstimator_Iface::objLLH()
{
    // [in] D: double vector of diffusion constants to estimate LLH for
    // [out] llh: double vector of log-likelihood for each given D value.
    checkNumArgs(1,1);
    auto D=getVec<FloatT>();
    auto llh=makeVec<FloatT>(D.n_elem);
    if(D.n_elem==1) {
        llh(0)=obj->LLH(D(0));
    } else {
        obj->LLH(D, llh); //parallelized
    }
}

void DEstimator_Iface::objLLHdim()
{
    // [in] D: double vector of diffusion constants to estimate LLH for
    // [in] dim: scalar integer giving index of dimension to estimate LLH for.  dim is zero-indexed:  0<=dim<Ndim;
    // [out] llh: double vector of log-likelihood for each given D value only considering dimension dim
    checkNumArgs(1,2);
    auto D=getVec<FloatT>();
    auto dim=getInt();
    auto llh=makeVec<FloatT>(D.n_elem);
    if(D.n_elem==1) {
        llh(0)=obj->LLHdim(D(0),dim);
    } else {
        obj->LLHdim(D, dim, llh); //parallelized
    }
}


void DEstimator_Iface::objStaticLLH(const std::string &method)
{
    // This method will be exposed as static methods LLH_laplace1D, LLH_recursive1D, LLH_markov1D,
    // Each uses the same arguments but a different underlying algorithm
    // [in] D: double vector of diffusion constants to estimate LLH for
    // [in] Obs: N x 1 double.  Vector of observed 1D positions
    // [in] T: N x 1 double.  Vector of observation times
    // [in] SE: N x 1 double.  Vector of 1D observation standard errors
    // [in] exposureT: scalar double.  The duration over which an observation is made
    // [out] llh: double vector of log-likelihood for each given D value.
    checkNumArgs(1,5);
    auto D = getVec<FloatT>();
    auto Obs = getVec<FloatT>();
    auto T = getVec<FloatT>();
    auto SE = getVec<FloatT>();
    auto exposureT = static_cast<FloatT>(getDouble());
    auto llh = makeVec<FloatT>(D.n_elem);
    if(!method.compare("laplace")){
        DEstimator<FloatT>::LLH_laplace1D(D, Obs, T, SE, exposureT, llh);
    } else if(!method.compare("recursive")){
        DEstimator<FloatT>::LLH_recursive1D(D, Obs, T, SE, exposureT, llh);
    } else if(!method.compare("markov")){
        DEstimator<FloatT>::LLH_markov1D(D, Obs, T, SE, exposureT, llh);
    }
}
