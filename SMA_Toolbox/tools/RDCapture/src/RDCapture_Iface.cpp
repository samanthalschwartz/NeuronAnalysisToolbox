/** @file RDCapture_Iface.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief The class declaration and inline and templated functions for RDCapture_Iface.
 */
#include "RDCapture_Iface.h"
#include "Handle.h"
#include "bessel.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    RDCapture_Iface iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}

RDCapture_Iface::RDCapture_Iface() 
    : Mex_Iface("RDCapture_Iface") 
{
    methodmap["survivalProb"]=boost::bind(&RDCapture_Iface::objSurvivalProb, this);
    methodmap["mu"]=boost::bind(&RDCapture_Iface::objMu, this);
    methodmap["nu"]=boost::bind(&RDCapture_Iface::objNu, this);
    methodmap["I0"]=boost::bind(&RDCapture_Iface::objI0, this);
    methodmap["I1"]=boost::bind(&RDCapture_Iface::objI1, this);
    methodmap["K0"]=boost::bind(&RDCapture_Iface::objK0, this);
    methodmap["K1"]=boost::bind(&RDCapture_Iface::objK1, this);
    methodmap["logI0"]=boost::bind(&RDCapture_Iface::objlogI0, this);
    methodmap["logI1"]=boost::bind(&RDCapture_Iface::objlogI1, this);
    methodmap["logK0"]=boost::bind(&RDCapture_Iface::objlogK0, this);
    methodmap["logK1"]=boost::bind(&RDCapture_Iface::objlogK1, this);
    methodmap["I0K0"]=boost::bind(&RDCapture_Iface::objI0K0, this);
    methodmap["I0K1"]=boost::bind(&RDCapture_Iface::objI0K1, this);
    methodmap["I1K0"]=boost::bind(&RDCapture_Iface::objI1K0, this);
    methodmap["I1K1"]=boost::bind(&RDCapture_Iface::objI1K1, this);
    methodmap["simulate"]=boost::bind(&RDCapture_Iface::objSimulate, this);
}

void RDCapture_Iface::objConstruct() 
{
    // [in] D: double.  Diffusion const [um^2/s]
    // [in] rho: double.  Capture Radius [um]
    // [in] lambda: double.  Capture rate [1/s]
    checkNumArgs(1,3);
    auto D=getDouble();
    auto rho=getDouble();
    auto lambda=getDouble();
    auto *obj = new RDCapture(D,rho,lambda);
    outputMXArray(Handle<RDCapture>::makeHandle(obj));
}

void RDCapture_Iface::objDestroy()
{
    checkNumArgs(0,1);
    Handle<RDCapture>::destroyObject(rhs[0]);
}

void RDCapture_Iface::getObjectFromHandle(const mxArray *mxhandle)
{
    obj=Handle<RDCapture>::getObject(mxhandle);
}

void RDCapture_Iface::objSurvivalProb()
{
    // [in] r0: Vector size:[N] or scalar - initial starting positions um.
    // [in] t: Vector size:[N] or scalar - times >0
    // If r0 and t are both vectors, they should be the same length and r0(n) is computed at time t(n).
    // [out] Q: Vector size:[N] Survival probability vector is of size N.
    checkNumArgs(1,2);
    auto r0=getVec<FloatT>();
    auto t=getVec<FloatT>();
    if(r0.n_elem>1 && t.n_elem>1 && r0.n_elem != t.n_elem) {
        mexErrMsgIdAndTxt("RDCapture_Iface:BadParameters","r0 and t must have the same length if both are vectors");
    }
    auto Q=makeVec<FloatT>(std::max(r0.n_elem,t.n_elem));
    if (t.n_elem ==1){
        if(r0.n_elem==1) Q(0) = obj->survivalProb(r0(0), t(0));
        else obj->survivalProb_parallel(Q, r0, t(0));
    } else {
        if(r0.n_elem==1) obj->survivalProb_parallel(Q, r0(0), t);
        else obj->survivalProb_parallel(Q, r0, t);
    }
}

void RDCapture_Iface::objMu()
{
    // [in] r0: Vector size:[N] or scalar - initial starting positions um.
    // [in] t: Vector size:[N] or scalar - times >0
    // If r0 and t are both vectors, they should be the same length and r0(n) is computed at time t(n).
    // [out] mu: Vector size:[N] mu computation.
    checkNumArgs(1,2);
    auto r0=getVec<FloatT>();
    auto t=getVec<FloatT>();
    if(r0.n_elem>1 && t.n_elem>1 && r0.n_elem != t.n_elem) {
        mexErrMsgIdAndTxt("RDCapture_Iface:BadParameters","r0 and t must have the same length if both are vectors");
    }
    auto mu=makeVec<FloatT>(std::max(r0.n_elem,t.n_elem));
    if (t.n_elem ==1){
        if(r0.n_elem==1) mu(0) = obj->mu(r0(0), t(0));
        else obj->mu_parallel(mu, r0, t(0));
    } else {
        if(r0.n_elem==1) obj->mu_parallel(mu, r0(0), t);
        else obj->mu_parallel(mu, r0, t);
    }
}

void RDCapture_Iface::objNu()
{
    // [in] t: Vector size:[N] or scalar - times >0
    // [out] nu: Vector size:[N] nu computation.
    checkNumArgs(1,1);
    auto t=getVec<FloatT>();
    auto nu=makeVec<FloatT>(t.n_elem);
    if (t.n_elem ==1) nu(0) = obj->nu(t(0));
    else obj->nu_parallel(nu, t);
}

void RDCapture_Iface::objSimulate()
{
    checkNumArgs(1,4);
    
    auto r0 = getVec<double>();
    auto ts = getVec<double>();
    IdxT N = getScalar<IdxT>();
    double max_dt = getScalar<double>();
    IdxT Nts = ts.n_elem;
    IdxT Nr0 = r0.n_elem;
    auto p = makeMat<double>(Nts,Nr0);
    obj->simulate(r0,ts,N,max_dt,p);
}


void RDCapture_Iface::objI0()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I0(x(n));
}

void RDCapture_Iface::objI1()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I1(x(n));
}

void RDCapture_Iface::objK0()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::K0(x(n));
}

void RDCapture_Iface::objK1()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::K1(x(n));
}

void RDCapture_Iface::objlogI0()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::logI0(x(n));
}

void RDCapture_Iface::objlogI1()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::logI1(x(n));
}

void RDCapture_Iface::objlogK0()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::logK0(x(n));
}

void RDCapture_Iface::objlogK1()
{
    checkNumArgs(1,1);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::logK1(x(n));
}

void RDCapture_Iface::objI0K0()
{
    checkMinNumArgs(1,1);
    checkMaxNumArgs(1,2);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    if(nrhs==2) {
        auto y = getVec<double>();
        if(y.n_elem != x.n_elem) throw std::invalid_argument("x and y must have the same length");
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I0K0(x(n),y(n));
    } else {
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I0K0(x(n));
    }
}


void RDCapture_Iface::objI0K1()
{
    checkMinNumArgs(1,1);
    checkMaxNumArgs(1,2);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    if(nrhs==2) {
        auto y = getVec<double>();
        if(y.n_elem != x.n_elem) throw std::invalid_argument("x and y must have the same length");
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I0K1(x(n),y(n));
    } else {
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I0K1(x(n));
    }
}

void RDCapture_Iface::objI1K0()
{
    checkMinNumArgs(1,1);
    checkMaxNumArgs(1,2);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    if(nrhs==2) {
        auto y = getVec<double>();
        if(y.n_elem != x.n_elem) throw std::invalid_argument("x and y must have the same length");
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I1K0(x(n),y(n));
    } else {
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I1K0(x(n));
    }
}


void RDCapture_Iface::objI1K1()
{
    checkMinNumArgs(1,1);
    checkMaxNumArgs(1,2);
    auto x = getVec<double>();
    arma::uword N = x.n_elem;
    auto val = makeDVec(N);
    if(nrhs==2) {
        auto y = getVec<double>();
        if(y.n_elem != x.n_elem) throw std::invalid_argument("x and y must have the same length");
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I1K1(x(n),y(n));
    } else {
        for(arma::uword n=0; n<N; n++) val(n)=Bessel::BE::I1K1(x(n));
    }
}




