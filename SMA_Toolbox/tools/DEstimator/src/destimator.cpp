/**
 * @file destimator.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @author Peter Relich (physx.grad\@gmail.com)
 * @date 08-05-2014
 * @brief The class definitions for DEstimator.
 */
#include <exception>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdio>

#include "destimator.h"
#include "DCompLib.h"

/* Static member variables */
template<class FloatT>
const FloatT DEstimator<FloatT>::log2pi=log(2*arma::Datum<FloatT>::pi);

template<class FloatT>
const FloatT DEstimator<FloatT>::min_variance=std::numeric_limits<FloatT>::epsilon();

template<class FloatT>
const FloatT DEstimator<FloatT>::min_laplace_variance=1e-8;

/* Constructor */
template<class FloatT>
DEstimator<FloatT>::DEstimator(const MatT &Obs_, const VecT &T_, const MatT &SE_, FloatT exposureT_)
    : Ndim(static_cast<int>(Obs_.n_cols)), N(static_cast<int>(Obs_.n_rows)),
      Obs(Obs_), T(T_), SE(SE_), exposureT(exposureT_)
{
    if(Obs.n_rows!=T.n_elem) throw std::length_error("Obs and T length missmatch.");
    if(Obs.n_rows!=SE.n_rows) throw std::length_error("Obs and SE length missmatch.");
    if(Obs.n_cols!=SE.n_cols) throw std::length_error("Obs and SE dim missmatch.");
    dT = vec_delta(T);
    vM.set_size(N);
    vD.set_size(N-1);
}

/* Public Methods */
template<class FloatT>
FloatT DEstimator<FloatT>::LLH(FloatT D)
{
    FloatT llh=0;
    for(int dim=0; dim<Ndim; dim++){
        FloatT *pObs=Obs.memptr()+dim*N; //Select out the "dim" column from the matrix of observations
        computeVariance(D,dT,SE.col(dim),exposureT,min_variance,vD,vM);
        llh += coreLLH(N,pObs,dT.memptr(),vD.memptr(),vM.memptr());
    }
    return llh;
}

template<class FloatT>
FloatT DEstimator<FloatT>::LLHdim(FloatT D, int dim)
{
    if(dim<0 || dim>=Ndim) throw std::domain_error("dim>=Ndim || dim<0");
    computeVariance(D,dT,SE.col(dim),exposureT,min_variance,vD,vM);
    FloatT *pObs = Obs.memptr()+N*dim; //Select out the "dim" column from the matrix of observations
    return coreLLH(N,pObs,dT.memptr(),vD.memptr(),vM.memptr());
}

template<class FloatT>
void DEstimator<FloatT>::LLH(const VecT &D, VecT &LLH)
{
    int ND = static_cast<int>(D.n_elem);
    LLH.zeros();
    #pragma omp parallel
    {
        VecT local_vD(N-1);
        VecT local_vM(N);
        #pragma omp for
        for(int k=0; k<ND; k++){
            for(int dim=0; dim<Ndim; dim++) {
                FloatT *pObs = Obs.memptr()+dim*N; //Select out the "dim" column from the matrix of observations
                computeVariance(D(k),dT,SE.col(dim),exposureT,min_variance,local_vD,local_vM);
                LLH(k) += coreLLH(N,pObs,dT.memptr(),local_vD.memptr(),local_vM.memptr());
            }
        }
    }
}

template<class FloatT>
void DEstimator<FloatT>::LLHdim(const VecT &D, int dim, VecT &LLH)
{
    int ND = static_cast<int>(D.n_elem);
    if(dim<0 || dim>=Ndim) throw std::domain_error("dim>=Ndim || dim<0");
    FloatT *pObs=Obs.memptr()+N*dim; //Select out the "dim" column from the matrix of observations
    LLH.zeros();
    #pragma omp parallel
    {
        VecT local_vD(N-1);
        VecT local_vM(N);
        #pragma omp for
        for(int k=0; k<ND; k++){
            computeVariance(D(k),dT,SE.col(dim),exposureT,min_variance,local_vD,local_vM);
            LLH(k) = coreLLH(N,pObs,dT.memptr(),local_vD.memptr(),local_vM.memptr());
        }
    }
}

/* Public Static Methods */
template<class FloatT>
void DEstimator<FloatT>::LLH_recursive1D(const VecT &D, const VecT &Obs, 
                        const VecT &T, const VecT &SE, FloatT exposureT,  VecT &LLH)
{
    int N = static_cast<int>(Obs.n_elem);
    int ND = static_cast<int>(D.n_elem);
    checkLLHargs(D,Obs,T,SE,LLH);
    VecT dT = vec_delta(T);
    #pragma omp parallel
    {
        VecT vD(N-1);
        VecT vM(N);
        #pragma omp for
        for(int k=0;k<ND;k++){
            computeVariance(D(k),dT,SE,exposureT,min_variance,vD,vM);
            LLH(k) = coreLLH(N,Obs.memptr(),dT.memptr(),vD.memptr(),vM.memptr());
        }
    }
}

template<class FloatT>
void DEstimator<FloatT>::LLH_laplace1D( const VecT &D, const VecT &Obs, const VecT &T, const VecT &SE, FloatT exposureT, VecT &LLH)
{
    int N = static_cast<int>(Obs.n_elem);
    int ND = static_cast<int>(D.n_elem);
    checkLLHargs(D,Obs,T,SE,LLH);
    
    VecT dT = vec_delta(T);
    LLH.zeros();
    #pragma omp parallel
    {
        VecT vD(N-1), vM(N);
        VecT ivD(N-1), ivM(N);
        VecT theta(N);
        SymTriDiag<FloatT> hess(N);
        #pragma omp for
        for(int k=0;k<ND;k++){
            computeVariance(D(k),dT,SE,exposureT,min_laplace_variance,vD,vM);
            int vMsign = 1;
            for(int n=0;n<N;n++) vMsign *= sgn(vM(n));
            //Compute hessian
            ivD = 1./vD;
            ivM = 1./vM;
            hess.a = -ivD;
            hess.b(0) = ivM(0)+ivD(0);
            for(int n=1; n<N-1; n++) hess.b(n) = ivM(n)+ivD(n)+ivD(n-1);
            hess.b(N-1) = ivM(N-1)+ivD(N-2);
            hess.solve(Obs%ivM, theta);//Solve for most likely true particle positions, theta

            //Compute log likelihood
            for(int n=0;n<N;n++) LLH(k) += square(Obs(n)-theta(n))*ivM(n);
            for(int n=0;n<N-1;n++) LLH(k) += square(theta(n+1)-theta(n))*ivD(n);
            LLH(k) += (N-1)*log2pi;
            LLH(k) += logprod(vM);
            LLH(k) += logprod(vD);
            LLH(k) += vMsign*hess.logdet();
            LLH(k) *= -0.5;
        }
    }
}


template<class FloatT>
void DEstimator<FloatT>::LLH_markov1D( const VecT &D, const VecT &Obs, const VecT &T, const VecT &SE, FloatT exposureT, VecT &LLH)
{
    int N = static_cast<int>(Obs.n_elem);
    int ND = static_cast<int>(D.n_elem);
    checkLLHargs(D,Obs,T,SE,LLH);
    VecT dT = vec_delta(T);
    LLH.zeros();

    #pragma omp parallel
    {
        VecT vD(N-1), vM(N), theta(N-1);
        VecT dObs = vec_delta(Obs);
        SymTriDiag<FloatT> cov(N-1);
        #pragma omp for
        for(int k=0;k<ND;k++){
            computeVariance(D(k),dT,SE,exposureT,min_variance,vD,vM);
            //Compute covariance tri-diagonal symmetric matrix
            cov.a = -vM.subvec(1,N-1); //Off diagonal
            for(int n=0;n<N-1;n++) cov.b(n) = vD(n)+vM(n)+vM(n+1); //Central diagonal
            cov.solve(dObs, theta);
            LLH(k) += (N-1)*log2pi;
            LLH(k) += arma::dot(dObs, theta);//problem in mkl avx code for 2013b(?)
            LLH(k) += cov.logdet();
            LLH(k) *= -0.5;
        }
    }
}

/* Private Static Methods */
template<class FloatT>
FloatT DEstimator<FloatT>::coreLLH(int N, const FloatT Obs[], const FloatT dT[], const FloatT vD[], const FloatT vM[])
{
    VecT alpha(N-1);
    FloatT eta = vD[0]+vM[0];
    FloatT mu = Obs[0];
    FloatT LLH = 0;
    for(int n=1;n<N-1;n++){
        FloatT temp_alpha = vM[n]+eta;
        alpha[n-1] = temp_alpha;
        LLH += square(Obs[n]-mu)/temp_alpha;
        mu = (mu*vM[n]+Obs[n]*eta)/temp_alpha;
        eta = vM[n]*eta/temp_alpha+vD[n];
    }
    alpha[N-2] = vM[N-1]+eta;
    LLH += (N-1)*log2pi;
    LLH += logprod(alpha);
    LLH += square(Obs[N-1]-mu)/alpha[N-2];
    LLH *= -0.5; //Factored out -0.5 for everything above
    return LLH;
}

template<class FloatT>
void DEstimator<FloatT>::computeVariance(FloatT D, const VecT &dT, const VecT &SE,
                                         FloatT exposureT ,FloatT min_variance, VecT &vD, VecT &vM)
{
    int N = static_cast<int>(SE.n_elem);
    D=fabs(D);
    vD = 2*D*dT;
    for(int n=0; n<N-1; n++) if(vD(n)<min_variance) vD(n)=min_variance;
    vM = SE%SE - 2*D*exposureT/6; //variance is square of SE
    for(int n=0; n<N; n++) {
        if(fabs(vM(n))<min_variance)
            vM(n) = sgn(vM(n))<0 ? -min_variance : min_variance;
    }
}

template<class FloatT>
inline void DEstimator<FloatT>::checkLLHargs(const VecT &D, const VecT &Obs, const VecT &T, const VecT &SE,  VecT &LLH)
{
    if(D.n_elem<1) throw std::length_error("len(D)<1.");
    if(Obs.n_elem<2) throw std::length_error("len(Obs)<2.");
    if(Obs.n_rows!=T.n_elem) throw std::length_error("Obs and T length missmatch.");
    if(Obs.n_elem!=SE.n_elem) throw std::length_error("Obs and SE length missmatch.");
    if(D.n_elem!=LLH.n_elem) throw std::length_error("D and LLH length missmatch.");
}

/* Explicit Template Instantiation */
/* These ensure the compiler emits code for both the double and float versions of DEstimator */
template class DEstimator<float>;
template class DEstimator<double>;
