/** @file unscentedtransform.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief 
 */


#include "unscentedtransform.h"
#include <stdexcept>


UnscentedTransform::UnscentedTransform(double _kappa, double _alpha, double _beta) 
{
    setParams(_kappa,_alpha,_beta);
}

UnscentedTransform::UnscentedTransform(double _kappa, double _alpha, double _beta, IdxT _Ndim) 
{
    setParams(_kappa,_alpha,_beta);
    setNdim(_Ndim);
}

UnscentedTransform::UnscentedTransform(IdxT _Ndim) 
{
    setNdim(_Ndim);
}

void UnscentedTransform::setParams( double _kappa, double _alpha, double _beta)
{
    if(_kappa<0) throw std::range_error("kappa should be non-negative.");
    if(_alpha<0  || 1<_alpha) throw std::range_error("alpha should be in closed interval [0,1]");
    if(_beta<0) throw std::range_error("beta should be non-negative.");
    kappa = _kappa;
    alpha = _alpha;
    beta = _beta;
}

void UnscentedTransform::setNdim(IdxT _Ndim)
{
    Ndim = _Ndim;
    if(Ndim>0) {
        lambda = alpha*alpha*(Ndim+kappa)-Ndim;
        centerMeanWeight = lambda/(Ndim+lambda);
        double centerCovWeight = centerMeanWeight + 1 + beta - alpha*alpha;
        generalWeight = 1./(2*(Ndim+lambda));
        covWeights = VecT(1+2*Ndim);
        covWeights.fill(generalWeight);
        covWeights(0) = centerCovWeight;
    }
}
