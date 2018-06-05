/** @file unscentedtransform.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief The class declaration and inline and templated functions UnscentedTransform.
 */
#ifndef _UNSCENTEDTRANSFORM_H
#define _UNSCENTEDTRANSFORM_H
#include <armadillo>
#include <iostream>

class UnscentedTransform {
public:
    using IdxT = arma::uword;
    using VecT = arma::vec;
    using MatT = arma::mat;
    
    UnscentedTransform() {};
    UnscentedTransform(double _kappa, double _alpha, double _beta);
    UnscentedTransform(double _kappa, double _alpha, double _beta, IdxT _Ndim);
    UnscentedTransform(IdxT _Ndim);
    
    void setParams( double _kappa, double _alpha, double _beta);
    void setNdim(IdxT _Ndim);

    template<typename FuncT>
    void transform(FuncT &&func, const VecT &mean, const MatT &cov, VecT &tmean, MatT &tcov) const;

private:
    IdxT Ndim=0;
    double kappa = 0;
    double alpha = 0.01;
    double beta = 2;
    
    //Computed from Ndim, kappa, alpha, and beta
    double lambda;
    double centerMeanWeight;
    double generalWeight;
    VecT covWeights;
};

template<typename FuncT>
void UnscentedTransform::transform(FuncT &&func, const VecT &mean, const MatT &cov, VecT &tmean, MatT &tcov) const
{
    if(Ndim<=0 || mean.n_elem!=Ndim) throw std::runtime_error("call setNdim to set number of dimensions");

    auto V = arma::chol((Ndim+lambda)*cov).eval();
    MatT tSigmaPoints(Ndim,2*Ndim+1);
    tSigmaPoints.col(0) = func(mean);
    tmean = centerMeanWeight*tSigmaPoints.col(0);
    for(IdxT n=0;n<Ndim;n++){
        tSigmaPoints.col(n+1) = func(mean+V.col(n));
        tmean += generalWeight * tSigmaPoints.col(n+1);
        tSigmaPoints.col(n+Ndim+1) = func(mean-V.col(n));
        tmean += generalWeight * tSigmaPoints.col(n+Ndim+1);
    }
    tSigmaPoints.each_col() -= tmean;
    tcov = tSigmaPoints*diagmat(covWeights)*tSigmaPoints.t();
 }
 
#endif /* _UNSCENTEDTRANSFORM_H */
