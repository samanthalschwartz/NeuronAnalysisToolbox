/**
 * @file DCompLib.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-21-2014
 * @brief Some fast routines for diffusion estimation computations
 */
#include "DCompLib.h"

template<class FloatT>
FloatT logprod(const arma::Col<FloatT> &X)
{
    static const FloatT max=std::numeric_limits<FloatT>::max();
    static const FloatT min=std::numeric_limits<FloatT>::min();
    FloatT prod=1;
    FloatT Lprod=0;
    for(unsigned n=0; n<X.n_elem; n++) {
        FloatT x=fabs(X(n));
        FloatT temp_prod=prod*x;
        if(temp_prod<min || temp_prod>max) {
            Lprod+=log(prod);
            prod=x;
        } else {
            prod=temp_prod;
        }
    }
    Lprod+=log(prod);
    return Lprod;
}

template<class FloatT>
FloatT SymTriDiag<FloatT>::logdet() const
{
    VecT det_factors(N);
    FloatT B=b(0);
    FloatT D=-a(0)/B;
    FloatT Sdet=sgn(B);
    det_factors(0)=fabs(B);
    for(int n=1;n<N-1;n++){
        B=b(n)+a(n-1)*D;
        D=-a(n)/B;
        Sdet*=sgn(B);
        det_factors(n)=fabs(B);
    }
    B=b(N-1)+a(N-2)*D;
    Sdet*=sgn(B);
    det_factors(N-1)=fabs(B);
    return Sdet*logprod(det_factors);
}

template<class FloatT>
void SymTriDiag<FloatT>::solve(const VecT &Y, VecT &X) const
{
    VecT D(N-1), E(N-1);
    D(0)=-a(0)/b(0);
    E(0)=Y(0)/b(0);
    for(int n=1;n<N-1;n++){
        FloatT B=b(n)+a(n-1)*D(n-1);
        D(n)=-a(n)/B;
        E(n)=(Y(n)-a(n-1)*E(n-1))/B;
    }
    X(N-1)=(Y(N-1)-a(N-2)*E(N-2))/(b(N-1)+a(N-2)*D(N-2));
    for(int n=N-2;n>=0;n--) X(n)=E(n)+D(n)*X(n+1);
}

/* Explicit Template Intatiations */
template class SymTriDiag<float>;
template class SymTriDiag<double>;
template float logprod<float>(const arma::Col<float> &X);
template double logprod<double>(const arma::Col<double> &X);
