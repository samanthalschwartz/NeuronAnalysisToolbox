/** @file laplacetform.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief The class declaration and inline for the numerical laplace transform computations
 *
 * Numerical Laplace inverse routines.
 * LaplaceInverseGS - GaverStehfast routine.  Parrallel with openmp.
 * 
 */

#ifndef _LAPLACETFORM_H
#define _LAPLACETFORM_H

#include <armadillo>
#include <limits> 

class LaplaceInverseGS{
public:
    using VecT = arma::vec;
 
    static const int MinOrder=2;
    static const int MaxOrder=14;
    static const int DefaultOrder=7;
    int order;
    VecT weights;

    LaplaceInverseGS(int order=DefaultOrder);
    
    template<typename Func>
    double invert(Func &&F, double t) const;
    
    template<typename Func>
    double invertParallel(Func &&F, double t) const;
    
    void computeWeights(int order, VecT &weights);

private:
    static constexpr double log2=log(2);
    double factorial(double n) const;
};


template<typename Func>
double LaplaceInverseGS::invert(Func &&F, double t) const
{
    VecT vals(2*order);
    double alpha = log2/t;
    for(int k=1; k<=2*order; k++){
        vals(k-1) = weights(k-1)*F(k*alpha);
    }
    return alpha*arma::sum(vals);
}


template<typename Func>
double LaplaceInverseGS::invertParallel(Func &&F, double t) const
{
    VecT vals(2*order);
    double alpha = log2/t;
    #pragma omp parallel for
    for(int k=1; k<=2*order; k++){
        vals(k-1) = weights(k-1)*F(k*alpha);
    }
    return alpha*arma::sum(vals);
}

inline double LaplaceInverseGS::factorial(double n) const {
    if(n>=18) {
        return tgamma(n+1);
    } else {
        double val=1;
        for(int k=2;k<=n;k++) val*=k;
        return val;
    }
}


#endif /* _LAPLACETFORM_H */
