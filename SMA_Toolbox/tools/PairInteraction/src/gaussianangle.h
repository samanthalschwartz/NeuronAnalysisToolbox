/** @file gaussianangle.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 02-2017
 * @brief class for computation of p(phi | mu, sigma).  The pdf of the angle to the origin of a point
 *        drawn from 2D multivariate normal (bivariate normal) distribution.
 */

#ifndef _GAUSSIANANGLE_H
#define _GAUSSIANANGLE_H

#include "computations.h"

namespace pair_int {

    class GaussianAngleDist 
    {
    public:
        using VecT=arma::vec;
        using MatT=arma::mat;
        GaussianAngleDist() {};
        GaussianAngleDist(const VecT &mu_, const MatT &sigma_);
        GaussianAngleDist(const VecT &mu_, double sigmaX_, double sigmaY_, double rho_);
        GaussianAngleDist(const VecT &mu_, double sigmaX_, double sigmaY_);
        void setSigma(const MatT &sigma_);
        void setSigma(double sigmaX_, double sigmaY_, double rho_);
        void setSigma(double sigmaX_, double sigmaY_);
        void setMean(const VecT &mu);
        
        double computePDF(double theta) const;
        double sample() const;
        
    protected:
        bool sigma_is_set=false;
        double rho,sigmaX,sigmaY;
        double det_sigma, sqrt_det_sigma;
        double sigma00, sigma11, sigma01; // elements of covarance matrix [sigma00, sigma01;sigma10,sigma11]
        bool mu_is_set=false;
        VecT mu;
        double pdf_mu;
        static const double constexpr sqrt2 = sqrt(2);
        static const double constexpr sqrt2pi = sqrt(2*arma::datum::pi);
    };
    
} /* pair_int */

#endif /* _GAUSSIANANGLE_H */
