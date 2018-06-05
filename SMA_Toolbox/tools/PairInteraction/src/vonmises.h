/** @file vonmises.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief vonMises Distribution class for estimating von Mises distribution
 *        parameters for a difference of gaussians.  Also includes static methods for LLH and other computations
 */

#ifndef _VONMISES_H
#define _VONMISES_H

#include <armadillo>

#include "computations.h"

namespace pair_int {
    
    double vonMisesLLH(double d,double kappa);
    
    double sampleVonMises(double mu, double kappa);

    VecT sampleVonMises(double mu, double kappa, IdxT N);
    
    double estimateVonMisesKappa(const arma::vec &alphas);
    
    void estimateGaussianAngleVonMisesParams(IdxT N, const arma::vec &center, const arma::mat &cholCov,
                                                    double &mu, double &kappa);
     
    inline
    double vonMisesLLH(double d, double kappa)
    {
        return cos(d)*kappa - LOG2PI - Bessel::BE::logI0(kappa);
    }

    
} /* namespace pair_int */
#endif /* _VONMISES_H */
