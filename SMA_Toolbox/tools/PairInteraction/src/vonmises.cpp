/** @file vonmises.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 02-2017
 * @brief vonMises Distribution class for estimating von Mises distribution
 *        parameters for a difference of gaussians.  Also includes static methods for LLH and other computations
 */

#include "vonmises.h"
#include "rngmanager.h"

namespace pair_int {

    double estimateVonMisesKappa(const arma::vec &alphas)
    {
        IdxT N = alphas.n_elem;
        double kappa;
        std::complex<double> i(0,1);
        double R = std::abs(arma::sum(arma::exp(i*alphas)))/double(N);
        if(R<0.53) {
            kappa = 2*R+ pow(R,3) + (5/6)*pow(R,5);
        } else if(R<0.85) {
            kappa = -.4 + 1.39*R+0.43/(1-R);
        } else {
            kappa = 1/(pow(R,3)-4*pow(R,2)+3*R);
        }
        if(N<15) {
            kappa = std::max(kappa-2/(N*kappa), 0.);
        } else {
            kappa = pow(N-1,3)*kappa/(pow(N,3)+N);
        }
        return kappa;
    }
    
    double sampleVonMises(double mu, double kappa)
    {
        if(kappa<1e-6) return RNG.randu()*2*arma::datum::pi; //Uniform distribution
        double alpha = 1 + sqrt(1+4*kappa*kappa);
        double beta = (alpha - sqrt(2*alpha))/(2*kappa);
        double r = (1+beta*beta)/(2*beta);
        double zeta,f,c,u;
        do {
            zeta = cos(arma::datum::pi*RNG.randu());
            f = (1+r*zeta)/(r+zeta);
            c = kappa*(r-f);
            u = RNG.randu();
        } while ((u > c*(2-c)) && (c > log(c)-log(u)+1));
        double v = mu+copysign(1, RNG.randu()-.5)*std::acos(f);
        std::complex<double> vi{0,v};
        return fmod(arg(exp(vi)), 2*arma::datum::pi);
    }
    
    
    VecT sampleVonMises(double mu, double kappa, IdxT N)
    {
        if(kappa<1e-6) return RNG.randu(N)*2*arma::datum::pi; //Uniform distribution
        double alpha = 1 + sqrt(1+4*kappa*kappa);
        VecT samp(N);
        double beta = (alpha - sqrt(2*alpha))/(2*kappa);
        double r = (1+beta*beta)/(2*beta);
        double zeta,f,c,u;
        for(IdxT n=0; n<N; n++) {
            do {
                zeta = cos(arma::datum::pi*RNG.randu());
                f = (1+r*zeta)/(r+zeta);
                c = kappa*(r-f);
                u = RNG.randu();
            } while( (u > c*(2-c)) && (c > log(c)-log(u)+1));
            double v = mu + copysign(1, RNG.randu()-.5)*acos(f);
            samp(n) = fmod(std::arg(std::exp(std::complex<double>{0,v})) ,2*arma::datum::pi);
        }
        return samp;
    }
    
    void estimateGaussianAngleVonMisesParams(IdxT N, const arma::vec &center, const arma::mat &cholCov,
                                                    double &mu, double &kappa)
    {
        arma::mat samp = cholCov*RNG.randn(2,N);
        samp.each_col() += center;
        mu = atan2(center(1),center(0));
        arma::vec angles = arma::atan2(samp.row(1),samp.row(0));
        kappa = estimateVonMisesKappa(angles);
    }
    
} /* namespace pair_int */
