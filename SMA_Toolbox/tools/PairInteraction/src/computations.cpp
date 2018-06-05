/** @file computations.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief Computational helper functions
 */

#include "computations.h"

namespace pair_int {
    
    void gaussianBayes(const VecT &yt_mean, const MatT &yt_cov, 
                       const VecT &xt_mean, const VecT &xt_cov, 
                       VecT &new_mean,      MatT &new_cov)
    {
        auto new_cov_inv = arma::inv(arma::diagmat(xt_cov)) + arma::inv_sympd(yt_cov);
        new_cov = arma::inv(new_cov_inv);
        new_mean = arma::solve(new_cov_inv, arma::solve(arma::diagmat(xt_cov),xt_mean)+arma::solve(yt_cov,yt_mean));
    }

    /* for axis-aligned covariance matricies as vectors */
    void gaussianBayes(const VecT &yt_mean, const VecT &yt_cov, 
                       const VecT &xt_mean, const VecT &xt_cov, 
                             VecT &new_mean,      VecT &new_cov)
    {
        new_cov = 1./(1./xt_cov + 1./yt_cov);
        new_mean = new_cov%(xt_mean/xt_cov + yt_mean/yt_cov);
    }
    
    void transV(double ax, double ay, double bx, double by, 
                double &cx, double &cy, double &r, double &phi)
    {
        double dx = ax-bx;
        double dy = ay-by;
        cx = .5*(ax+bx);
        cy = .5*(ay+by);
        r = sqrt(square(dx) + square(dy));
        phi = fmod(atan2(dy,dx),2*arma::datum::pi);
    }

    arma::vec sampleOU(const arma::vec &times,double mu,double sigma,double gamma, double r0)
    {
        IdxT N = times.n_elem;
        arma::vec samp(N);
        if(times(0)<=0) throw std::invalid_argument("Times must all be greater than 0.");
        double E = exp(-gamma*times(0));
        double muR = mu*(1-E) + r0*E;
        double sigmaR = sigma*sqrt(1-square(E));
        samp(0) = muR + sigmaR*RNG.randn();
        for(IdxT n=1;n<N;n++) {
            double dt = times(n)-times(n-1);
            E = exp(-gamma*dt);
            muR = mu*(1-E) + samp(n-1)*E;
            sigmaR = sigma*sqrt(1-square(E));
            samp(n) = muR + sigmaR*RNG.randn();
        }
        return samp;
    }
    
    arma::vec sampleOURegular(IdxT N, double dt,double mu,double sigma,double gamma,double r0)
    {
        arma::vec samp(N);
        double E=exp(-gamma*dt);
        double sigmaR= sigma*sqrt(1-square(E));
        samp(0) = mu*(1-E) + r0*E + sigmaR*RNG.randn();
        for(IdxT n=1;n<N;n++) samp(n) = mu*(1-E) + samp(n-1)*E + sigmaR*RNG.randn();
        return samp;
    }
    
    arma::vec sampleOUEquilibrium(const arma::vec &times,double mu,double sigma,double gamma)
    {
        IdxT N = times.n_elem;
        arma::vec samp(N);
        if(times(0)<=0) throw std::invalid_argument("Times must all be greater than 0.");
        double E,muR, sigmaR;
        samp(0) = mu + sigma*RNG.randn();
        for(IdxT n=1;n<N;n++) {
            double dt = times(n)-times(n-1);
            E = exp(-gamma*dt);
            muR = mu*(1-E) + samp(n-1)*E;
            sigmaR = sigma*sqrt(1-square(E));
            samp(n) = muR + sigmaR*RNG.randn();
        }
        return samp;
    }
    
    arma::vec sampleOUEquilibriumRegular(IdxT N, double dt,double mu,double sigma,double gamma)
    {
        arma::vec samp(N);
        samp(0) = mu + sigma*RNG.randn();
        double E=exp(-gamma*dt);
        double sigmaR= sigma*sqrt(1-square(E));
        for(IdxT n=1;n<N;n++) samp(n) = mu*(1-E) + samp(n-1)*E + sigmaR*RNG.randn();
        return samp;
    }
    
    arma::vec logNormalize(const arma::vec &log_weights)
    {
        arma::vec ws = log_weights - arma::max(log_weights);
        const double MIN_LOG = -1.0e2; //Smaller than this is treated as
        double sum=0.;
        for(arma::uword n=0; n<ws.n_elem; n++) {
            if(ws(n)>MIN_LOG) {
                ws(n) = exp(ws(n));
                sum+=ws(n);
            } else {
                ws(n)=0;
            }
        }
        ws/=sum;  //Normalize
        return ws;
    }
    
    VecT logNormalizeBinary(const VecT &logS0, const VecT &logS1)
    {
        return 1./(1+exp(logS1-logS0));
    }
    
    /**
     * Compute log(sum(ws)) from vector of log(ws)
     * 
     * 
     */
    double logSum(const arma::vec &ws)
    {
        double C = arma::max(ws);
        const double MIN_LOG = -1.0e2; //Smaller than this is treated as
        double sum=0.;
        for(arma::uword n=0; n<ws.n_elem; n++) if(ws(n)-C > MIN_LOG)  sum += exp(ws(n)-C);
        return C + log(sum);
    }
    

}
