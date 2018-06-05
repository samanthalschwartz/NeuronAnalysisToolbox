/** @file computations.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief Computational helper functions
 */

#ifndef _COMPUTATIONS_H
#define _COMPUTATIONS_H

#include <armadillo>
#include <cmath>
#include <complex>
#include "rngmanager.h"

#include "bessel.h"

namespace pair_int {
    static constexpr const double LOG2PI = log(2*arma::datum::pi);
    using IdxT = arma::uword;
    using VecT = arma::vec;
    using MatT = arma::mat;
    
    template<typename NumT>
    NumT square(NumT x);
    
    double gaussianLLH(double dist,double var);

    double gaussianLLH(const VecT &dist,const VecT &var);
    double gaussianLLH(const VecT &x, const VecT &mean, const VecT &var);
    
    VecT gaussianLLH_all(const VecT &dist,const VecT &var);
    
    double multiGaussianLLH(const VecT &d, const MatT &sigma);
    double multiGaussianLLH(const VecT &p, const VecT &mean, const MatT &sigma);
    
    VecT sampleGaussian(const VecT &mean, const VecT &sigma);
    VecT sampleGaussianFull(const VecT &mean, const MatT &sigma_chol);
    
    double sampleGaussian(double mean, double sigma);
    
    void gaussianBayes(const VecT &yt_mean, const MatT &yt_cov, 
                       const VecT &xt_mean, const VecT &xt_cov, 
                             VecT &new_mean,      MatT &new_cov);
    void gaussianBayes(const VecT &yt_mean, const VecT &yt_cov, 
                       const VecT &xt_mean, const VecT &xt_cov, 
                             VecT &new_mean,      VecT &new_cov);
    

    
    
    double distance(double ax, double ay, double bx, double by);
    
    double distance(const VecT &y);
    
    void transV(double ax, double ay, double bx, double by, 
                double &cx, double &cy, double &r, double &phi);
    VecT transV(const VecT &y);
    VecT transVinv(const VecT &v);
    
    void sampleOU(double t,double mu,double sigma,double gamma,double r0, double &r, double &llh);
    void sampleOUEquilibrium(double t,double mu,double sigma, double &r, double &llh);
    arma::vec sampleOU(const arma::vec &times,double mu,double sigma,double gamma,double r0);
    arma::vec sampleOURegular(IdxT N, double dt,double mu,double sigma,double gamma,double r0);
    arma::vec sampleOUEquilibrium(const arma::vec &times,double mu,double sigma,double gamma);
    arma::vec sampleOUEquilibriumRegular(IdxT N, double dt,double mu,double sigma,double gamma);
    
    arma::vec logNormalize(const arma::vec &log_weights);
    VecT logNormalizeBinary(const VecT &logS0, const VecT &logS1);
    double logSum(const arma::vec &log_weights);
    
    /* Inlined Methods */
    template<typename NumT>
    inline 
    NumT square(NumT x)
    {
        return x*x;
    }
    
    inline 
    double gaussianLLH(double dist,double var)
    {
        return -.5*(LOG2PI + log(var) + dist*dist/var);
    }
    
    inline 
    double gaussianLLH(const VecT &dist,const VecT &var)
    {
        double llh=0;
        for(IdxT n=0; n<dist.n_elem; n++) llh+=gaussianLLH(dist(n),var(n));
        return llh;
    }
    
    inline 
    double gaussianLLH(const VecT &x, const VecT &mean, const VecT &var)
    {
        double llh=0;
        for(IdxT n=0; n<x.n_elem; n++) llh+=gaussianLLH(x(n)-mean(n),var(n));
        return llh;
    }
    
    inline
    VecT gaussianLLH_all(const VecT &dist,const VecT &var)
    {
        IdxT N = dist.n_elem;
        arma::vec llhAll(N);
        for(IdxT n=0; n<dist.n_elem; n++) llhAll(n) = gaussianLLH(dist(n),var(n));
        return llhAll;
    }
    
    inline
    double multiGaussianLLH(const VecT &d, const MatT &sigma)
    {
        arma::uword N = d.n_elem;
        double llh = LOG2PI*N + log(arma::det(sigma)) + arma::sum(arma::dot(d,inv(sigma)*d));
        llh *=-.5;
        return llh;
    }
    
    inline
    double multiGaussianLLH(const VecT &p, const VecT &mean, const MatT &sigma)
    {
        arma::uword N = p.n_elem;
        auto d = p-mean;
        double llh = LOG2PI*N + log(arma::det(sigma)) + arma::sum(arma::dot(d,inv(sigma)*d));
        llh *=-.5;
        return llh;
    }
    
    inline 
    double distance(double ax, double ay, double bx, double by)
    {
        return sqrt(square(ax-bx)+square(ay-by));
    }
    
    inline
    double distance(const VecT &y)
    {
        return sqrt(square(y(0)-y(2))+square(y(1)-y(3)));
    }
    
    inline
    arma::vec transV(const VecT &y)
    {
        double dx = y(0)-y(2); //a_x-b_x
        double dy = y(1)-y(3); //a_y-b_y
        double cx = .5*(y(0)+y(2));
        double cy = .5*(y(1)+y(3));
        double r = sqrt(square(dx) + square(dy));
        double phi = fmod(atan2(dy,dx),2*arma::datum::pi); //Angle of \vec{a} - \vec{b} normalized to [0,2*pi)
        arma::vec tmp{cx,cy,r,phi};
        return tmp;
    }
    
    inline
    VecT transVinv(const VecT &v)
    {
        double cosphi = cos(v(3));
        double sinphi = sin(v(3));
        VecT tmp{v(0)+.5*v(2)*cosphi, v(1)+.5*v(2)*sinphi,
                  v(0)-.5*v(2)*cosphi, v(1)-.5*v(2)*sinphi};
        return tmp;
    }
    
    inline 
    VecT sampleGaussian(const VecT &mean, const VecT &sigma)
    {
        return mean+sigma%RNG.randn(mean.n_elem);
    }
    
    inline 
    VecT sampleGaussianFull(const VecT &mean, const MatT &sigma_chol)
    {
        return mean + sigma_chol*RNG.randn(mean.n_elem);
    }
    
    inline
    double sampleGaussian(double mean, double sigma)
    {
        return mean + sigma*RNG.randn();
    }
    
    
    inline
    void sampleOU(double t,double mu,double sigma,double gamma,double r0, double &r, double &llh)
    {
        double E = exp(-gamma*t);
        double muR = mu*(1-E) + r0*E;
        double varR = sigma*sigma*(1-square(E));
        r = muR + sqrt(varR)*RNG.randn();
        llh = gaussianLLH(r-muR, varR);
    }
    
    inline
    void sampleOUEquilibrium(double t,double mu,double sigma, double &r, double &llh)
    {
        r = mu + sigma*RNG.randn();
        llh = gaussianLLH(r-mu, sigma*sigma);
        
    }
    

} /* namespace pair_int */

#endif /* _COMPUTATIONS_H */
