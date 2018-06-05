/** @file gaussianangle.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 02-2017
 * @brief vclass for computation of p(phi | mu, sigma).  The pdf of the angle to the origin of a point
 *        drawn from 2D multivariate normal (bivariate normal) distribution.
 */

#include "gaussianangle.h"

namespace pair_int {

    GaussianAngleDist::GaussianAngleDist(const VecT &mu_, const MatT &sigma_)
    {
        setSigma(sigma_);
        setMean(mu_);
    }
    
    GaussianAngleDist::GaussianAngleDist(const VecT &mu_, double sigmaX_, double sigmaY_, double rho_)
    {
        setSigma(sigmaX_,sigmaY_,rho_);
        setMean(mu_);
    }
    
    GaussianAngleDist::GaussianAngleDist(const VecT &mu_, double sigmaX_, double sigmaY_)
    {
        setSigma(sigmaX_,sigmaY_);
        setMean(mu_);
    }
    
    void GaussianAngleDist::setSigma(const MatT &sigma_)
    {
        if(sigma_.n_rows != 2 || sigma_.n_cols !=2) throw std::logic_error("Covariance matrix must be 2x2");
        sigma00 = sigma_(0,0);
        sigma01 = sigma_(0,1);
        sigma11 = sigma_(1,1);
        sigmaX = sqrt(sigma00);
        sigmaY = sqrt(sigma11);
        rho = sigma01/(sigmaX*sigmaY);
        sqrt_det_sigma = sigmaX*sigmaY*sqrt(1-square(rho));
        det_sigma = sqrt_det_sigma*sqrt_det_sigma;
        sigma_is_set = true;
        if(mu_is_set) setMean(mu); //recompute pdf_mu
    }
    
    void GaussianAngleDist::setSigma(double sigmaX_, double sigmaY_, double rho_)
    {
        if(sigmaX_<= 0 || sigmaY_ <= 0) throw std::logic_error("Sigma must be positive");
        if(rho< -1 ||  rho > 1) throw std::logic_error("rho out of range");
        sigmaX = sigmaX_;
        sigmaY = sigmaY_;
        rho = rho_;
        sigma00 = square(sigmaX);
        sigma01 = rho*sigmaX*sigmaY;
        sigma11 = square(sigmaY);
        sqrt_det_sigma = sigmaX*sigmaY*sqrt(1-square(rho));
        det_sigma = sqrt_det_sigma*sqrt_det_sigma;
        sigma_is_set = true;
        if(mu_is_set) setMean(mu); //recompute pdf_mu
    }

    void GaussianAngleDist::setSigma(double sigmaX_, double sigmaY_)
    {
        if(sigmaX_<= 0 || sigmaY_ <= 0) throw std::logic_error("Sigma must be positive");
        sigmaX = sigmaX_;
        sigmaY = sigmaY_;
        rho = 0;
        sigma00 = square(sigmaX);
        sigma01 = 0;
        sigma11 = square(sigmaY);
        sqrt_det_sigma = sigmaX*sigmaY;
        det_sigma = sqrt_det_sigma*sqrt_det_sigma;
        sigma_is_set = true;
        if(mu_is_set) setMean(mu); //recompute pdf_mu
    }
    
    void GaussianAngleDist::setMean(const VecT &mu_)
    {
        mu = mu_;
        if(sigma_is_set) {
            double H = square(sigmaY*mu(0)) - 2*sigmaX*sigmaY*rho*mu(0)*mu(1) + square(sigmaX*mu(1));
            pdf_mu = exp(-H/(2*det_sigma))/(2*arma::datum::pi*sqrt_det_sigma); //norm2D_PDF(H,mu,sigma);
        }
        mu_is_set = true;
    }
    
    double GaussianAngleDist::computePDF(double theta) const
    {
        if(!mu_is_set) throw std::logic_error("mu is not set");
        if(!sigma_is_set) throw std::logic_error("sigma is not set");
        double costheta = cos(theta);
        double sintheta = sin(theta);
        double C = (sigma11*square(costheta) - sigma01*sin(2*theta) + sigma00*square(sintheta)) / det_sigma;
        double sqrtC = sqrt(C);
        double D = (mu(0)*(sigma11*costheta - sigma01*sintheta) + mu(1)*(sigma11*sintheta - sigma01*costheta)) / (sqrtC*det_sigma);
        double G = (mu(0)*sintheta-mu(1)*costheta)/(sqrtC*sqrt_det_sigma);
        double cdfD = .5*(1+erf(D/sqrt2)); //norm1D_CDF(D;0,1)
        double pdfG = exp(-.5*square(G))/sqrt2pi; //norm1D_PDF(G;0,1)
        double Q = pdf_mu + D*cdfD*pdfG/sqrt_det_sigma;
        return Q/C;
    }
    
    double GaussianAngleDist::sample() const
    {
        VecT p = RNG.randn(2);
        double X = mu(0) + sigmaX*p(0);
        double Y = mu(1) + sigmaY*(p(0)*rho + p(1)*sqrt(1-square(rho)));
        return fmod(atan2(Y,X), 2*arma::datum::pi);
    }
} /* namespace pair_int */
