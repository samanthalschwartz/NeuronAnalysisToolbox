/** @file Gauss2DMLE.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-25-2014
 * @brief The class definition and template Specializations for Gauss2DMLE
 */
#include <algorithm>
#include <memory>

#include "Gauss2DMLE.h"
#include "cGaussMLE/cGaussMLE.h"
#include "cGaussMLE/GaussLib.h"

/* Constant model estimator names: These are the estimator names we have defined for this class */

const std::vector<std::string> Gauss2DMLE::hyperparameter_names(
    { "I_min", "I_max", "bg_min", "bg_max" });

/** @brief Create a new Model for 2D Gaussian Point Emitters with known PSF under uniform priors.
 * @param[in] size The width and hegiht of the image in pixels
 * @param[in] PSFSigma The standard deviation of the Gaussian PSF
 * Also initializes internal precomputed computational stencils and seed rng.
 */
Gauss2DMLE::Gauss2DMLE(const IVecT &size, const VecT &psf_sigma)
    : Gauss2DModel(size,psf_sigma),
      pos_dist(UniformRNG(0,1)),
      I_dist(UniformRNG(I_min,I_max)),
      bg_dist(UniformRNG(bg_min,bg_max))
{
    candidate_eta_I=(0.5*(I_max-I_min))*candidate_sample_dist_ratio;
    candidate_eta_bg=(0.5*(bg_max-bg_min))*candidate_sample_dist_ratio;
}

Gauss2DMLE::StatsT Gauss2DMLE::get_stats() const
{
    StatsT stats=Gauss2DModel::get_stats();
    stats["hyperparameter.I_min"]=I_min;
    stats["hyperparameter.I_max"]=I_max;
    stats["hyperparameter.BG_min"]=bg_min;
    stats["hyperparameter.BG_max"]=bg_max;
    return stats;
}

void Gauss2DMLE::set_hyperparameters(const VecT &hyperparameters)
{
    I_min=hyperparameters(0);
    I_max=hyperparameters(1);
    bg_min=hyperparameters(2);
    bg_max=hyperparameters(3);
    //Reset distributions
    I_dist.a(I_min);
    I_dist.b(I_max);
    bg_dist.a(bg_min);
    bg_dist.b(bg_max);
}

bool Gauss2DMLE::theta_in_bounds(const ParamT &theta) const
{
    bool xOK =  (theta(0)>=0)      && (theta(0)<=size(0));
    bool yOK =  (theta(1)>=0)      && (theta(1)<=size(1));
    bool IOK =  (theta(2)>=I_min)  && (theta(2)<=I_max);
    bool bgOK = (theta(3)>=bg_min) && (theta(3)<=bg_max);
    return xOK && yOK && IOK && bgOK;
}

void Gauss2DMLE::bound_theta(ParamT &theta) const
{
    theta(0)=restrict_value_range(theta(0), 0, size(0));       // Prior: Uniform on [0,sizeX]
    theta(1)=restrict_value_range(theta(1), 0, size(1));       // Prior: Uniform on [0,sizeY]
    theta(2)=restrict_value_range(theta(2), I_min, I_max);  // Prior: Uniform on [I_min,I_max]
    theta(3)=restrict_value_range(theta(3), bg_min, bg_max);// Prior: Uniform on [bg_min,bg_max]
}


double Gauss2DMLE::prior_log_likelihood(const Stencil &s) const
{
    return 0;
}

double Gauss2DMLE::prior_relative_log_likelihood(const Stencil &s) const
{
    return 0;
}

Gauss2DMLE::ParamT
Gauss2DMLE::prior_grad(const Stencil &s) const
{
    return ParamT(arma::fill::zeros);
}

Gauss2DMLE::ParamT
Gauss2DMLE::prior_grad2(const Stencil &s) const
{
    return ParamT(arma::fill::zeros);
}

Gauss2DMLE::ParamT
Gauss2DMLE::prior_cr_lower_bound(const Stencil &s) const
{
    return ParamT(arma::fill::zeros);
}



/* Template Specializations */
template<>
Gauss2DMLE::Stencil
CGaussHeuristicMLE<Gauss2DMLE>::compute_estimate(const ImageT &im, const ParamT &theta_init)
{
    Gauss2DMLE::ParamT theta_est(arma::fill::zeros);
    if(model.size(0)==model.size(1) && model.psf_sigma(0)==model.psf_sigma(1)){ //only works for square images and iso-tropic psf
        float Nmax;
        arma::fvec4 ftheta_est;
        //Convert from double
        arma::fmat fim=arma::conv_to<arma::fmat>::from(im);
        //Compute
        CenterofMass2D(model.size(0), fim.memptr(), &ftheta_est[0], &ftheta_est[1]);
        GaussFMaxMin2D(model.size(0), model.psf_sigma(0), fim.memptr(), &Nmax, &ftheta_est[3]);
        ftheta_est[2]=std::max(0., (Nmax-ftheta_est[3]) * 2 * arma::datum::pi * model.psf_sigma(0) * model.psf_sigma(0));
        //Back to double
        theta_est=arma::conv_to<arma::mat>::from(ftheta_est);
        //Swap x/y and add .5 tp convert from CGauss to mappel coordinates
        float temp=theta_est(0)+.5;
        theta_est(0)=theta_est(1)+.5;
        theta_est(1)=temp;
    }
    return model.make_stencil(theta_est);
}

template<>
void
CGaussMLE<Gauss2DMLE>::compute_estimate(const ImageT &im, const ParamT &theta_init, ParamT &theta, ParamT &crlb, double &llh)
{
    if(model.size(0)==model.size(1) && model.psf_sigma(0)==model.psf_sigma(1)){//only works for square images and iso-tropic psf
        float fllh;
        arma::fvec4 ftheta, fcrlb;
        //Convert from double
        arma::fmat fim=arma::conv_to<arma::fmat>::from(im);
        //Compute
        MLEFit(fim.memptr(), model.psf_sigma(0), model.size(0), max_iterations, ftheta.memptr(), fcrlb.memptr(), &fllh);
        //Back to double
        theta=arma::conv_to<arma::vec>::from(ftheta);
        crlb=arma::conv_to<arma::vec>::from(fcrlb);
        //Swap x/y and add .5 tp convert from CGauss to mappel coordinates
        float temp=theta(0)+.5;
        theta(0)=theta(1)+.5;
        theta(1)=temp;
        llh=log_likelihood(model,im,model.make_stencil(theta));
    } else {
        theta.zeros();
        crlb.zeros();
        llh=0.0;
    }
}


