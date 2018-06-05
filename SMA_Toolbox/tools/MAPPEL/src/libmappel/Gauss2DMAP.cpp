/** @file Gauss2DMAP.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-25-2014
 * @brief The class definition and template Specializations for Gauss2DMAP
 */
#include <cmath>
#include <algorithm>

#include "Gauss2DMAP.h"
#include "cGaussMLE/cGaussMLE.h"
#include "cGaussMLE/GaussLib.h"


using arma::datum;

const std::vector<std::string> Gauss2DMAP::hyperparameter_names(
    { "Beta_pos", "Mean_I", "Kappa_I", "Mean_bg", "Kappa_bg" });

/** @brief Create a new Model for 2D Gaussian Point Emitters with known PSF under uniform priors.
 * @param[in] size The width and hegiht of the image in pixels
 * @param[in] PSFSigma The standard deviation of the Gaussian PSF
 * Also initializes internal precomputed computational stencils and seed rng.
 */
Gauss2DMAP::Gauss2DMAP(const IVecT &size, const VecT &psf_sigma)
    : Gauss2DModel(size,psf_sigma),
      pos_dist(BetaRNG(beta_pos,beta_pos)),
      I_dist(GammaRNG(kappa_I,mean_I/kappa_I)),
      bg_dist(GammaRNG(kappa_bg,mean_bg/kappa_bg)),
      log_prior_pos_const(log_prior_beta_const(beta_pos)),
      log_prior_I_const(log_prior_gamma_const(kappa_I,mean_I)),
      log_prior_bg_const(log_prior_gamma_const(kappa_bg,mean_bg))
{
    candidate_eta_I=mean_I*candidate_sample_dist_ratio;
    candidate_eta_bg=mean_bg*candidate_sample_dist_ratio;
}

Gauss2DMAP::StatsT Gauss2DMAP::get_stats() const
{
    StatsT stats=Gauss2DModel::get_stats();
    stats["hyperparameter.Beta_pos"]=beta_pos;
    stats["hyperparameter.Mean_I"]=mean_I;
    stats["hyperparameter.kappa_I"]=kappa_I;
    stats["hyperparameter.Mean_bg"]=mean_bg;
    stats["hyperparameter.Kappa_bg"]=kappa_bg;
    return stats;
}


void Gauss2DMAP::set_hyperparameters(const VecT &hyperparameters)
{
    // Params are {beta_pos, mean_I, kappa_I, mean_bg, kappa_bg}
    beta_pos=check_unit_hyperparameter(hyperparameters(0));
    mean_I=check_positive_hyperparameter(hyperparameters(1));
    kappa_I=check_positive_hyperparameter(hyperparameters(2));
    mean_bg=check_positive_hyperparameter(hyperparameters(3));
    kappa_bg=check_positive_hyperparameter(hyperparameters(4));
    log_prior_pos_const=log_prior_beta_const(beta_pos);
    log_prior_I_const=log_prior_gamma_const(kappa_I,mean_I);
    log_prior_bg_const=log_prior_gamma_const(kappa_bg,mean_bg);
    //Reset distributions
    pos_dist.set_params(beta_pos, beta_pos);
    I_dist.kappa(kappa_I);
    I_dist.theta(mean_I/kappa_I);
    bg_dist.kappa(mean_bg);
    bg_dist.theta(mean_bg/kappa_bg);
}

bool Gauss2DMAP::theta_in_bounds(const ParamT &theta) const
{
    bool xOK = (theta(0)>=prior_epsilon) && (theta(0)<=size(0)-prior_epsilon);
    bool yOK = (theta(1)>=prior_epsilon) && (theta(1)<=size(1)-prior_epsilon);
    bool IOK = (theta(2)>=prior_epsilon);
    bool bgOK = (theta(3)>=prior_epsilon);
    return xOK && yOK && IOK && bgOK;
}

void Gauss2DMAP::bound_theta(ParamT &theta) const
{
    theta(0)=restrict_value_range(theta(0), prior_epsilon, size(0)-prior_epsilon); // Prior: Support on [0,size]
    theta(1)=restrict_value_range(theta(1), prior_epsilon, size(1)-prior_epsilon); // Prior: Support on [0,size]
    theta(2)=std::max(prior_epsilon,theta(2));// Prior: Support on [0, inf)
    theta(3)=std::max(prior_epsilon,theta(3));// Prior: Support on [0, inf)
}

double Gauss2DMAP::prior_log_likelihood(const Stencil &s) const
{
    double rllh=prior_relative_log_likelihood(s);
    return rllh+ 2*log_prior_pos_const + log_prior_I_const + log_prior_bg_const;
}

double Gauss2DMAP::prior_relative_log_likelihood(const Stencil &s) const
{
    double xrllh=rllh_beta_prior(beta_pos, s.x(), size(0));
    double yrllh=rllh_beta_prior(beta_pos, s.y(), size(1));
    double Irllh=rllh_gamma_prior(kappa_I, mean_I, s.I());
    double bgrllh=rllh_gamma_prior(kappa_bg, mean_bg, s.bg());
    return xrllh+yrllh+Irllh+bgrllh;
}

Gauss2DMAP::ParamT
Gauss2DMAP::prior_grad(const Stencil &s) const
{
    ParamT grad=make_param();
    grad(0)=beta_prior_grad(beta_pos, s.x(), size(0));
    grad(1)=beta_prior_grad(beta_pos, s.y(), size(1));
    grad(2)=gamma_prior_grad(kappa_I, mean_I, s.I());
    grad(3)=gamma_prior_grad(kappa_bg, mean_bg, s.bg());
    return grad;
}

Gauss2DMAP::ParamT
Gauss2DMAP::prior_grad2(const Stencil &s) const
{
    ParamT grad2=make_param();
    grad2(0)= beta_prior_grad2(beta_pos, s.x(), size(0));
    grad2(1)= beta_prior_grad2(beta_pos, s.y(), size(1));
    grad2(2)= gamma_prior_grad2(kappa_I, s.I());
    grad2(3)= gamma_prior_grad2(kappa_bg, s.bg());
    return grad2;
}

Gauss2DMAP::ParamT
Gauss2DMAP::prior_cr_lower_bound(const Stencil &s) const
{
    //TODO complete these calculations
    ParamT pcrlb=make_param();
    pcrlb.zeros();
    return pcrlb;
}

/* Template Specializations */
template<>
Gauss2DMAP::Stencil
CGaussHeuristicMLE<Gauss2DMAP>::compute_estimate(const ImageT &im, const ParamT &theta_init)
{
    ParamT theta_est(arma::fill::zeros);
    if(model.size(0)==model.size(1) && model.psf_sigma(0)==model.psf_sigma(1)){ //only works for square images and iso-tropic psf
        float Nmax;
        arma::fvec4 ftheta_est;
        //Convert from double
        arma::fmat fim=arma::conv_to<arma::fmat>::from(im);
        //Compute
        CenterofMass2D(model.size(0), fim.memptr(), &ftheta_est[0], &ftheta_est[1]);
        GaussFMaxMin2D(model.size(0), model.psf_sigma(0), fim.memptr(), &Nmax, &ftheta_est[3]);
        ftheta_est[2]=std::max(0., (Nmax-ftheta_est[3])*2*arma::datum::pi*model.psf_sigma(0)*model.psf_sigma(0));
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
CGaussMLE<Gauss2DMAP>::compute_estimate(const ImageT &im, const ParamT &theta_init, ParamT &theta, ParamT &crlb, double &llh)
{
    if(model.size(0)==model.size(1) && model.psf_sigma(0)==model.psf_sigma(1)){//only works for square images and iso-tropic psf
        float fllh;
        arma::fvec4 fcrlb, ftheta;
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
