
/** @file Gauss2DMLE.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-22-2014
 * @brief The class declaration and inline and templated functions for Gauss2DMLE.
 */

#ifndef _GAUSS2DMLE_H
#define _GAUSS2DMLE_H

#include "Gauss2DModel.h"


/** @brief A 2D Likelyhood model for a point emitter localization with 
 * Symmetric Gaussian PSF and Poisson Noise, using a  Uniform Prior over the 
 * parameter vectors
 *
 * This model matches the model used in cGaussMLE.  
 * So we can use this as comparison.
 * 
 */
class Gauss2DMLE : public Gauss2DModel {
private:
    /* Theta prior parameters */
    double I_min=1e2; /**< The minimum intensity for our Uniform prior */
    double I_max=1e5; /**< The maximum intensity for our Uniform prior */
    double bg_min=1.0e-6; /**< The minimum bg for our Uniform prior (estimating bg=0 is bad for the numerics) */
    double bg_max=1e2; /**< The maximum bg for our Uniform prior */
    UniformRNG pos_dist;
    UniformRNG I_dist;
    UniformRNG bg_dist;
public:
    static const std::vector<std::string> hyperparameter_names;

    /* Constructor/Destructor */
    Gauss2DMLE(const IVecT &size, const VecT &psf_sigma);

    /* Model values setting and information */
    std::string name() const {return "Gauss2DMLE";}
    StatsT get_stats() const;

    /* Sample from Theta Prior */
    ParamT sample_prior(RNG &rng);
    void set_hyperparameters(const VecT &hyperparameters);
    void bound_theta(ParamT &theta) const;
    bool theta_in_bounds(const ParamT &theta) const;

    double prior_log_likelihood(const Stencil &s) const;
    double prior_relative_log_likelihood(const Stencil &s) const;
    ParamT prior_grad(const Stencil &s) const;
    ParamT prior_grad2(const Stencil &s) const;
    ParamT prior_cr_lower_bound(const Stencil &s) const;
};

/* Template Specialization Declarations */

/* Inlined Methods */
inline
Gauss2DMLE::ParamT
Gauss2DMLE::sample_prior(RNG &rng)
{
    ParamT theta;
    theta(0)=size(0)*pos_dist(rng); // Prior: Uniform on [0,size]
    theta(1)=size(1)*pos_dist(rng); // Prior: Uniform on [0,size]
    theta(2)=I_dist(rng);   // Prior: Uniform on [I_min,I_max]
    theta(3)=bg_dist(rng);  // Prior: Uniform on [bg_min,bg_max]
    return theta;
}



#endif /* _GAUSS2DMLE_H */
