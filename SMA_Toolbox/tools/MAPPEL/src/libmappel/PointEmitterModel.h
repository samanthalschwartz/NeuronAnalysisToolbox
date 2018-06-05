/** @file PointEmitterModel.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-13-2014
 * @brief The class declaration and inline and templated functions for PointEmitterModel.
 *
 * The base class for all point emitter localization models
 */

#ifndef _POINTEMITTERMODEL_H
#define _POINTEMITTERMODEL_H

#include <fstream>
#include <string>
#include <map>

#include <armadillo>

#include "util.h"
#include "numerical.h"
#include "stencil.h"
#include "display.h"
#include "stackcomp.h"
#include "estimator.h"
#include "mcmc.h"


extern const std::vector<std::string> model_names;


/** @brief A Base type for point emitter localization models.
 *
 * We don't assume much here, so that it is possible to have a wide range of 2D and 3D models.
 *
 * 
 * 
 */
class PointEmitterModel {
public:
    int num_candidate_sampling_phases=0; /**< The number of different sampling phases for candidate selection MCMC.  Each phase changes a different subset of variables.*/

    typedef arma::Col<int> IVecT; /**< A type to represent integer data arrays */
//     typedef arma::vec VecT; /**< A type to represent floating-point data arrays */
    typedef arma::mat MatT; /**< A type to represent floating-point data matricies */
    typedef arma::cube CubeT; /**< A type to represent floating-point data cubes */
    typedef arma::mat ParamVecT; /**< A Vector of parameter values */
    
    typedef std::map<std::string,double> StatsT;  /**< A convenient form for reporting dictionaries of named FP data to matlab */

    /* Constant model parameter information */
    const int num_params;

    PointEmitterModel(int num_params);
    virtual ~PointEmitterModel() {}

    virtual std::string name() const =0;
    virtual StatsT get_stats() const =0;

    ParamVecT make_param_vec(int n) const;
    MatT make_param_mat() const;

    friend std::ostream& operator<<(std::ostream &out, PointEmitterModel &model);
protected:
    double prior_epsilon=1E-6; /**< The amount to keep away parameter values from the singularities of the prior distribtions */
    double candidate_sample_dist_ratio=1./30.; /**< Controls the candidate distribution spread for MCMC stuff */
    double candidate_eta_x; /**< The standard deviation for the normally distributed pertebation to theta_x in the random walk MCMC sampling */
    double candidate_eta_y; /**< The standard deviation for the normally distributed pertebation to theta_y in the random walk MCMC sampling */
    double candidate_eta_I; /**< The standard deviation for the normally distributed pertebation to theta_I in the random walk MCMC sampling */
    double candidate_eta_bg; /**< The standard deviation for the normally distributed pertebation to theta_bg in the random walk MCMC sampling */
};

/* Inline member function definitions */
inline
PointEmitterModel::ParamVecT PointEmitterModel::make_param_vec(int n=0) const
{
    return ParamVecT(num_params, n);
}

inline
PointEmitterModel::MatT PointEmitterModel::make_param_mat() const
{
    return MatT(num_params, num_params);
}

/* Template function definitions */
template<class Model>
inline
typename Model::ImageT
model_image(const Model &model, const typename Model::ParamT &theta) 
{
    return model_image(model, model.make_stencil(theta,false));
}

template<class Model, class rng_t>
inline
typename Model::ImageT
simulate_image(const Model &model, const typename Model::ParamT &theta, rng_t &rng) 
{
    return simulate_image(model, model.make_stencil(theta,false), rng);
}

template<class Model>
inline
double
log_likelihood(const Model &model, const typename Model::ImageT &data_im, 
               const typename Model::ParamT &theta)
{
    return log_likelihood(model, data_im, model.make_stencil(theta,false));
}

template<class Model>
inline
double
relative_log_likelihood(const Model &model, const typename Model::ImageT &data_im, 
               const typename Model::ParamT &theta)
{
    return relative_log_likelihood(model, data_im, model.make_stencil(theta,false));
}

template<class Model>
inline
typename Model::ParamT
model_grad(const Model &model, const typename Model::ImageT &data_im, 
               const typename Model::ParamT &theta)
{
    auto grad=model.make_param();
    model_grad(model, data_im, model.make_stencil(theta), grad);
    return grad;
}

template<class Model>
inline
typename Model::MatT
model_hessian(const Model &model, const typename Model::ImageT &data_im, 
               const typename Model::ParamT &theta)
{
    auto grad=model.make_param();
    auto hess=model.make_param_mat();
    model_hessian(model, data_im, model.make_stencil(theta), grad, hess);
    copy_Usym_mat(hess);
    return hess;
}

template<class Model>
inline
typename Model::MatT
model_positive_hessian(const Model &model, const typename Model::ImageT &data_im, 
              const typename Model::ParamT &theta)
{
    auto grad=model.make_param();
    auto hess=model.make_param_mat();
    model_hessian(model, data_im, model.make_stencil(theta), grad, hess);
    hess = -hess;
    copy_Usym_mat(hess);
    modified_cholesky(hess);
    cholesky_convert_full_matrix(hess); //convert from internal format to a full (poitive definite) matrix
    return hess;
}


/** @brief Calculate the Cramer-Rao lower bound at the given paramters
 * @param[in] theta The parameters to evaluate the CRLB at
 * @param[out] crlb The calculated parameters
 */
template<class Model>
typename Model::ParamT
cr_lower_bound(const Model &model, const typename Model::Stencil &s)
{
    return arma::pinv(fisher_information(model,s)).eval().diag();
}

template<class Model>
inline
typename Model::ParamT
cr_lower_bound(const Model &model, const typename Model::ParamT &theta) 
{
    return cr_lower_bound(model,model.make_stencil(theta));
}

template<class Model>
inline
typename Model::MatT
fisher_information(const Model &model, const typename Model::ParamT &theta) 
{
    return fisher_information(model,model.make_stencil(theta));
}


#endif /* _POINTEMITTERMODEL_H */
