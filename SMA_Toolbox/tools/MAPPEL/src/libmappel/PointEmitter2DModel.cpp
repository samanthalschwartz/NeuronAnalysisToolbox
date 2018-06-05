/** @file PointEmitter2DModel.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-26-2014
 * @brief The class definition and template Specializations for PointEmitter2DModel
 */

#include "PointEmitter2DModel.h"
#include "util.h"

const std::vector<std::string> PointEmitter2DModel::estimator_names(
    { "HeuristicMLE", "CGaussHeuristicMLE", "CGaussMLE", "NewtonMLE", "NewtonRaphsonMLE", "QuasiNewtonMLE", "SimulatedAnnealingMLE"});

PointEmitter2DModel::PointEmitter2DModel(int num_params, const IVecT &size, const VecT &psf_sigma)
    : PointEmitterModel(num_params),
      size(size),
      psf_sigma(psf_sigma)
{
    assert(size.n_elem==2);
    assert(psf_sigma.n_elem==2);
    candidate_eta_x=size(0)*candidate_sample_dist_ratio;
    candidate_eta_y=size(1)*candidate_sample_dist_ratio;
}

PointEmitter2DModel::StatsT
PointEmitter2DModel::get_stats() const
{
    StatsT stats;
    stats["dimensions"]=2;
    stats["sizeX"]=size(0);
    stats["sizeY"]=size(1);
    stats["psfSigmaX"]=psf_sigma(0);
    stats["psfSigmaY"]=psf_sigma(1);
    stats["numParams"]=num_params;
    stats["Xmin"]=0;
    stats["Xmax"]=size(0);
    stats["Ymin"]=0;
    stats["Ymax"]=size(1);
    stats["candidate.etaX"]=candidate_eta_x;
    stats["candidate.etaY"]=candidate_eta_y;
    stats["candidate.etaI"]=candidate_eta_I;
    stats["candidate.etabg"]=candidate_eta_bg;
    return stats;
}
