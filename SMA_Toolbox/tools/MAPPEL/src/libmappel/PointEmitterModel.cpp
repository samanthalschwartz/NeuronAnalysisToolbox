/** @file PointEmitterModel.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-13-2014
 * @brief The class definition and template Specializations for PointEmitterModel
 */

#include "PointEmitterModel.h"
#include "util.h"
#include <omp.h>

const std::vector<std::string>
model_names={"Gauss2DMAP", "Gauss2DMLE", "Gauss2DsMAP", "Gauss2DsMLE", "Blink2DsMAP",
             "GaussHSMAP", "GaussHSsMAP", "Blink2DsMAP"};

PointEmitterModel::PointEmitterModel(int num_params)
    : num_params(num_params)
{
//     enable_all_cpus();
//     omp_set_num_threads(1);
}

std::ostream& operator<<(std::ostream &out, PointEmitterModel &model)
{
    auto stats=model.get_stats();
    out<<"["<<model.name()<<":";
    for(auto it=stats.cbegin(); it!=stats.cend(); ++it) out<<" "<<it->first<<"="<<it->second;
    out<<"]";
    return out;
}
