/** @file pairinteractionsmc.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief Sequential Monte Carlo sampling of pair interaction model
 *
 */
#include <cmath>
#include <random>

#include "pairinteractionsmc.h"
#include "unscentedtransform.h"
#include "vonmises.h"
#include "gaussianangle.h"

using namespace pair_int;

PairInteractionSMC::PairInteractionSMC(const ParamsT &params)
{
    initDefaultPrior();
    setParams(params);
}

PairInteractionSMC::PairInteractionSMC(const ParamsT &params, const MatListT &data)
{
    initDefaultPrior();
    setParams(params);
    addData(data);
}

PairInteractionSMC::IdxT PairInteractionSMC::addData(const MatListT &data)
{
    IdxT N = data.size();
    for(IdxT n=0;n<N; n++) if(data[n].n_rows!=9) std::length_error("Data should have 9 rows");
    clearParticles();
    pair_data.insert(pair_data.end(), data.begin(),data.end());
    IdxT Ndata = pair_data.size();
    pair_Nobs.set_size(Ndata);
    for(IdxT n=0;n<Ndata;n++) pair_Nobs(n) = pair_data[n].n_cols;
    return Ndata;
}

PairInteractionSMC::IdxT PairInteractionSMC::addData(const CubeT &data)
{
    IdxT N = data.n_slices;
    for(IdxT n=0;n<N; n++) if(data.slice(n).n_rows!=9) std::length_error("Data should have 9 rows");
    clearParticles();
    for(IdxT n=0;n<N; n++) pair_data.emplace_back(data.slice(n));
    IdxT Ndata = pair_data.size();
    pair_Nobs.set_size(Ndata);
    for(IdxT n=0;n<Ndata;n++) pair_Nobs(n) =pair_data[n].n_cols;
    return Ndata;
}

PairInteractionSMC::IdxT PairInteractionSMC::addData(const MatT &data)
{
    if(data.n_rows!=9) std::length_error("Data should have 9 rows");
    clearParticles();
    pair_data.emplace_back(data);
    IdxT Ndata = pair_data.size();
    pair_Nobs.set_size(Ndata);
    for(IdxT n=0;n<Ndata;n++) pair_Nobs(n) = pair_data[n].n_cols;
    return Ndata;
}

void PairInteractionSMC::clearData()
{
    clearParticles();
    pair_data.clear();
    pair_Nobs.reset();
}

PairInteractionSMC::MatListT PairInteractionSMC::getData()
{
    return pair_data;
}

void PairInteractionSMC::setParams(const ParamsT &params)
{
    for(auto &param: params) {
        auto name = param.first;
        if     (name == "DA" && param.second>0)          DA = param.second;
        else if(name == "DB" && param.second>0)          DB = param.second;
        else if(name == "DAB" && param.second>0)         DAB = param.second;
        else if(name == "Dphi" && param.second>0)        Dphi = param.second;
        else if(name == "rho" && param.second>0)         rho = param.second;
        else if(name == "sigma_rho" && param.second>0)   sigma_rho = param.second;
        else if(name == "gamma" && param.second>0)       gamma = param.second;
        else if(name == "rho_bind" && param.second>0)    rho_bind = param.second;
        else if(name == "lambda_bind" && param.second>0) lambda_bind = param.second;
        else if(name == "koff" && param.second>0)        koff = param.second;
    }
    valid_params = (DA>0 && DB>0 && DAB>0 && Dphi>0 && rho>0 && sigma_rho>0 && gamma>0 && rho_bind>0 && lambda_bind>0 && koff>0);
    if (valid_params){
        rd.setParams(DA+DB,rho_bind,lambda_bind);
    }
    sigma_rho_sq = square(sigma_rho);
}

double PairInteractionSMC::updateParamsDefault(ParamsT &params, const std::string &key, double default_val) const
{
    auto it = params.find(key);
    if(it==params.end()) {
        params.emplace(key,default_val);
        return default_val;
    } else {
        return it->second;
    }
}

PairInteractionSMC::ParamsT PairInteractionSMC::getParams() const
{
    ParamsT params;
    params.emplace("DA",DA);
    params.emplace("DB",DB);
    params.emplace("DAB",DAB);
    params.emplace("Dphi",Dphi);
    params.emplace("rho",rho);
    params.emplace("sigma_rho",sigma_rho);
    params.emplace("gamma",gamma);
    params.emplace("rho_bind",rho_bind);
    params.emplace("lambda_bind",lambda_bind);
    params.emplace("koff",koff);
    params.emplace("valid_params",valid_params);
    return params;
}

PairInteractionSMC::ParamsT PairInteractionSMC::getPriorParams() const
{
    ParamsT params =  pos_prior->getParams();
    params["priorState0"]=state_prior(0);
    params["priorState1"]=state_prior(1);
    return params;
}

void PairInteractionSMC::initDefaultPrior()
{
    VecT rect{0,0,1,1};
    pos_prior = std::make_shared<RectangularPositionPrior>(rect);
    log_state_prior = arma::log(state_prior);
}

void PairInteractionSMC::setPriorState(const VecT &_state_prior) 
{
    if(_state_prior.n_elem != 2) throw std::invalid_argument("state prior must have size 2.");
    double s = arma::sum(_state_prior);
    double tol = 1e-8;
    if(fabs(s-1)>tol) throw std::invalid_argument("state prior must have sum==1");
    state_prior = _state_prior;
    log_state_prior = arma::log(state_prior);
}

void PairInteractionSMC::setPriorPositionGaussian(const VecT &mean, const MatT &cov)
{
    pos_prior = std::make_shared<GaussianPositionPrior>(mean, cov);
}

void PairInteractionSMC::setPriorPositionRectangular(const VecT &rect)
{
    pos_prior = std::make_shared<RectangularPositionPrior>(rect);
}

void PairInteractionSMC::setPriorPositionCircular(const VecT &center, double radius)
{
    pos_prior = std::make_shared<CircularPositionPrior>(center, radius);
}

void PairInteractionSMC::clearParticles()
{
    valid_particles=false;
    Nparticles=0;
    pair_state.clear();
    pair_ancestor.clear();
    pair_weights.reset();
    pair_llh.reset();
    pair_obs_llh.reset();
}

void PairInteractionSMC::initializeParticles(IdxT N)
{
    if(!valid_params) throw std::logic_error("Invalid parameters");
    if(pair_data.empty()) throw std::logic_error("No valid Data");
    if(N<1) throw std::out_of_range("Must have at least 1 particle.");
    Nparticles=N;
    IdxT Ndata = pair_data.size();
    
    //Pre-allocate these so we can parallelize
    pair_state.resize(Ndata);
    pair_ancestor.resize(Ndata);
    
    pair_weights.set_size(Nparticles, Ndata);
    pair_weights.zeros();

    pair_llh.set_size(Nparticles, Ndata);
    pair_llh.zeros();
    
    pair_obs_llh.set_size(Ndata);
    pair_obs_llh.zeros();
    
    valid_particles=false;
}

PairInteractionSMC::VecT PairInteractionSMC::obsLLH() const
{
    if(!valid_particles) throw std::logic_error("Particle SMC status:invalid");
    return pair_obs_llh;
}

void PairInteractionSMC::sampleParticle(MatListT &L, VecT &llh)
{
    if(!valid_particles) throw std::logic_error("Particle SMC status:invalid");
    IdxT Ndata = pair_data.size();
    L.clear();
    llh.set_size(Ndata);
    for(IdxT d=0; d<Ndata; d++) {
        IdxT pidx = drawParticleSample(d);
        L.emplace_back(reconstructParticle(d,pidx));
        llh(d) = pair_llh(pidx,d);
    }
}

PairInteractionSMC::VecT PairInteractionSMC::sampleParticleLLH()
{
    if(!valid_particles) throw std::logic_error("Particle SMC status:invalid");
    IdxT Ndata = pair_data.size();
    VecT llh(Ndata);
    for(IdxT d=0; d<Ndata; d++) llh(d) = pair_llh(drawParticleSample(d),d);
    return llh;
}

void PairInteractionSMC::getAllParticles(int obsIdx, CubeListT &L, MatT &weights, MatT &llh) const
{
    if(!valid_particles) throw std::logic_error("Particle SMC status:invalid");
    L.clear();
    IdxT Ndata = pair_data.size();
    weights.set_size(Nparticles,Ndata);
    llh.set_size(Nparticles,Ndata);
    for(IdxT d=0; d<Ndata; d++) {
        int T = obsIdx;
        int Nobs = static_cast<int>(pair_Nobs(d));
        if(T<0 || T>=Nobs) T=Nobs-1; //Default to all observations
        CubeT particles(5,T+1,Nparticles);
        for(IdxT n=0; n<Nparticles; n++) particles.slice(n) = reconstructParticle(d,n,T);
        L.emplace_back(particles);
    }
    weights = pair_weights;
    llh = pair_llh;
}

bool PairInteractionSMC::runParticleFilter(int Nparticles, const std::string &propTypeStr)
{
    initializeParticles(Nparticles);
    IdxT Ndata = pair_data.size();
    auto propType = getProposalType(propTypeStr);
    #pragma omp parallel 
    {
        VecT ws,llh;
        IMatT ancestor;
        CubeT state;
        MatT data;
        double oLLH;
        #pragma omp for
        for(IdxT n=0;n<Ndata;n++) {
            data = pair_data[n];
            runParticleFilter(propType, data, Nparticles, state, ancestor,ws,llh,oLLH);
            pair_weights.col(n) = ws;
            pair_llh.col(n) = llh;
            pair_obs_llh(n) = oLLH;
            pair_state[n] = state;
            pair_ancestor[n] = ancestor;
        }
    }
    if(Ndata>0) {
        valid_particles=true;
    }
    return valid_particles;
}

/**
 * 
 * 
 * @param[in] obsIdx Should be an index between 0 and Nobs-1 indicating the last observation to sample at.  [-1 to get a particle for all observations]
 * 
 */
PairInteractionSMC::MatT PairInteractionSMC::reconstructParticle(IdxT data_idx, IdxT pidx, int obsIdx) const
{
    auto &state = pair_state[data_idx];
    auto &anc = pair_ancestor[data_idx];
    int Nobs = static_cast<int>(pair_Nobs(data_idx));
    if(obsIdx<0 || obsIdx>=Nobs) obsIdx=Nobs-1; //Default to all observations
    MatT particle(5,obsIdx+1);
    for(int t=obsIdx; t>=0; t--) {
        particle.col(t) = state.slice(pidx).col(t);
        if(t>0) pidx = anc(pidx,t-1);
    }
    return particle;
}

PairInteractionSMC::ProposalType PairInteractionSMC::getProposalType(const std::string &propTypeStr) const
{
    ProposalType typ;
    if (propTypeStr=="Transition") {
        typ = ProposalType::TRANSITION;
    } else if(propTypeStr=="Observation"){
        typ = ProposalType::OBSERVATION;
    } else if(propTypeStr=="Combined"){
        typ = ProposalType::COMBINED;
    } else {
        throw std::runtime_error("Unknown proposal type");
    }
    return typ;
}

bool PairInteractionSMC::runParticleFilter(ProposalType propType, const MatT &O, IdxT N, CubeT &state, IMatT &ancestor, 
                                           VecT &weights, VecT& llh, double &obsLLH) const
{
    VecT dt = diff(O.row(0).t());
    MatT xs = O.rows(1,4);
    MatT ws = O.rows(5,8);
    IdxT T = O.n_cols; //number of time indexes (observations)
    state.set_size(5,T,N); //Ordering most efficient for our traversal strategy
    ancestor.set_size(N,T-1);
    weights.set_size(N);
    llh.set_size(N);
    
    VecT llhP(N), llhQ(N);
    VecT llh_G(N), llh_H(N);
    //Temporaries
    VecT zt(N),z0(N); 
    MatT yt(4,N),y0(4,N);
    //Sample the prior proposal distribution for state given true positions
    //  (z0,y0) ~ q(z0, y0 | x0 w0) = q(z0 | y0) q(y0 | x0 w0).
    sampleQ0(xs.col(0), ws.col(0), yt, zt,  llhQ);
    state.tube(0,0) = zt;
    state.tube(1,0,4,0) = yt;
    for(IdxT n=0;n<N;n++) llhP(n) = llhG0(yt.col(n), zt(n)) + llhH(xs.col(0), ws.col(0), yt.col(n));

    llh = llhP; //Set log likelihood
    VecT logw = llhP-llhQ;
    weights = logNormalize(logw);

    for(IdxT t=1; t<T; t++) {
        double ess = 1./arma::sum(weights%weights);
//         std::cout<<"[T: "<<t<<"] --";
        if(ess < std::max(MIN_ESS, N*MIN_ESS_RATIO)) {
//             std::cout<<"<Resampling>\n * Ancestors preserved:";
            ancestor.col(t-1) = RNG.resample_dist(weights,N);
//             std::cout<<arma::unique(ancestor.col(t-1)).eval().n_elem<<"/"<<N<<"\n";
            VecT tmp_logw = logw(ancestor.col(t-1));
            VecT tmp_llh = llh(ancestor.col(t-1));
            logw = tmp_logw;
            llh = tmp_llh;
            //Establish y0 and z0
            for(IdxT n=0;n<N; n++) {
                IdxT anc = ancestor(n,t-1);
                z0(n) = state(0,t-1,anc);
                for(IdxT k=0;k<4; k++) y0(k,n) = state(k+1,t-1,anc);
            }
        } else { //No resample
//             std::cout<<"\n";
            for(IdxT n=0;n<N;n++) {
                ancestor(n,t-1) = n;
                z0(n) = state(0,t-1,n);
                for(IdxT k=0;k<4; k++) y0(k,n) = state(k+1,t-1,n);
            }
//             z0 = state.tube(0,t-1);
//             y0 = state.tube(1,t-1,4,t-1);
        }
//         std::cout<<" * Initially Bound: "<<arma::sum(z0)<<"/"<<N<<"\n";
        sampleProposal(propType, xs.col(t), ws.col(t), y0, z0, dt(t-1), yt, zt, llhQ);

        if(propType == ProposalType::TRANSITION) { 
            // if G()==Q() is proposal density we already have llhG == llhQ
            llhP = llhQ + llhH_vec(xs.col(t), ws.col(t), yt); 
        } else {
            llh_G = llhG_vec(yt, zt, y0, z0, dt(t-1));
            llh_H = llhH_vec(xs.col(t), ws.col(t), yt);
            llhP = llh_G+llh_H;
        } 
        state.tube(0,t) = zt;
        state.tube(1,t,4,t) = yt;
        llh += llhP; //Add to log likelihood
        logw += llhP-llhQ;
        weights = logNormalize(logw);
//         std::cout<<"llhG: "<<llh_G.t();
//         std::cout<<"llhH: "<<llh_H.t();
//         std::cout<<"llhP: "<<llhP.t();
//         std::cout<<"llhP~: "<<llhP.t()-arma::max(llhP);
//         std::cout<<"llh: "<<llh.t();
//         std::cout<<"llh~: "<<llh.t()-arma::max(llh);
//         std::cout<<"llhQ: "<<llhQ.t();
//         std::cout<<"logA: "<<(llhP-llhQ).eval().t();
//         std::cout<<"Weights: "<<weights.t();
//         std::cout<<"  zt: "<<zt.t();
//         std::cout<<"-----------------------------------------------------------\n";
    }
    obsLLH = logSum(logw)-log(N);
    return true;
}


PairInteractionSMC::VecT PairInteractionSMC::computeLLH(const MatListT &obs_data_arr, const MatListT &particle_arr) const
{
    if(particle_arr.size()!=obs_data_arr.size()) throw std::length_error("Particles array must have same length as observation array");
    if(!valid_params) throw std::logic_error("Invalid parameters");
    IdxT N = obs_data_arr.size();
    for(IdxT n=0;n<N;n++) {
        if(particle_arr[n].n_rows!=5) throw std::length_error("Particles must have 5 rows");
        if(obs_data_arr[n].n_rows!=9) throw std::length_error("Observations must have 9 rows");
        if(obs_data_arr[n].n_cols!=particle_arr[n].n_cols) throw std::length_error("Observations and particles must have same number of columns");
    }
    VecT llh(N);
    #pragma omp parallel for
    for(IdxT n=0;n<N;n++) llh(n) = computeLLH(obs_data_arr[n],particle_arr[n]);
    return llh;
}

double PairInteractionSMC::computeLLH(IdxT data_idx, const MatT &particle) const
{
    if(!valid_params) throw std::logic_error("Invalid parameters");
    if(particle.n_rows!=5) throw std::length_error("Particles must have 5 rows");
    return computeLLH(pair_data[data_idx],particle);
}

PairInteractionSMC::VecT PairInteractionSMC::computeLLH(const MatListT &particle_arr) const
{
    if(particle_arr.size()!=pair_data.size()) throw std::length_error("Particles array must have same length of NData");
    if(!valid_params) throw std::logic_error("Invalid parameters");
    IdxT N=pair_data.size();
    for(auto &p: particle_arr) if(p.n_rows!=5) throw std::length_error("Particles must have 5 rows");
    VecT llh(N);
    #pragma omp parallel for
    for(IdxT n=0;n<N;n++) llh(n) = computeLLH(pair_data[n],particle_arr[n]);
    return llh;
}

PairInteractionSMC::VecT PairInteractionSMC::computeLLH(const CubeT &particle_arr) const
{
    if(particle_arr.n_slices!=pair_data.size()) throw std::length_error("Particles array cube must have same #slices as NData");
    if(!valid_params) throw std::logic_error("Invalid parameters");
    IdxT N=particle_arr.n_slices;
    if(particle_arr.n_rows!=5) throw std::length_error("Particles must have 5 rows");
    VecT llh(N);
    #pragma omp parallel for
    for(IdxT n=0;n<N;n++) llh(n) = computeLLH(pair_data[n],particle_arr.slice(n));
    return llh;
}

double PairInteractionSMC::computeLLH(const MatT &obs_data, const MatT &particle) const
{
    IdxT Nobs = obs_data.n_cols;
    if(!valid_params) throw std::logic_error("Invalid parameters");
    if(obs_data.n_cols!=particle.n_cols) throw std::length_error("Observation and particle have different lengths");
    if(obs_data.n_rows!=9) throw std::length_error("Observations must have 9 rows");
    if(particle.n_rows!=5) throw std::length_error("Particles must have 5 rows");
    VecT dt=arma::diff(obs_data.row(0).t());
    double llh = 0.;
    double z0 = particle(0,0);
    MatT y0 = particle.col(0).rows(1,4);
    
    auto x0 = obs_data.col(0).rows(1,4);
    auto w0 = obs_data.col(0).rows(5,8);
    llh += llhG0(y0,z0);
    llh += llhH(x0,w0,y0);
    
    for(IdxT n=1;n<Nobs;n++){
        z0 = particle(0,n-1);
        y0 = particle.col(n-1).rows(1,4);
        auto zn = particle(0,n);
        auto yn = particle.col(n).rows(1,4);
        auto xn = obs_data.col(n).rows(1,4);
        auto wn = obs_data.col(n).rows(5,8);
        llh += llhG(yn,zn,y0,z0,dt(n-1));
        llh += llhH(xn,wn,yn);
    }
    return llh;
}


PairInteractionSMC::MatT 
PairInteractionSMC::computeLLH_debug(const MatT &obs_data, const MatT &particle) const
{
    using arma::span;
    IdxT Nobs = obs_data.n_cols;
    if(!valid_params) throw std::logic_error("Invalid parameters");
    if(obs_data.n_cols!=particle.n_cols) throw std::length_error("Observation and particle have different lengths");
    if(obs_data.n_rows!=9) throw std::length_error("Observations must have 9 rows");
    if(particle.n_rows!=5) throw std::length_error("Particles must have 5 rows");
    VecT dt=arma::diff(obs_data.row(0).t());
    MatT llh(9,Nobs);
    auto z0 = particle(0,0);
    auto y0 = particle.col(0).rows(1,4);
    auto x0 = obs_data.col(0).rows(1,4);
    auto w0 = obs_data.col(0).rows(5,8);
    llh.col(0).rows(0,4) = llhG0_all(y0,z0);
    llh.col(0).rows(5,8) = llhH_all(x0,w0,y0);
    
    for(IdxT n=1;n<Nobs;n++){
        z0 = particle(0,n-1);
        y0 = particle.col(n-1).rows(1,4);
        auto zn = particle(0,n);
        auto yn = particle.col(n).rows(1,4);
        auto xn = obs_data.col(n).rows(1,4);
        auto wn = obs_data.col(n).rows(5,8);
        llh.col(n).rows(0,4) = llhG_all(yn,zn,y0,z0,dt(n-1));
        llh.col(n).rows(5,8) = llhH_all(xn,wn,yn);
    }
    return llh;
}

void PairInteractionSMC::simulate(IdxT Nparticles, IdxT Nsteps, double tmax, const std::string &simType, 
                                  ParamsT &params, CubeT &obs, CubeT &particles, CubeT &llh) const
{
    params["Nparticles"] = Nparticles;
    params["Nsteps"] = Nsteps;
    params["tmax"] = tmax;
    if(Nparticles<1) throw std::invalid_argument("Nparticles must be >=1");
    if(Nsteps<1) throw std::invalid_argument("Nsteps must be >=1");
    if(tmax<=0) throw std::invalid_argument("tmax must be >0");
    auto modelParams = getParams();
    params.insert(modelParams.begin(), modelParams.end());
    auto priorParams = getPriorParams();
    params.insert(priorParams.begin(), priorParams.end());
    if(Nparticles<1 || Nsteps<1 || tmax<=0) throw std::invalid_argument("Invalid simulation constants.");
    if(simType =="Ancestral") {
        params["simType"] = double(SimulationType::ANCESTRAL);
        simulateAncestral(Nparticles, Nsteps, tmax, params, obs, particles, llh);
    } else if(simType=="Blur") { 
        params["simType"] = double(SimulationType::BLUR);
        simulateBlur(Nparticles, Nsteps, tmax, params, obs, particles, llh);
    } else if (simType=="EC") {
        params["simType"] = double(SimulationType::ERBANCHAPMAN);
        simulateEC(Nparticles, Nsteps, tmax, params, obs, particles, llh);
    } else {
        throw std::invalid_argument("Unknown simType: "+simType);
    }
}

void PairInteractionSMC::simulateAncestral(IdxT Nparticles, IdxT Nsteps, double tmax, 
                                           ParamsT &params, CubeT &obs, CubeT &particles, CubeT &llh) const
{
    obs.set_size(9,Nsteps,Nparticles);//rows are: [t, obs_xa, obs_xb, obs_ya, obs_xb, w_xa, w_xb, w_ya, w_yb]
    particles.set_size(5,Nsteps,Nparticles); //rows are: [z, true_ax true_ay true_bx true_by]
    llh.set_size(9,Nsteps,Nparticles); //rows are: [z, yax yay ybx yby xax xay xbx xby]
    
    double dt = (Nsteps==1) ? 0 : tmax/(Nsteps-1);
    params["dt"] = dt;
    double wmean = updateParamsDefault(params,"wmean",rho);
    double wsigma = updateParamsDefault(params,"wsimga",0.15*wmean);
    double wgamma = updateParamsDefault(params,"wgamma",1./(10*dt));

    
    for(IdxT n=0; n<Nparticles; n++){
        //Sample W
        auto &O=obs.slice(n);
        auto &P=particles.slice(n);
        auto &L=llh.slice(n);
        O.row(5) = sampleOUEquilibriumRegular(Nsteps, dt, wmean, wsigma, wgamma).t();
        O.row(6) = O.row(5);
        O.row(7) = sampleOUEquilibriumRegular(Nsteps, dt, wmean, wsigma, wgamma).t();
        O.row(8) = O.row(7);

        auto Ot=O.col(0);
        auto Pt=P.col(0);
        auto Lt=L.col(0);
        
        double zt;
        VecT yt(4);
        VecT llh_Gt(5);
        sampleG0(yt, zt, llh_Gt);
        Pt(0) = zt;
        Pt.rows(1,4) = yt;
        Lt.rows(0,4) = llh_Gt;

        VecT xt(4);
        VecT llh_Ht(4);
        auto wt = O.col(0).rows(5,8);
        sampleH_all(wt, yt, xt, llh_Ht);
        Ot(0)=0;
        Ot.rows(1,4) = xt;
        Lt.rows(5,8) = llh_Ht;
        for(IdxT t=1;t<Nsteps;t++) {
            auto Ot=O.col(t);
            auto Pt=P.col(t);
            auto Lt=L.col(t);
            double z0 = P(0,t-1);
            auto y0 = P.col(t-1).rows(1,4);
            //Sample forward propogator G
            sampleG(y0, z0, dt, yt, zt, llh_Gt);
            Pt(0) = zt;
            Pt.rows(1,4) = yt;
            Lt.rows(0,4) = llh_Gt;
            //Sample observations H
            auto wt = Ot.rows(5,8);
            sampleH_all(wt, yt, xt, llh_Ht);
            Ot(0)=dt*t;
            Ot.rows(1,4) = xt;
            Lt.rows(5,8) = llh_Ht;
        }
    }
}


void PairInteractionSMC::simulateBlur(IdxT Nparticles, IdxT Nsteps, double tmax, 
                                      ParamsT &params,CubeT &obs, CubeT &particles, CubeT &llh) const
{
    
}

void PairInteractionSMC::simulateEC(IdxT Nparticles, IdxT Nsteps, double tmax, 
                                    ParamsT &params,CubeT &obs, CubeT &particles, CubeT &llh) const
{
    
}


void  PairInteractionSMC::sampleQ0(const VecT &x0, const VecT &w0, MatT &y0, VecT &z0, VecT &llh) const
{
    IdxT N = z0.n_elem;
    if(y0.n_cols!=N || y0.n_rows!=4) throw std::invalid_argument("y0 must be size:[4xN]");
    if(llh.n_elem!=N) throw std::invalid_argument("llh must be size:[N]");
    llh.fill(0);
    VecT llhZ(2);
    for(IdxT n=0; n<N; n++){
        y0.col(n) = sampleGaussian(x0,w0);
        llh(n) += llhH(x0,w0,y0.col(n));// p(y0|x0,w0) which happens to be the exact same as p(x0|y0,w0).
        double p_y0_free = exp(llhG0(y0.col(n),0));  // p(y_0 | z_0=F)*p(z_0=F)
        double p_y0_bound = exp(llhG0(y0.col(n),1)); // p(y_0 | z_0=B)*p(z_0=B)
        //Normalize probabilities
        p_y0_free /= (p_y0_free + p_y0_bound);
        if(RNG.randu() <= p_y0_free) { //free
            z0(n) = 0;
            llh(n) += log(p_y0_free);
        } else { //bound
            z0(n) = 1;
            llh(n) += log(1-p_y0_free);
        }
    }
}

void PairInteractionSMC::sampleG0(VecT &y0, double& z0, VecT &llhAll) const
{
    if (RNG.randu() < state_prior(0)) { // Free
        z0=0;
        llhAll(0)=log_state_prior(0);
        y0.subvec(0,1) = pos_prior->sample_pos();
        llhAll(1) = .5*pos_prior->llh_pos(y0.subvec(0,1));
        llhAll(2) = llhAll(1);
        y0.subvec(2,3) = pos_prior->sample_pos();
        llhAll(3) = .5*pos_prior->llh_pos(y0.subvec(2,3));
        llhAll(4) = llhAll(3);
    } else {
        z0=1;
        llhAll(0)=log_state_prior(1);
        //sample in v-space
        y0.subvec(0,1) = pos_prior->sample_pos();
        llhAll(1) = .5*pos_prior->llh_pos(y0.subvec(0,1));
        llhAll(2) = llhAll(1);
        y0(2) = rho+RNG.randn()*sigma_rho;
        llhAll(3) = gaussianLLH(y0(2)-rho, sigma_rho_sq)-log(y0(2)); // r
        y0(3) = RNG.randu()*2*arma::datum::pi;
        llhAll(4) = -LOG2PI; // theta
        y0 = transVinv(y0); //Translate back to y space
    }
}

/**
 * @param[in] xt The current observed position
 * @param[in] wt The current observed localization error
 * @param[in] y0 The previous estimate of the position.  Size [4,Nparticles].  Each particle is a column
 * @param[in] z0 The previous estimate of the state. Size Nparticles.
 * @param[in] t  The time elapsed since previous measurement
 * @param[out] yt The new estimate of position. Size [4,Nparticles]
 * @param[out] zt The new estimate of the state. Size [Nparticles]
 * @param[out] llh The log-likelihood of choosing [yt,zt].  This is llhQ of the proposal distribution.
 */
void PairInteractionSMC::sampleQcomb(const VecT &xt, const VecT &wt, const MatT &y0, const VecT &z0, 
                                     double t, MatT &yt, VecT &zt, VecT &llh) const
{
    IdxT N=z0.n_elem;

    UnscentedTransform ut(0,.5,2,4);
    double varA = 2*t*DA;
    double varB = 2*t*DB;
    double E = exp(-gamma*t);
    double varC = 2*t*DAB;
    double varR = sigma_rho_sq*(1-square(E));
    double sigmaR = sqrt(varR);
    double varPhi = std::min(2*arma::datum::pi, 2*Dphi*t);

    VecT yt_var_free = {varA,varA,varB,varB};
    VecT vt_var_bound = {varC,varC,varR,varPhi};
    VecT vt_var_c = {varC,varC};
    auto vt_cov_bound = arma::diagmat(vt_var_bound);
    double p_unbind = 1-exp(-koff*t);
    VecT yt_var(4); //Est. yt variance using gaussian bayes' rule
    
    VecT yt_ut_mean(4); //Est. Mean for yt given z0=Bound using UT
    MatT yt_ut_cov(4,4); //Est. covariance for yt given z0=Bound using UT
    VecT yt_gb_mean(4); //Est. Mean for yt using Gaussian Bayes rule
    VecT yt_gb_var(4); //Est. variance for yt using Gaussian Bayes rule
    VecT ct_mean(2),ct_var(2);
    VecT xct = {.5*(xt(0)+xt(2)), .5*(xt(1)+xt(3))};
    VecT wtsq = wt%wt;
    VecT xct_var = {.25*(wtsq(0)+wtsq(2)), .25*(wtsq(1)+wtsq(3))};
    
    //These parameters are for the estimation of angle phi from vector x_a - x_b.  The vector from b to a.
    //Only need to estimate VM parameters once for all particles here since the estimate is based on xt only.
    //For most realistic frame rates phi angle distribution is mainly dependent on xt more so that y_t-1 since
    //The lateral and angular diffusions are large frame-to-frame.
    VecT delta_xt = {xt(0)-xt(2), xt(1)-xt(3)}; //xt_a - xt_b: vec from b to a
    VecT sigma_delta_xt = {sqrt(wtsq(0)+wtsq(2)),sqrt(wtsq(1)+wtsq(3))}; //w_a^2 + w_b^2.  Apparent variance of a-b from a.
    GaussianAngleDist phiDist(delta_xt,sigma_delta_xt(0),sigma_delta_xt(1));
    for(IdxT n=0; n<N;n++) {
        //Step1: sample zt ~ p(zt | xt, wt, z0, y0)
        //Step1a: p(xt | zt=F, wt, y0) = \int p(xt | yt, wt) p(yt | zt, z0, y0) = N(
        //p(xt|zt=F,wt,y0) ~ N(xt|y0, varF + varObs)
        double llh_xt_free = 0;
        for(IdxT k=0; k<2; k++) llh_xt_free += gaussianLLH(xt(k)-y0(k,n),varA+ wtsq(k));
        for(IdxT k=2; k<4; k++) llh_xt_free += gaussianLLH(xt(k)-y0(k,n),varB+ wtsq(k));
        double p_xt_free = exp(llh_xt_free);
        
        //p(xt|zt=B,wt,y0)  Use unscented transform to estimate p(yt|zt=B,y0) ~ N(yt|yt_mean,yt_cov).
        VecT v0 = transV(y0.col(n));
        VecT vt_mean_bound = v0;
        vt_mean_bound(2) = E*vt_mean_bound(2) + (1-E)*rho;
        ut.transform([](const VecT &v) {return transVinv(v);}, vt_mean_bound, vt_cov_bound, yt_ut_mean, yt_ut_cov);
        double llh_xt_bound = multiGaussianLLH(xt, yt_ut_mean, (yt_ut_cov+arma::diagmat(wtsq)).eval());
        double p_xt_bound = exp(llh_xt_bound);
        //p (zt=F| z0, y0)
        double p_zt_free = (z0(n)==0.) ? rd.survivalProb(v0(2),t) : p_unbind;
        double p_zt_bound = 1.-p_zt_free;
        // Uses Bayes' rule to get p(zt | xt,wt,y0,z0) = p(xt | zt, wt, y0, z0)*p(zt|y0,z0) / Norm
        double p_zt_free_given_xt = p_xt_free*p_zt_free/(p_xt_free*p_zt_free + p_xt_bound*p_zt_bound);
        std::cout<<"\n[Particle:"<<n<<"] -- z0:"<<z0(n)<<" v0: "<<v0.t()<<" vt_mean:"<<vt_mean_bound.t()<<"\n";
        std::cout<<"llh(xt|free):"<<llh_xt_free<<" llh(xt|bound):"<<llh_xt_bound<<"\n";
        std::cout<<"p(xt|free):"<<p_xt_free<< " p(zt=f):"<<p_zt_free<<" p(xt|bound):"<<p_xt_bound<< " p(zt=b):"<<p_zt_bound<<"\n";
        std::cout<<"p(zt=0| x_t)="<<p_zt_free_given_xt;
        double rnd = RNG.randu();
        if(rnd<p_zt_free_given_xt){
            //Next state is free
//             std::cout<<" rnd:"<<rnd<<" sampled->zt=0\n";
            zt(n)=0.;
            llh(n)=log(p_zt_free_given_xt);
            //Sample ys using bayes' rule for gaussians since yt_var_free and wtsq are diagonal so is yt_var
            gaussianBayes(y0.col(n), yt_var_free, xt, wtsq, yt_gb_mean, yt_gb_var);
//             std::cout<<"Free prediction using Gaussian Bayes' rule>>>\n";
//             std::cout<<"Y0: "<<y0.col(n).t()<<"yt sigma: "<<arma::sqrt(yt_var_free).eval().t()<<"xt: "<<xt.t()<<"xt sigma: "<<arma::sqrt(wtsq).eval().t()<<"Est yt_mean: "<<yt_gb_mean.t()<<"Est yt sigma: "<<arma::sqrt(yt_gb_var).eval().t();
            yt.col(n) = sampleGaussian(yt_gb_mean, arma::sqrt(yt_gb_var).eval()); //sample new yt
            double freeLLH_delta = gaussianLLH(yt.col(n),yt_gb_mean, yt_gb_var);
//             std::cout<<"Sampled yt:"<<yt.col(n).t();
//             std::cout<<"q(zt|xt,wt,y0,z0):"<<llh(n)<<" q(yt|zt,xt,wt,y0,z0): "<<freeLLH_delta<<" LLH Q total:"<<freeLLH_delta+llh(n)<<std::endl;
            llh(n) += freeLLH_delta; //record llh
            
        } else {
            //Next state is bound
//             std::cout<<" rnd:"<<rnd<<" sampled->zt=1\n";
            zt(n)=1.;
            double llh_zt = log(1-p_zt_free_given_xt);
            llh(n)=llh_zt;
            //Sample cs using bayes' rule for gaussians
            gaussianBayes(v0.subvec(0,1), vt_var_c, xct, xct_var, ct_mean, ct_var);
//             std::cout<<"GaussianBayes: c0:["<<v0(0)<<","<<v0(1)<<"] varC:"<<varC<<" xct:["<<xct(0)<<","<<xct(1)<<"]: xctVar:["<<xct_var(0)<<","<<xct_var(1)
//                      <<"] ct_mean:["<<ct_mean(0)<<","<<ct_mean(1)<<"] ctVar:["<<ct_var.t();
            VecT ct = sampleGaussian(ct_mean, arma::sqrt(ct_var).eval());
            yt(0,n)=ct(0);
            yt(1,n)=ct(1);
            llh(n) += gaussianLLH(ct, ct_mean, ct_var);
            //Sample r using OU model
//             std::cout<<"LLH zt=1:"<<log(1-p_zt_free_given_xt)<<" LLH Qc:"<<gaussianLLH(ct, ct_mean, ct_var);
            if(z0(n)==0.){ //Were free, now bound, sample from equilibrium dist
                yt(2,n) = sampleGaussian(rho,sigma_rho);
                llh(n) += gaussianLLH(yt(2,n)-rho,sigma_rho_sq);
//                 std::cout<<" LLH Qr:"<<gaussianLLH(yt(2,n)-rho,sigma_rho_sq);
            } else {
                double r_mean = E*v0(2) + (1-E)*rho;
                yt(2,n) = sampleGaussian(r_mean, sigmaR);
                llh(n) += gaussianLLH(yt(2,n)-r_mean, varR);
//                 std::cout<<" LLH Qr:"<<gaussianLLH(yt(2,n)-r_mean, varR);
            }
            //Sample phi using observation probability for x_b-x_a 
            yt(3,n) = phiDist.sample();
            llh(n) += log(phiDist.computePDF(yt(3,n)));
//                std::cout<<" LLH Qphi:"<<vonMisesLLH(yt(3,n)-mu, kappa)<<"\n";
//             std::cout<<"q(zt|xt,wt,y0,z0):"<<llh_zt<<" q(yt|zt,xt,wt,y0,z0): "<<llh(n)-llh_zt<<" LLH Q total:"<<llh(n)<<std::endl;
            //Translate back to y-space
//             std::cout<<" v0: "<<v0.t();
//             std::cout<<" vt: "<<yt.col(n).t();
            yt.col(n) = transVinv(yt.col(n));
//             std::cout<<" y0: "<<y0.col(n).t();
//             std::cout<<" yt: "<<yt.col(n).t();
//             std::cout<<" xt: "<<xt.t();
//             std::cout<<"delt:"<<((yt.col(n)-xt)).eval().t();
//             std::cout<<" wt: "<<wt.t();
        }
    }
//     std::cout<<"Number bound(0): "<<arma::sum(z0)<<"/"<<N<<" bound(t): "<<arma::sum(zt)<<"/"<<N<<"\n";
//     std::cout<<"LLH: "<<llh.t()<<"State: "<<zt.t()<<"\n";
}

void PairInteractionSMC::sampleQobs(const VecT &xt, const VecT &wt, const MatT &y0, const VecT &z0, double t, MatT &yt, VecT &zt, VecT &llh) const
{
    IdxT N=yt.n_cols;
    sampleH_inv(xt,wt,yt,llh);
    VecT ztFree(N,arma::fill::zeros);
    VecT ztBound(N,arma::fill::ones);
    VecT llhGfree = llhG_vec(yt,ztFree,y0,z0,t);
    VecT llhGbound = llhG_vec(yt,ztBound,y0,z0,t);
    VecT pFree = logNormalizeBinary(llhGfree, llhGbound);
//     std::cout<<"llhGfree: "<<llhGfree.t()<<"\n";
//     std::cout<<"llhGbound: "<<llhGbound.t()<<"\n";
//     std::cout<<"pFree: "<<pFree.t()<<"\n";
    for(IdxT n=0; n<N; n++){
        if(RNG.randu() <= pFree(n)) {
            zt(n)=0;
            llh(n)+=log(pFree(n));
        } else {
            zt(n)=1;
            llh(n)+=log(1-pFree(n));
        }
    }
//     std::cout<<"LLh: "<<llh.t()<<"\n";
//     std::cout<<"#Free z0:"<<arma::sum(z0)<<" #Free zt:"<<arma::sum(zt)<<"\n";
}

void PairInteractionSMC::sampleH_inv(const VecT &xt, const VecT &wt, MatT &yt, VecT &llh) const
{
    IdxT N=yt.n_cols;
    VecT wsquare = arma::square(wt);
    llh.zeros();
    if(llh.n_elem != N) throw std::logic_error("uninitialized LLH");
    for(IdxT n=0; n<N;n++) for(IdxT k=0; k<4; k++) {
        yt(k,n) = xt(k)+RNG.randn()*wt(k);
        llh(n) += gaussianLLH(xt(k)-yt(k,n), wsquare(k));
    }
}

void PairInteractionSMC::sampleG(const MatT &y0, const VecT &z0, double t, MatT &yt, VecT &zt, VecT &llh) const
{
    IdxT N=z0.n_elem;
    double varA = 2*t*DA;
    double varB = 2*t*DB;
    double sqrtVarA = sqrt(varA);
    double sqrtVarB = sqrt(varB);
    double varC = 2*t*DAB;
    double sqrtVarC = sqrt(varC);
    double unbindProb = 1-exp(-koff*t);
    double log_unbindProb = log(unbindProb);
    double log_1m_unbindProb = -koff*t;
    double E = exp(-gamma*t);
    double varR = sigma_rho_sq*(1-square(E));
    double sqrtVarR = sqrt(varR);
    double phiKappa = 1./(2*Dphi*t);
    VecT v0(4),vt(4);
    double diff_sigma = sqrt(varA+varB);
    GaussianAngleDist phiDist;
    phiDist.setSigma(diff_sigma,diff_sigma); 
    for(IdxT n=0;n<N;n++){
        llh(n)=0;
        if(z0(n)==0.) { //was free
            // log(p(z=F | z0=F, y0))
            double survivalProb = rd.survivalProb(distance(y0.col(n)),t);
            if(RNG.randu() <= survivalProb) { //F -> F
                //Set state (Z)
                zt(n) = 0.; //Free
                llh(n) = log(survivalProb);
                //Sample (Y) p(yt | y0, zt=Free)
                for(auto k=0;k<2;k++) {
                    yt(k,n) = sampleGaussian(y0(k,n), sqrtVarA);
                    llh(n) += gaussianLLH(yt(k,n)-y0(k,n), varA);
                }
                for(auto k=2;k<4;k++) {
                    yt(k,n) = sampleGaussian(y0(k,n), sqrtVarB);
                    llh(n) += gaussianLLH(yt(k,n)-y0(k,n),varB);
                }
            } else { //F->B Capture 
                //Set state (Z)
                zt(n) = 1.; //Bound
                llh(n) = log(1-survivalProb);
                //Sample (V)  p(vt | v0, zt=B)
                VecT v0 = transV(y0.col(n));
                for(IdxT k=0; k<2; k++) {
                    vt(k) = sampleGaussian(v0(k), sqrtVarC);
                    llh(n) += gaussianLLH(vt(k)-v0(k), varC);
                }
                vt(2) = sampleGaussian(rho,sigma_rho);
                llh(n) += gaussianLLH(vt(2)-rho, sigma_rho_sq)-log(vt(2)); // r
                
                // draw from
                VecT diff_mean = y0.submat(0,n,1,n) - y0.submat(2,n,3,n); //y0_a - y0_b: vec from b to a
                phiDist.setMean(diff_mean);
                vt(3) = phiDist.sample();
                llh(n) += log(phiDist.computePDF(vt(3)));

                yt.col(n) = transVinv(vt); // Translate back to y-space
            }
        } else { //was bound
            if(RNG.randu() <= unbindProb) { //B -> F (Unbinding)
                //Set state (Z)
                zt(n) = 0.; //Free
                llh(n) = log_unbindProb;
                //Sample (Y) p(yt | y0, zt=Free)
                for(auto k=0;k<2;k++) {
                    yt(k,n) = sampleGaussian(y0(k,n), sqrtVarA);
                    llh(n) += gaussianLLH(yt(k,n)-y0(k,n),varA);
                }
                for(auto k=2;k<4;k++) {
                    yt(k,n) = sampleGaussian(y0(k,n), sqrtVarB);
                    llh(n) += gaussianLLH(yt(k,n)-y0(k,n),varB);
                }
            } else { //B -> B
                //Set state (Z)
                zt(n) = 1.; //Bound
                llh(n) = log_1m_unbindProb;
                
                //Sample (V)  p(vt | v0, zt=B)
                VecT v0 = transV(y0.col(n));
                //Sample C
                for(IdxT k=0; k<2; k++) {
                    vt(k) = sampleGaussian(v0(k), sqrtVarC);
                    llh(n) += gaussianLLH(vt(k)-v0(k),varC);
                }
                //Sample R
                double muR = rho*(1-E) + v0(2)*E;
                vt(2) = sampleGaussian(muR, sqrtVarR);
                llh(n) += gaussianLLH(vt(2)-muR, varR) - log(vt(2));
                //Sample Phi
                vt(3) = sampleVonMises(v0(3),phiKappa);
                llh(n) += vonMisesLLH(vt(3)-v0(3), phiKappa);
                yt.col(n) = transVinv(vt); // Translate back to y-space
            }
        }
    }
}


/*
 * 
 * We could make this more efficient by eliminating the sqrt calls.  Evaluate later to see if this is necessary
 */
void PairInteractionSMC::sampleG(const VecT &y0, double z0, double t, VecT &yt, double &zt, VecT &llhAll) const
{
    if(z0==0.) { //was free
        double survivalProb = rd.survivalProb(distance(y0),t);  // log(p(z=F | z0=F, y0))
        if(RNG.randu() <= survivalProb) { //F -> F
            //Set state (Z)
            zt = 0.; //Free
            llhAll(0) = log(survivalProb);

            double varA = 2*t*DA;
            double varB = 2*t*DB;
            double sqrtVarA = sqrt(varA);
            double sqrtVarB = sqrt(varB);
            for(auto n=0;n<2;n++) {
                yt(n) = y0(n) + RNG.randn()*sqrtVarA;
                llhAll(n+1) = gaussianLLH(yt(n)-y0(n),varA);
            }
            for(auto n=2;n<4;n++) {
                yt(n) = y0(n) + RNG.randn()*sqrtVarB;
                llhAll(n+1) = gaussianLLH(yt(n)-y0(n),varB);
            }
        } else { //F->B Capture 
            //Set state (Z)
            zt = 1.; //Bound
            llhAll(0) = log(1-survivalProb);
 
            VecT v0 = transV(y0);
            double varC = 2*t*DAB;
            double sqrtVarC = sqrt(varC);
            for(IdxT n=0; n<2; n++) {
                yt(n) = v0(n) + RNG.randn()*sqrtVarC;
                llhAll(1+n) = gaussianLLH(yt(n)-v0(n),varC);
            }
            yt(2) = RNG.randn()*sigma_rho+rho;
            llhAll(3) = gaussianLLH(yt(2)-rho, sigma_rho_sq)-log(yt(2)); // r
            
            // draw from p_VM(phi | mu, kappa_est(c0, SigmaA+SigmaB))
            double diff_sigma = sqrt(2*t*(DA+DB));
            VecT diff_mean = y0.subvec(0,1) - y0.subvec(2,3); //y0_a - y0_b: vec from b to a
            GaussianAngleDist phiDist(diff_mean,diff_sigma,diff_sigma); 
            yt(3) = phiDist.sample();
            llhAll(4) = log(phiDist.computePDF(yt(3)));
            yt = transVinv(yt); // Translate back to y-space
        }
    } else { //was bound
        double unbindProb = 1-exp(-koff*t);
        if(RNG.randu() <= unbindProb) { //B -> F (Unbinding)
            //Set state (Z)
            zt = 0.; //Free
            llhAll(0) = log(unbindProb);
            
            double varA = 2*t*DA;
            double varB = 2*t*DB;
            double sqrtVarA = sqrt(varA);
            double sqrtVarB = sqrt(varB);
            for(auto n=0;n<2;n++) {
                yt(n) = y0(n) + RNG.randn()*sqrtVarA;
                llhAll(n+1) = gaussianLLH(yt(n)-y0(n),varA);
            }
            for(auto n=2;n<4;n++) {
                yt(n) = y0(n) + RNG.randn()*sqrtVarB;
                llhAll(n+1) = gaussianLLH(yt(n)-y0(n),varB);
            }
        } else { //B -> B
            //Set state (Z)
            zt = 1.; //Bound
            llhAll(0) = log(1-unbindProb);
            
            VecT v0 = transV(y0);
            double varC = 2*t*DAB;
            double sqrtVarC = sqrt(varC);
            for(IdxT n=0; n<2; n++) {
                yt(n) = v0(n) + RNG.randn()*sqrtVarC;
                llhAll(1+n) = gaussianLLH(yt(n)-v0(n),varC);
            }
            double E = exp(-gamma*t);
            double muR = rho*(1-E) + v0(2)*E;
            double varR = sigma_rho_sq*(1-square(E));
            yt(2) = muR + sqrt(varR)*RNG.randn();
            llhAll(3) = gaussianLLH(yt(2)-muR, varR) - log(yt(2)); // r
            double phiKappa = 1./(2*Dphi*t);
            yt(3) = sampleVonMises(v0(3),phiKappa);
            llhAll(4) = vonMisesLLH(yt(3)-v0(3), phiKappa); // theta
            yt = transVinv(yt); // Translate back to y-space
        }
    }
}

PairInteractionSMC::VecT PairInteractionSMC::llhG0_all(const VecT &y0, double z0) const
{
    VecT llhAll(5);
    if(z0==0.){ //Free
        llhAll(0) = log_state_prior(0);
        llhAll(1) = .5*pos_prior->llh_pos(y0.subvec(0,1));
        llhAll(2) = llhAll(1);
        llhAll(3) = .5*pos_prior->llh_pos(y0.subvec(2,3));
        llhAll(4) = llhAll(3);
    } else { //Bound
        VecT v0 = transV(y0);
        llhAll(0) = log_state_prior(1);
        llhAll(1) = .5*pos_prior->llh_pos(v0.subvec(0,1)); //c_x, c_y
        llhAll(2) = llhAll(1);
        llhAll(3) = gaussianLLH(v0(2)-rho, sigma_rho_sq)-log(v0(2)); // r
        llhAll(4) = -LOG2PI; // theta
    }
    return llhAll;
}


VecT PairInteractionSMC::llhG_vec(const MatT &yt, const VecT &zt, const MatT &y0, const VecT &z0, double t) const
{
    IdxT N = yt.n_cols;
    double varA = 2*t*DA;
    double varB = 2*t*DB;
    double varC = 2*t*DAB;
    double E = exp(-gamma*t);
    double varR = sigma_rho_sq*(1-square(E));
    double phiKappa = 1./(2*Dphi*t);
    double diff_sigma = sqrt(2*t*(DA+DB)); //variance of mutual diffusion in each dimension
    GaussianAngleDist phiDist; 
    phiDist.setSigma(diff_sigma,diff_sigma);
    
    double log_unbindProb = log1p(-exp(-koff*t));
    double log_1m_unbindProb = -koff*t;
    VecT v0(4),vt(4);
    VecT llh(N);
    for(IdxT n=0;n<N;n++){
        if(zt(n)==0){ //Is now free
            if(z0(n)==0) llh(n) = rd.survivalLogProb(distance(y0.col(n)),t);  // Prob to survive capture, log(p(z=F | z0=F, y0))
            else         llh(n) = log_unbindProb;  //Prob to unbind log(p(z=F | z0=B, y0))
            for(IdxT k=0; k<2; k++) llh(n) += gaussianLLH(yt(k,n)-y0(k,n),varA);
            for(IdxT k=2; k<4; k++) llh(n) += gaussianLLH(yt(k,n)-y0(k,n),varB);
        } else { //Is now bound
            v0 = transV(y0.col(n));
            vt = transV(yt.col(n));
            if(z0(n)==0) { // F-> B
                llh(n) = rd.captureLogProb(v0(2),t);  // log(p(z=B | z0=F, y0))
                for(IdxT k=0; k<2; k++) llh(n) += gaussianLLH(vt(k)-v0(k),varC);
                llh(n) += gaussianLLH(vt(2)-rho, sigma_rho_sq)-log(vt(2)); // r
                VecT diff_mean = y0.submat(0,n,1,n)-y0.submat(2,n,3,n);//y0_a - y0_b: vec from b to a
                phiDist.setMean(diff_mean);
                llh(n) += log(phiDist.computePDF(vt(3)));
            } else { // B-> B
                llh(n) = log_1m_unbindProb; //Prob to remain bound log(p(z=B | z0=B, y0))
                for(IdxT k=0; k<2; k++) llh(n) += gaussianLLH(vt(k)-v0(k),varC);
                double muR = rho*(1-E)+v0(2)*E;
                llh(n) += gaussianLLH(vt(2)-muR, varR) - log(vt(2)); // r
                llh(n) += vonMisesLLH(vt(3)-v0(3), phiKappa); // theta
            }
        }
    }
    return llh;
}

VecT PairInteractionSMC::llhH_vec(const VecT &x, const VecT &w, const MatT &y) const
{
    IdxT N = y.n_cols;
    VecT wsquare = w%w;
    VecT llh(N,arma::fill::zeros);
    for(IdxT n=0; n<N; n++) {
        for(IdxT k=0; k<4;k++) llh(n) += gaussianLLH(x(k)-y(k,n),wsquare(k));
    }
    return llh;
}


PairInteractionSMC::VecT PairInteractionSMC::llhG_all(const VecT &yt, double zt, const VecT &y0, double z0, double t) const
{
    VecT llhAll(5);
    if(z0==0.){ //was free
        if(zt==0.) { // F-> F
            llhAll(0) = rd.survivalLogProb(distance(y0),t);  // log(p(z=F | z0=F, y0))
            double varA = 2*t*DA;
            for(IdxT n=0; n<2; n++) llhAll(1+n) = gaussianLLH(yt(n)-y0(n),varA);
            double varB = 2*t*DB;
            for(IdxT n=2; n<4; n++) llhAll(1+n) = gaussianLLH(yt(n)-y0(n),varB);
        } else { // F-> B
            VecT v0 = transV(y0);
            VecT v = transV(yt);
            llhAll(0) = rd.captureLogProb(v0(2),t);  // log(p(z=B | z0=F, y0))
            double varC = 2*t*DAB;
            for(IdxT n=0; n<2; n++) llhAll(1+n) = gaussianLLH(v(n)-v0(n),varC);
            llhAll(3) = gaussianLLH(v(2)-rho, sigma_rho_sq)-log(v(2)); // r
            double diff_sigma = sqrt(2*t*(DA+DB));
            VecT diff_mean = y0.subvec(0,1)-y0.subvec(2,3); //y0_a - y0_b: vec from b to a
            GaussianAngleDist phiDist(diff_mean,diff_sigma,diff_sigma); 
            llhAll(4) = log(phiDist.computePDF(v(3)));
        }
    } else { //was bound
        if(zt==0.) { // B-> F
            llhAll(0) = log1p(-exp(-koff*t)); //Prob to unbind log(p(z=F | z0=B, y0))
            double varA = 2*t*DA;
            for(IdxT n=0; n<2; n++) llhAll(1+n) = gaussianLLH(yt(n)-y0(n),varA);
            double varB = 2*t*DB;
            for(IdxT n=2; n<4; n++) llhAll(1+n) = gaussianLLH(yt(n)-y0(n),varB);
        } else { // B-> B
            llhAll(0) = -koff*t; //Prob to remain bound log(p(z=B | z0=B, y0))
            VecT v0 = transV(y0);
            VecT v = transV(yt);
            double varC = 2*t*DAB;
            for(IdxT n=0; n<2; n++) llhAll(1+n) = gaussianLLH(v(n)-v0(n),varC);
            double E = exp(-gamma*t);
            double muR = rho*(1-E)+v0(2)*E;
            double varR = sigma_rho_sq*(1-square(E));
            llhAll(3) = gaussianLLH(v(2)-muR, varR) - log(v(2)); // r
            double phiKappa = 1./(2*Dphi*t);
            llhAll(4) = vonMisesLLH(v(3)-v0(3), phiKappa); // theta
        }
    }
    return llhAll;
}


void PairInteractionSMC::llhG0_parallel(const MatT &y0, const VecT &z0, VecT &llh) const
{
    IdxT N = z0.n_elem;
    if(y0.n_rows!=4 || y0.n_cols!=N) throw std::invalid_argument("y0 must be 4xN where N is length of z0");
    if(llh.n_elem!=N) throw std::invalid_argument("llh should have same number of elements as z0");
    #pragma omp parallel for
    for(IdxT n=0; n<N; n++) llh(n) = llhG0(y0.col(n),z0(n));
}

void PairInteractionSMC::llhG_parallel(const MatT &yt, const VecT &zt, const MatT &y0, const VecT &z0, double t, VecT &llh) const 
{
    IdxT N = z0.n_elem;
    if(y0.n_rows!=4 || y0.n_cols!=N) throw std::invalid_argument("y0 must be 4xN where N is length of z0");
    if(zt.n_elem!=N) throw std::invalid_argument("zt be same length as z0");
    if(yt.n_rows!=4 || yt.n_cols!=N) throw std::invalid_argument("yt must be 4xN where N is length of z0");
    if(llh.n_elem!=N) throw std::invalid_argument("llh should have same number of elements");
    #pragma omp parallel for
    for(IdxT n=0; n<N; n++) llh(n) = llhG(yt.col(n),zt(n), y0.col(n), z0(n),t);
}

void PairInteractionSMC::llhH_parallel(const MatT &x, const MatT &w, const MatT &y, VecT &llh) const
{
    IdxT N = x.n_cols;
    if(x.n_rows!=4) throw std::invalid_argument("x must be 4xN");
    if(w.n_rows!=4 || w.n_cols!=N) throw std::invalid_argument("w must be 4xN");
    if(y.n_rows!=4 || y.n_cols!=N) throw std::invalid_argument("y must be 4xN");
    if(llh.n_elem!=N) throw std::invalid_argument("llh should have same number of elements as x");
    #pragma omp parallel for
    for(IdxT n=0; n<N; n++) llh(n) = llhH(x.col(n),w.col(n),y.col(n));
}



GaussianPositionPrior::GaussianPositionPrior(const VecT &_mean, const MatT &_cov)
    : mean(_mean), cov(_cov)
{
    if(mean.n_elem!=2) throw std::invalid_argument("mean point should be size:[2]");
    if(cov.n_rows!=2 || cov.n_cols!=2) throw std::invalid_argument("covariance matrix should be size:[2,2]");
    try {
        chol_cov = arma::chol(cov);
    } catch (std::runtime_error &e) {
        throw std::runtime_error("Invalid prior covariance matrix");
    }
    VecT lambda = arma::eig_sym(cov);
    double omega = 2*arma::sum(lambda);
    double q_sq = lambda(1)/lambda(0);
    double q = sqrt(q_sq);
    double q_sq_p1 = 1+q_sq;
    double denom = 4*q_sq*omega;
    A = q_sq_p1/(q*omega); // (1+q^2)/(q*omega)
    B = -q_sq_p1*q_sq_p1/denom; // -(1+q^2)^2/(4*q^2*omega)
    C = (1-q_sq*q_sq)/denom; // (1-q^4)/(4*q^2*omega);
}

GaussianPositionPrior::ParamsT GaussianPositionPrior::getParams() const
{
    ParamsT params;
    params["priorPosType"] = double(DistType::GAUSSIAN);
    params["priorMeanX"] = mean(0);
    params["priorMeanY"] = mean(1);
    params["priorCov11"] = cov(0,0);
    params["priorCov12"] = cov(0,1);
    params["priorCov21"] = cov(1,0);
    params["priorCov22"] = cov(1,1);
    return params;
}

RectangularPositionPrior::RectangularPositionPrior(const arma::vec &_rect) : rect(_rect)
{
    width=rect(2)-rect(0);
    height=rect(3)-rect(1);
    if(width>height) std::swap(width,height); // Ensure width<=height
    area = width * height;
    llh = -log(area); //log likelihood of a position
    width_sq = width*width;
    height_sq = height*height;
}

double RectangularPositionPrior::p_r0(double r) const
{
    double r_sq = r*r;
    if(r <= width_sq) {
        return r_sq + arma::datum::pi*area-2*r*(height-width);
    } else if(r <= height_sq) {
        return 2*area*asin(width/r) + 2*height*sqrt(r_sq-width_sq) - width_sq - 2*r*height;
    } else {
        return 2*area*asin(width/r) + 2*height*sqrt(r_sq-width_sq) - width_sq 
               + 2*area*asin(height/r) + 2*width*sqrt(r_sq-height_sq) - height_sq 
               - arma::datum::pi*area-r_sq;
    }
}

RectangularPositionPrior::ParamsT RectangularPositionPrior::getParams() const
{
    ParamsT params;
    params["priorPosType"] = double(DistType::RECTANGULAR);
    params["priorXmin"] = rect(0);
    params["priorYmin"] = rect(1);
    params["priorXmax"] = rect(2);
    params["priorYmax"] = rect(3);
    return params;
}


CircularPositionPrior::CircularPositionPrior(const arma::vec &_center, double _radius) : center(_center), radius(_radius)
{
    area = arma::datum::pi*radius*radius; //pi*r^2
    llh = -log(area);  //log likelihood of a position
}

double CircularPositionPrior::p_r0(double r) const
{
    double omega = r/(2*radius);
    return (4*r/area)*(acos(omega)-omega*sqrt(1-omega*omega));
}


CircularPositionPrior::ParamsT CircularPositionPrior::getParams() const
{
    ParamsT params;
    params["priorPosType"] = double(DistType::CIRCULAR);
    params["priorCenterX"] = center(0);
    params["priorCenterY"] = center(1);
    params["priorRadius"] = radius;
    return params;
}
