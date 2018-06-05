/** @file PairInteractionSMC_Iface.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief The class declaration and inline and templated functions for PairInteractionSMC.
 */
#include "PairInteractionSMC_Iface.h"
#include "Handle.h"
using namespace pair_int;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    PairInteractionSMC_Iface iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}

PairInteractionSMC_Iface::PairInteractionSMC_Iface() 
    : Mex_Iface("PairInteractionSMC_Iface") 
{
    methodmap["addData"]=boost::bind(&PairInteractionSMC_Iface::addData, this);
    methodmap["clearData"]=boost::bind(&PairInteractionSMC_Iface::clearData, this);
    methodmap["getData"]=boost::bind(&PairInteractionSMC_Iface::getData, this);
    methodmap["setParams"]=boost::bind(&PairInteractionSMC_Iface::setParams, this);
    methodmap["getParams"]=boost::bind(&PairInteractionSMC_Iface::getParams, this);
    methodmap["setPriorState"]=boost::bind(&PairInteractionSMC_Iface::setPriorState, this);
    methodmap["setPriorPositionGaussian"]=boost::bind(&PairInteractionSMC_Iface::setPriorPositionGaussian, this);
    methodmap["setPriorPositionRectangular"]=boost::bind(&PairInteractionSMC_Iface::setPriorPositionRectangular, this);
    methodmap["setPriorPositionCircular"]=boost::bind(&PairInteractionSMC_Iface::setPriorPositionCircular, this);
    methodmap["runParticleFilter"]=boost::bind(&PairInteractionSMC_Iface::runParticleFilter, this);
    methodmap["obsLLH"]=boost::bind(&PairInteractionSMC_Iface::obsLLH, this);
    methodmap["sampleParticle"]=boost::bind(&PairInteractionSMC_Iface::sampleParticle, this);
    methodmap["sampleParticleLLH"]=boost::bind(&PairInteractionSMC_Iface::sampleParticleLLH, this);
    methodmap["getAllParticles"]=boost::bind(&PairInteractionSMC_Iface::getAllParticles, this);
    methodmap["computeLLH"]=boost::bind(&PairInteractionSMC_Iface::computeLLH, this);
    methodmap["computeLLH_debug"]=boost::bind(&PairInteractionSMC_Iface::computeLLH_debug, this);
    methodmap["simulate"]=boost::bind(&PairInteractionSMC_Iface::simulate, this);
    methodmap["llhG0"]=boost::bind(&PairInteractionSMC_Iface::llhG0, this);
    methodmap["llhG"]=boost::bind(&PairInteractionSMC_Iface::llhG, this);
    methodmap["llhH"]=boost::bind(&PairInteractionSMC_Iface::llhH, this);
}

void PairInteractionSMC_Iface::objConstruct() 
{
    // [in] params: parameters struct
    // [in] data: [optional] cell-array of obersevation matricies
    checkMinNumArgs(1,1);
    checkMaxNumArgs(1,2);
    auto params=getDoubleStruct();
    PairInteractionSMC *obj;
    if (nrhs==2) {
        auto data = getMatVector<double>();
        obj = new PairInteractionSMC(params,data);
    } else {
        obj = new PairInteractionSMC(params);
    }
    outputMXArray(Handle<PairInteractionSMC>::makeHandle(obj));
}

void PairInteractionSMC_Iface::objDestroy()
{
    checkNumArgs(0,1);
    Handle<PairInteractionSMC>::destroyObject(rhs[0]);
}

void PairInteractionSMC_Iface::getObjectFromHandle(const mxArray *mxhandle)
{
    obj=Handle<PairInteractionSMC>::getObject(mxhandle);
}

/* Exposed methods */

void PairInteractionSMC_Iface::addData()
{
    // Append new data to existing data
    // [in] data: cell-array of obersevation matricies
    // [out] Ndata: int - total number of observations in database
    checkNumArgs(1,1);
    outputInt(obj->addData(getMatVector<double>()));
}

void PairInteractionSMC_Iface::clearData()
{
    //Clear all data
    checkNumArgs(0,0);
    obj->clearData();
}

void PairInteractionSMC_Iface::getData()
{
    //Retrieve all data
    // [out] data: cell-array of obersevation matricies
    checkNumArgs(1,0);
    outputMatCellArray(obj->getData());
}

void PairInteractionSMC_Iface::setParams()
{
    // Set parameters
    // [in] params: struct of parameters to set
    checkNumArgs(0,1);
    obj->setParams(getDoubleStruct());
}

void PairInteractionSMC_Iface::getParams()
{
    // Get all parameters
    // [out] params: struct of parameters
    checkNumArgs(1,0);
    outputStatsToStruct(obj->getParams());
}

void PairInteractionSMC_Iface::setPriorState()
{
    // Set the prior to Gaussian
    // [in] state_prior - size:[2] - [pFree, pBound] must sum to 1    
    checkNumArgs(0,1);
    obj->setPriorState(getDVec());
}

void PairInteractionSMC_Iface::setPriorPositionGaussian()
{
    // Set the position prior to Gaussian
    // [in] center - [x,y] position of center of Gaussian
    // [in] cov - 2x2 covariance matrix for Gaussian
    checkNumArgs(0,2);
    auto center = getDVec();
    auto cov = getDMat();
    obj->setPriorPositionGaussian(center, cov);
}

void PairInteractionSMC_Iface::setPriorPositionRectangular()
{
    // Set the position prior to Gaussian
    // [in] rect - size:[4] - [min_x min_y max_x max_y];
    checkNumArgs(0,1);
    obj->setPriorPositionRectangular(getDVec());
}

void PairInteractionSMC_Iface::setPriorPositionCircular()
{
    // Set the position prior to Circular (Disk)
    // [in] center - [x,y] position of center of circle
    // [in] radius - scalar - radius of circle
    checkNumArgs(0,2);
    auto center = getDVec();
    double radius = getDouble();
    obj->setPriorPositionCircular(center,radius);
}
/*
void PairInteractionSMC_Iface::getPrior()
{
    // Get the state and position prior in a structure format
    checkNumArgs(1,0);
    
}*/

void PairInteractionSMC_Iface::runParticleFilter()
{
    // Run the particle filter
    // [in] Nparticles- int - number of particles to simulates
    // [in] propType- string - ["Transition", "Observation", "Combined"]
    // [out] success - boolean - true if sampling was successful
    checkNumArgs(1,2);
    auto Nparticles = getInt();
    auto propType = getString();
    outputBool(obj->runParticleFilter(Nparticles,propType));
}

void PairInteractionSMC_Iface::obsLLH()
{
    // Get the estimated LLH for each observation
    // [out] obsLLH - size:[Ndata] - log-likelihood of each observation 
    checkNumArgs(1,0);
    outputDVec(obj->obsLLH());
}

void PairInteractionSMC_Iface::sampleParticle()
{
    // Sample a particle for each data and report LLH
    // [out] particle - cell-array size:[Ndata] each element is size:[5,Nobs] assigment to unknown variables 
    // [out] llh [optional] - size:[Ndata] log-likelihood of each particle 
    checkMaxNumArgs(2,0);
    PairInteractionSMC::MatListT p;
    PairInteractionSMC::VecT llh;
    obj->sampleParticle(p,llh);
    outputMatCellArray<double>(p);
    outputDVec(llh);
}

void PairInteractionSMC_Iface::sampleParticleLLH()
{
    // Report SampleLLH of sampled particle for each data
    // [out] llh - size:[Ndata] log-likelihood of each sampled particle 
    checkNumArgs(1,0);
    outputDVec(obj->sampleParticleLLH());
}

void PairInteractionSMC_Iface::getAllParticles()
{
    // Return all particles and the weights
    // [in] obsIdx - int [optional] 
    // [out] particles - cell-array size:[Ndata] each element is [5,Nobs,Nparticles] 
    // [out] weights - size:[Nparticles, Ndata]
    // [out] llh - size:[Nparticles, Ndata]
    checkMaxNumArgs(3,1);
    int obsIdx = -1;
    if(nrhs==1) obsIdx = getInt();
    PairInteractionSMC::CubeListT ps;
    PairInteractionSMC::MatT weight,llh;
    obj->getAllParticles(obsIdx,ps,weight,llh);
    outputCubeCellArray<double>(ps);
    outputDMat(weight);
    outputDMat(llh);
}

void PairInteractionSMC_Iface::computeLLH()
{
    // Sample report LLH of a sampled particle for each data
    // [in] ObsData [optional] - cell-array of obersevation matricies [default:internal data]
    // [in] particle - cell-array of particle matricies
    // [out] llh - size:[Ndata] llh of particle for each data
    checkMaxNumArgs(1,2);
    if(nrhs==2) {
        auto obs = getMatVector<double>();
        auto p = getMatVector<double>();
        outputDVec(obj->computeLLH(obs,p));
    } else {
        auto p = getMatVector<double>();
        outputDVec(obj->computeLLH(p));
    }
}

void PairInteractionSMC_Iface::computeLLH_debug()
{
    // Sample report all LLH of a sampled particle for each data
    // [in] ObsData [optional] - single observation matrix 
    // [in] particle - single particle matrix
    // [out] llhAll - matrix: [9,Nobs] llh of each variable for each observation
    checkNumArgs(1,2);
    auto obs = getMat<double>();
    auto p = getMat<double>();
    outputMat<double>(obj->computeLLH_debug(obs,p));
}

void PairInteractionSMC_Iface::simulate()
{
    // Use the DBN definition and the prior to simulate directly from the network
    // [in] Nparticles [uint64] - number of particles to simulate
    // [in] Nsteps [uint64]- number of steps to simulate
    // [in] tmax [double] - maximum simulation time.  Timesetep is interpolated based in this value
    // [in] simType [string] - 'simType': {"Ancestral","Blur","EC"} 
    // [in] simParams - struct of names doubles
    //                  'wmean':  [optional] mean value of error in observation (w) [Default: rho] 
    //                  'wgamma': [optional] relaxation rate of error in observation (w)  [Default: 10*dt]
    //                  'wsigma': [optional] standard deviation of error in observation [Default: .15*wmean]
    //                  'oversample': Required for "Blur".  Integer >1 number of additional samples to average for each measurment
    //                  'rho_unbind': Optional for "EC".  Unbinding radius.
    // [out] obs - cube size: [9,Nsteps,Nparticles] rows are [t, obs_xa, obs_xb, obs_ya, obs_xb, w_xa, w_xb, w_ya, w_yb]
    // [out] particles - cube size: [5,Nsteps,Nparticles] rows are [z, true_ax true_ay true_bx true_by]
    // [out] llhAll - matrix: [9,Nobs,Nparticles] llh of each hidden variable for each observation [z, yax yay ybx yby xax xay xbx xby]
    // [out] simParams - returned simulation parameters as a structure array of named doubles
    checkNumArgs(4,5);
    auto Nparticles = getScalar<IdxT>();
    auto Nsteps = getScalar<IdxT>();
    auto tmax = getDouble();
    auto simType = getString();
    auto params = getDoubleStruct();
    arma::cube obs, particles, llhAll;
    obj->simulate(Nparticles,Nsteps,tmax,simType,params,obs,particles,llhAll);
    outputStack<double>(obs);
    outputStack<double>(particles);
    outputStack<double>(llhAll);
    outputStatsToStruct(params);
}



void PairInteractionSMC_Iface::llhG0()
{
    // Compute llh(g0(y0, z0)) in parallel for N particles
    // [in] y0 - double size:[4,N] true positions rows: [y0_ax, y0_ay, y0_bx, y0_by]
    // [in] z0 - double length:[N] state [0=free, 1=bound]
    // [out] llh - double length:[N] log-likelihood for each particle (y0,z0)
    checkNumArgs(1,2);
    auto y0 = getMat<double>();
    auto z0 = getVec<double>();
    auto llh = makeVec<double>(z0.n_elem);
    obj->llhG0_parallel(y0,z0,llh);
}

void PairInteractionSMC_Iface::llhG()
{
    // Compute llh(g(yt, zt | y0, z0)) in parallel for N particles
    // [in] yt - double size:[4,N] true positions rows: [yt_ax, yt_ay, yt_bx, yt_by]
    // [in] zt - double length:[N] state [0=free, 1=bound]
    // [in] y0 - double size:[4,N] true positions rows: [y0_ax, y0_ay, y0_bx, y0_by]
    // [in] z0 - double length:[N] state [0=free, 1=bound]
    // [in] t - elapsed time
    // [out] llh - double length:[N] log-likelihood for each particle (yt,zt)
    checkNumArgs(1,2);
    auto yt = getMat<double>();
    auto zt = getVec<double>();
    auto y0 = getMat<double>();
    auto z0 = getVec<double>();
    auto t = getDouble();
    auto llh = makeVec<double>(z0.n_elem);
    obj->llhG_parallel(yt,zt,y0,z0,t,llh);
}

void PairInteractionSMC_Iface::llhH()
{
    // Compute llh(h(xt | yt, wt)) in parallel for N particles
    // [in] xt - double size:[4,N] true positions rows: [xt_ax, xt_ay, xt_bx, xt_by]
    // [in] wt - double size:[4,N] true positions rows: [wt_ax, wt_ay, wt_bx, wt_by]
    // [in] yt - double size:[4,N] true positions rows: [yt_ax, yt_ay, yt_bx, yt_by]
    // [out] llh - double length:[N] log-likelihood for each observation xt
    checkNumArgs(1,2);
    auto xt = getMat<double>();
    auto wt = getMat<double>();
    auto yt = getMat<double>();
    auto llh = makeVec<double>(wt.n_rows);
    obj->llhH_parallel(xt,wt,yt,llh);
}




