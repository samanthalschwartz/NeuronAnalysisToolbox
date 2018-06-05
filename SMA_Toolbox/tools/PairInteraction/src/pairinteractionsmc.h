/** @file pairinteractionsmc.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-2016
 * @brief Sequential Monte Carlo sampling of pair interaction model
 */

#ifndef _PAIRINTERACTIONSMC_H
#define _PAIRINTERACTIONSMC_H

#include <armadillo>
#include <stdexcept>
#include <map>
#include <string>
#include <memory>

#include "rdcapture.h"
#include "bessel.h"
#include "rngmanager.h"
#include "computations.h"

namespace pair_int {

class PositionPrior
{
public:
    enum class DistType {GAUSSIAN=1, RECTANGULAR=2, CIRCULAR=3};
    using VecT = arma::vec;
    using MatT = arma::mat;
    using ParamsT = std::map<std::string,double>;
    PositionPrior() {}
    virtual ~PositionPrior() = default;
    virtual double p_r0(double r0) const =0;
    virtual inline double llh_r0(double r0) const {return log(p_r0(r0));}
    virtual double llh_pos(const VecT &r0) const =0;
    virtual VecT sample_pos() const = 0;
    virtual ParamsT getParams() const = 0;
};


class PairInteractionSMC{
public:
    enum class SimulationType {ANCESTRAL=1, BLUR=2, ERBANCHAPMAN=3};
    enum class ProposalType {TRANSITION=1, OBSERVATION=2, COMBINED=3};
    using IdxT = arma::uword;
    using MatT = arma::mat;
    using VecT = arma::vec;
    using CubeT = arma::cube;
    using IVecT = arma::uvec; //Index vector
    using IMatT = arma::umat; //Index matrix
    using VecListT = std::vector<arma::vec>;
    using MatListT = std::vector<arma::mat>;
    using IMatListT = std::vector<arma::umat>; //
    using CubeListT = std::vector<arma::cube>; //
    using ParamsT = std::map<std::string,double>;
    
    double DA=-1;   //Diffusion constant for free A. [um^2/s]
    double DB=-1;   //Diffusion constant for free B. [um^2/s]
    double DAB=-1;  //Diffusion constant for bound AB. [um^2/s]
    double Dphi=-1; //Angular diffusion constant for bound AB. [rad^2/s]
    double rho=-1;  //Bound distance mean [um]
    double sigma_rho=-1; // Bound distance variance [um^2]
    double gamma=-1;     // Bound distance OU relaxation rate [1/s]
    double rho_bind=-1;  // RDCapture binding distance [um]
    double lambda_bind=-1; //RDCapture binding rate [1/s]
    double koff=-1; //Unbinding rate [1/s]
    
    //Uniform = uninformative prior
    //Gaussian = starting positions are expected to be Gaussian distributed
    
    
    PairInteractionSMC(const ParamsT &params);
    PairInteractionSMC(const ParamsT &params, const MatListT &data);
    
    /* Set and add data */
    IdxT addData(const MatListT &data);
    IdxT addData(const CubeT &data);
    IdxT addData(const MatT &data);
    void clearData();
    MatListT getData();
    IdxT getNData() const { return pair_data.size();}
    
    /* Set and get parameters structure */
    void setParams(const ParamsT &params);
    ParamsT getParams() const;
    
    /* Set prior */
    void setPriorState(const VecT &state_prior);
    void setPriorPositionGaussian(const VecT &center, const MatT &cov);
    void setPriorPositionRectangular(const VecT &rect);
    void setPriorPositionCircular(const VecT &center, double radius);
    ParamsT getPriorParams() const;
    
    /* Run a new particle filter with current params and prior */
    bool runParticleFilter(int Nparticles, const std::string &propType);
    
    /* Data retrieval.  These methods are only valid if the particle filter is valid */
    VecT obsLLH() const;
    void sampleParticle(MatListT &p, VecT &llh);
    VecT sampleParticleLLH();
    void getAllParticles(CubeListT &p, MatT &weights, MatT &llh) const; //For debugging
    void getAllParticles(int obsIdx, CubeListT &p, MatT &weights, MatT &llh) const; //For debugging, get up to any point in time
    
    /* Auxiliary routines */
    /* compute the LLH of one or more particles against given data */
    double computeLLH(const MatT &obs_data, const MatT &particle) const;
    VecT computeLLH(const MatListT &obs_data_arr, const MatListT &particle_arr) const;
    
    /* compute LLH of particles against internal data*/
    double computeLLH(IdxT data_idx, const MatT &particle) const;
    VecT computeLLH(const MatListT &particle_arr) const;
    VecT computeLLH(const CubeT &particle_arr) const;
    
    /* Debugging output */
    MatT computeLLH_debug(const MatT &obs_data, const MatT &particle) const;
    
    void simulate(IdxT Nparticles, IdxT Nsteps, double tmax, const std::string &simType, ParamsT &params, CubeT &obs, CubeT &particles, CubeT &llh) const;
    
    /* Parallel Versions for calls from Matlab*/
    //     void sample_Q(const VecT &x, const VecT &w, const MatT &y0, const VecT &z0, VecT &z, MatT &y, VecT &llh) const;
    void llhG0_parallel(const MatT &y0, const VecT &z0, VecT &llh) const;
    void llhG_parallel(const MatT &yt, const VecT &zt, const MatT &y0, const VecT &z0, double t, VecT &llh) const; 
    void llhH_parallel(const MatT &x, const MatT &w, const MatT &y, VecT &llh) const;
    
protected:
    static const constexpr double MIN_ESS=5;
    static const constexpr double MIN_ESS_RATIO=1/20.;
    
    static const constexpr double LOG2PI = log(2*arma::datum::pi);
    MatListT pair_data;  //Field size:[Ndata]: element size: [9,Nstates]: The observed data.
                     //rows are [t,ax,ay,bx,by,sax,say,sbx,sby]. 
    IVecT pair_Nobs; //size:[Ndata].  Number of data for each obs.
    
    VecT state_prior{.5,.5};
    VecT log_state_prior;
    double sigma_rho_sq; //square of sigma_rho
    std::shared_ptr<PositionPrior> pos_prior;
    
//     VecT captureCacheDist();
//     VecT captureCache
    
    /* Particle filter result variable */
    bool valid_particles=false; //.true if the filter data is from a valid run.
    IdxT Nparticles=0; // scalar - number of particles
    CubeListT pair_state; //Field size:[Ndata]: element size: [5,Nstates]:  Hidden states (shared amongst particles)
    IMatListT pair_ancestor; //Field size:[Ndata]: element size: [Nparticles,Nobs-1]:
                                //rows are particles, columns correspond to times.  col(i-1) is the time i-1 particle (row) index
                                //in pair_state_idx vector that corresponds to the previous particle value when evaluation is happending at time i
                                //
    MatT pair_weights;  //Size:[Nparticles,Ndata] The weight of each particle at final time.
    MatT pair_llh;   //Size:[Nparticles,Ndata].  The log-likelihood of each particle at final time.
    VecT pair_obs_llh;  //Size:[Ndata] The log-likelihood of the observation sequence as estimated by the particles

    bool valid_params;
    RDCapture rd; //RDCapture object for this class;
    
    void initDefaultPrior(); // Helper for constructor to set some sensible prior default
    void clearParticles();
    void initializeParticles(IdxT N);
    
    double computeLLH_particle(IdxT data_idx) const;
    MatT computeLLH_particle_debug(IdxT data_idx) const;
    
    bool runParticleFilter(ProposalType propType, const MatT &Obs, IdxT Nparticles, CubeT &state, IMatT &anscestor,
                           VecT &ws, VecT& llh, double &obsLLH) const;
    
    
    MatT reconstructParticle(IdxT data_idx, IdxT part_idx) const;
    MatT reconstructParticle(IdxT data_idx, IdxT part_idx, int obs_idx) const;
    IdxT drawParticleSample(IdxT data_idx);
    
    void simulateAncestral(IdxT Nparticles, IdxT Nsteps, double tmax, ParamsT &params, CubeT &obs, CubeT &particles, CubeT &llh) const;
    void simulateBlur(IdxT Nparticles, IdxT Nsteps, double tmax, ParamsT &params, CubeT &obs, CubeT &particles, CubeT &llh) const;
    void simulateEC(IdxT Nparticles, IdxT Nsteps, double tmax, ParamsT &params, CubeT &obs, CubeT &particles, CubeT &llh) const;
    
    
    void sampleQ0(const VecT &x0, const VecT &w0, MatT &y0, VecT &z0, VecT &llh) const;
    //Select either sampleQ or sampleG based on ProposalType
    void sampleProposal(ProposalType propType, const VecT &xt, const VecT &wt, const MatT &y0, const VecT &z0, 
                        double t, MatT &yt, VecT &zt, VecT &llh) const;
    ProposalType getProposalType(const std::string &propTypeStr) const;

    /* Vectorized */
    void sampleQcomb(const VecT &xt, const VecT &wt, const MatT &y0, const VecT &z0, double t, MatT &yt, VecT &zt, VecT &llh) const; 
    void sampleQobs(const VecT &xt, const VecT &wt, const MatT &y0, const VecT &z0, double t, MatT &yt, VecT &zt, VecT &llh) const; 
    void sampleG(const MatT &y0, const VecT &z0, double t, MatT &yt, VecT &zt, VecT &llh) const; 
    void sampleH_inv(const VecT &xt, const VecT &wt, MatT &yt, VecT &llh) const;
    VecT llhG_vec(const MatT &yt, const VecT &zt, const MatT &y0, const VecT &z0, double t) const; 
    VecT llhH_vec(const VecT &x, const VecT &w, const MatT &y) const;
    
                
    void sampleG0(VecT &y0, double& z0, VecT &llhAll) const;
    void sampleG(const VecT &y0, double z0, double t, VecT &yt, double &zt, VecT &llhAll) const; 
    void sampleH_all(const VecT &w, const VecT &y, VecT &x, VecT &llhAll) const;
    
    double llhG0(const VecT &y0, double z0) const;
    double llhG(const VecT &yt, double zt, const VecT &y0, double z0, double dt) const; 
    double llhH(const VecT &x, const VecT &w, const VecT &y) const;
    
    /* llhAll versions return llh of each element for debugging purposes */
    VecT llhG0_all(const VecT &y0, double z0) const;
    VecT llhG_all(const VecT &yt, double zt, const VecT &y0, double z0, double dt) const; 
    VecT llhH_all(const VecT &x, const VecT &w, const VecT &y) const;
    
    void sample_Q0(const VecT &x0, const VecT &w0, IdxT N, VecT &z0, MatT &y0, VecT &llh) const;

    double updateParamsDefault(ParamsT &params, const std::string &key, double default_val) const;
};


class GaussianPositionPrior : public PositionPrior 
{
    VecT mean;
    MatT cov;
    MatT chol_cov; //cholesky decomposition (square root)
    double A; // (1+q^2)/(q*omega)
    double B; // -(1+q^2)^2/(4*q^2*omega)
    double C; // (1-q^4)/(4*q^2*omega);
public:
    GaussianPositionPrior(const VecT &_mean, const MatT &_cov);
    double p_r0(double r0) const;
    double llh_r0(double r0) const;
    double p_pos(const VecT &a0) const;
    double llh_pos(const VecT &a0) const;
    VecT sample_pos() const;
    ParamsT getParams() const;
};

class RectangularPositionPrior : public PositionPrior 
{
    VecT rect;
    double width, height; //Do not rely on these matching x and y.  we swap them so width<=height always
    double area;
    double llh;
    double width_sq, height_sq;
public:
    RectangularPositionPrior(const arma::vec &_rect);//format is [left, bot, right, top] [x_min, y_min, x_max, y_max]
    double p_r0(double r0) const;
    inline double llh_pos(const VecT &a0) const {return llh;}
    VecT sample_pos() const;
    ParamsT getParams() const;
};

class CircularPositionPrior :public PositionPrior 
{
    VecT center;
    double radius;
    double area;
    double llh;
public:
    CircularPositionPrior(const arma::vec &center, double radius);
    double p_r0(double r0) const;
    inline double llh_pos(const VecT &a0) const {return llh;}
    VecT sample_pos() const;
    ParamsT getParams() const;
};


inline 
double GaussianPositionPrior::p_r0(double r) const
{
    double r_sq = r*r;
    return r*A*exp(B*r_sq)*Bessel::BE::I0(C*r_sq);
}

inline 
double GaussianPositionPrior::llh_r0(double r) const
{
    double r_sq = r*r;
    return log(r*A) + B*r_sq + Bessel::BE::logI0(C*r_sq);
}

inline 
double GaussianPositionPrior::llh_pos(const VecT &a0) const
{
    return pair_int::multiGaussianLLH((a0-mean).eval(),cov);
}

inline 
GaussianPositionPrior::VecT GaussianPositionPrior::sample_pos() const
{
    return chol_cov*RNG.randn(2);
}

inline 
RectangularPositionPrior::VecT  RectangularPositionPrior::sample_pos() const
{
    VecT tmp{rect(0) + (rect(2)-rect(0))*RNG.randu(), 
             rect(1) + (rect(3)-rect(1))*RNG.randu()};
    return tmp;
}

inline 
CircularPositionPrior::VecT CircularPositionPrior::sample_pos() const
{
    double r = sqrt(radius*radius*RNG.randu());
    double theta = 2*arma::datum::pi*RNG.randu();
    VecT tmp{center(0)+r*cos(theta), center(1)+r*sin(theta)};
    return tmp;
}


inline
double PairInteractionSMC::llhG0(const VecT &y0, double z0) const
{
    return arma::sum(llhG0_all(y0,z0));
}

inline
double PairInteractionSMC::llhG(const VecT &yt, double zt, const VecT &y0, double z0, double dt) const
{
    return arma::sum(llhG_all(yt, zt, y0, z0, dt));
}


inline 
double PairInteractionSMC::llhH(const VecT &x, const VecT &w, const VecT &y) const
{
    return pair_int::gaussianLLH((x-y).eval(),arma::square(w).eval());
}

inline
PairInteractionSMC::VecT PairInteractionSMC::llhH_all(const VecT &x, const VecT &w, const VecT &y) const
{
    return pair_int::gaussianLLH_all((x-y).eval(),arma::square(w).eval());
}

inline
void PairInteractionSMC::sampleH_all(const VecT &w, const VecT &y, VecT &x, VecT &llhAll) const
{
    x = y + w%RNG.randn(4);
    llhAll = pair_int::gaussianLLH_all((x-y).eval(),arma::square(w).eval());
}

inline
PairInteractionSMC::IdxT PairInteractionSMC::drawParticleSample(IdxT data_idx)
{
    return RNG.resample_dist( pair_weights.col(data_idx));  //Not thread safe
}


inline
void  PairInteractionSMC::sampleProposal(ProposalType propType, 
                const VecT &xt, const VecT &wt, const MatT &y0, const VecT &z0, 
                double t, MatT &yt, VecT &zt, VecT &llh) const
{
    switch(propType) {
        case ProposalType::TRANSITION:
            sampleG(y0, z0, t, yt, zt, llh);
            break;
        case ProposalType::COMBINED:
            sampleQcomb(xt,wt,y0, z0, t, yt, zt, llh);
            break;
        case ProposalType::OBSERVATION:
            sampleQobs(xt,wt,y0, z0, t, yt, zt, llh);
            break;
        default:
//             sampleH(xt,wt,y0, z0, t, yt, zt, llh);
            throw std::runtime_error("Bad proposal type");
    }
}
inline 
void PairInteractionSMC::getAllParticles(CubeListT &L, MatT &weights, MatT &llh) const
{
    getAllParticles(-1, L, weights, llh);
}
    
inline
PairInteractionSMC::MatT PairInteractionSMC::reconstructParticle(IdxT data_idx, IdxT pidx) const
{
    return reconstructParticle(data_idx,pidx,-1);
}

} /* namespace pair_int */
#endif /* _PAIRINTERACTIONSMC_H */
