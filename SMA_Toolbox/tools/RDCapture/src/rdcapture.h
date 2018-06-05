/** @file rdcapture.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief Simulate reaction diffusion to capture on a disk
 *
 */

#ifndef _RDCAPTURE_H
#define _RDCAPTURE_H

#include <stdexcept>
#include <random>

#include <armadillo>

#include "laplacetform.h"

class RDCapture{
public:
    using IdxT = arma::uword;
    using VecT = arma::vec; // Type for vectors in this class
    using MatT = arma::mat; // Type for vectors in this class
    double D=0;  //Diffusion constant (um^2/s)
    double rho=0; //Capture radius (um)
    double lambda=0; //Capture rate (1/s)
    
    RDCapture() {}; 
    RDCapture(double _D, double _rho, double _lambda); 

    void setParams(double _D, double _rho, double _lambda);
    /* Top Level factorized calls built around survivalProbLT  */
    double survivalProb(double r0, double t)  const;
    /* Vectorized */
    VecT survivalProb(const VecT &r0, const VecT &t)  const;
    VecT survivalProb(double r0, const VecT &t)  const;
    VecT survivalProb(const VecT &r0, double t)  const;
    /* Parallel vectorized versions with openMP */
    void survivalProb_parallel(VecT &Q, const VecT &r0, const VecT &t)  const;
    void survivalProb_parallel(VecT &Q, double r0, const VecT &t)  const;
    void survivalProb_parallel(VecT &Q, const VecT &r0, double t)  const;
    
    /* Compute mu or nu from LT */
    void mu_parallel(VecT &muVals, const VecT &r0, const VecT &t) const;
    void mu_parallel(VecT &muVals, double r0, const VecT &t) const;
    void mu_parallel(VecT &muVals, const VecT &r0, double t) const;
    
    void nu_parallel(VecT &nuVals, const VecT &t) const;
    
    double mu(double r0, double t) const;
    double nu(double t) const;
    
    /* Core computational routines */
    double survivalProbLT(double r0, double s)  const;
    double muLT(double r0, double s)  const;
    double nuLT(double s) const;

    double survivalLogProb(double r0, double t)  const;
    double captureLogProb(double r0, double t)  const;
    
    void simulate(const VecT &r0, const VecT &ts, IdxT N, double max_dt, MatT &survivalProb) const;
private:
    using RngT = std::mt19937;
    static constexpr double NSIGMA_CUTTOFF = 6.; //Number of sigmas out at which point to avoid computing survival prob
    double lambdaInv;
    LaplaceInverseGS LinvGS; /* Gavier-Stehfest inverter */
    bool captureFeasible(double r0, double t) const;
    void run_simulate(RngT & rng, double r0, const VecT &ts, IdxT N, double max_dt, VecT &survivalProb) const;
};

/* Inlined methods */

inline
bool RDCapture::captureFeasible(double r0, double t) const
{
    return (r0<=rho) || (r0-rho <= NSIGMA_CUTTOFF*sqrt(2*D*t));
}

inline
double RDCapture::survivalProb(double r0, double t) const
{
//     r0 = std::max(fabs(r0),1e-9); //Round away from 0
//     if (!(D>0 && rho>0 && lambda>0)) throw std::logic_error("RDCapture invalid parameter values.");
    if(captureFeasible(r0,t)) {
        double captureP = LinvGS.invert([r0,this] (double s) {return survivalProbLT(r0,s);},t);
        /* TODO: FIGURE THIS OUT */
        captureP = std::min(std::max(captureP,0.),1.);
        return captureP;
    } else {
        return 1.;
    }
}   

inline
double RDCapture::survivalLogProb(double r0, double t) const
{
    return log(survivalProb(r0,t));
}

inline
double RDCapture::captureLogProb(double r0, double t) const
{
    double p = survivalProb(r0,t);
    if(p==1.) return 0.;
    else if(std::isnan(p)) return -INFINITY;
    else if(p<1e-3) return log(1-p);
    else return log1p(-p);
}

#endif /* RDCAPTURE_H */
