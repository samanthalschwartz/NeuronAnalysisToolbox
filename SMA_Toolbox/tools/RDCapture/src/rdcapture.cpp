/** @file rdcapture.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief Simulate reaction diffusion to capture on a disk
 *
 */
#include "rdcapture.h"
#include "bessel.h"

RDCapture::RDCapture(double _D, double _rho, double _lambda): LinvGS()
{
    setParams(_D,_rho,_lambda);
}

void RDCapture::setParams(double _D, double _rho, double _lambda)
{
    D=_D;
    rho=_rho;
    lambda=_lambda;
    lambdaInv=1/lambda;
}

RDCapture::VecT RDCapture::survivalProb(const VecT &r0, const VecT &t) const
{
    const int N = static_cast<int>(r0.n_elem);
    assert(r0.n_elem==t.n_elem);
    VecT Q(N);
    for(int n=0; n<N; n++) Q(n) = survivalProb(r0(n), t(n));
    return Q;
}

RDCapture::VecT RDCapture::survivalProb(const VecT &r0, double t) const
{
    const int N = static_cast<int>(r0.n_elem);
    VecT Q(N);
    for(int n=0; n<N; n++) Q(n) = survivalProb(r0(n), t);
    return Q;
}

RDCapture::VecT RDCapture::survivalProb(double r0, const VecT &t) const
{
    const int N = static_cast<int>(t.n_elem);
    VecT Q(N);
    for(int n=0; n<N; n++)  Q(n) = survivalProb(r0, t(n));
    return Q;
}

void RDCapture::survivalProb_parallel(VecT &Q, const VecT &r0, const VecT &t) const
{
    const int N = static_cast<int>(Q.n_elem);
    assert(Q.n_elem==r0.n_elem);
    assert(Q.n_elem==t.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) Q(n) = survivalProb(r0(n),t(n));
}

void RDCapture::survivalProb_parallel(VecT &Q, const VecT &r0, double t) const
{
    const int N = static_cast<int>(Q.n_elem);
    assert(Q.n_elem==r0.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) Q(n) = survivalProb(r0(n),t);
}

void RDCapture::survivalProb_parallel(VecT &Q, double r0, const VecT &t) const
{
    const int N = static_cast<int>(Q.n_elem);
    assert(Q.n_elem==t.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) Q(n) = survivalProb(r0, t(n));
}

double RDCapture::mu(double r0, double t) const
{
    return LinvGS.invert([this, r0] (double s) {return muLT(r0,s);},t);
}

double RDCapture::nu(double t) const
{
    return LinvGS.invert([this] (double s) {return nuLT(s);},t);
}

void RDCapture::mu_parallel(VecT &muVal, const VecT &r0, const VecT &t) const
{
    const int N = static_cast<int>(muVal.n_elem);
    assert(muVal.n_elem==r0.n_elem);
    assert(muVal.n_elem==t.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) muVal(n) = LinvGS.invert([&] (double s) {return muLT(r0(n),s);},t(n));
}

void RDCapture::mu_parallel(VecT &muVal, double r0, const VecT &t) const
{
    const int N = static_cast<int>(muVal.n_elem);
    assert(muVal.n_elem==t.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) muVal(n) = LinvGS.invert([&] (double s) {return muLT(r0,s);},t(n));
}

void RDCapture::mu_parallel(VecT &muVal, const VecT &r0, double t) const
{
    const int N = static_cast<int>(muVal.n_elem);
    assert(muVal.n_elem==r0.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) muVal(n) = LinvGS.invert([&] (double s) {return muLT(r0(n),s);},t);
}


void RDCapture::nu_parallel(VecT &nuVal, const VecT &t) const
{
    const int N = static_cast<int>(nuVal.n_elem);
    assert(nuVal.n_elem==t.n_elem);
    #pragma omp parallel for
    for(int n=0; n<N; n++) nuVal(n) = LinvGS.invert([&] (double s) {return nuLT(s);},t(n));
}


double RDCapture::muLT(double r0, double s) const
{
	double omega=sqrt(s/D);
	if(r0<rho){
        return (1-rho*omega*Bessel::BE::I0K1(omega*r0,omega*rho))/s;
	} else {
        return (rho*omega*Bessel::BE::I1K0(omega*rho,omega*r0))/s;		
	}
}

double RDCapture::nuLT(double s) const
{
	double z=rho*sqrt(s/D);
    return (1-2*Bessel::BE::I1K1(z))/s;    
}

double RDCapture::survivalProbLT(double r0, double s) const
{
    return (1-muLT(r0,s)/(lambdaInv + nuLT(s)))/s;
}

void  RDCapture::simulate(const VecT &r0, const VecT &ts, IdxT N, double max_dt, MatT &survivalProb) const
{
    IdxT Nr0= r0.n_elem;
    IdxT Nts= ts.n_elem;
    #pragma omp parallel
    {
        RngT rng;
        VecT p_r0(Nts);
        #pragma omp for
        for(IdxT n=0; n<Nr0; n++)  {
            run_simulate(rng,r0(n),ts,N,max_dt,p_r0);
            survivalProb.col(n)=p_r0;
        }
    }
}

void  RDCapture::run_simulate(RngT &rng, double r0, const VecT &ts, IdxT N, double max_dt, VecT &survivalProb) const
{
    IdxT Nt = ts.n_elem;
    if(survivalProb.n_elem != Nt) throw std::invalid_argument(" survivalProb and ts must have the same length");
    arma::mat P(N,2);
    arma::vec capture_time(N,arma::fill::zeros);
    for(IdxT n=0; n<N; n++) {
        P(n,0)=r0;
        P(n,1)=0.;
    }
    double dt=max_dt;
    IdxT next_t_idx=0;
    double next_t = ts(next_t_idx);
    double t=0.;
    double sigma = sqrt(2*dt*D);
    std::normal_distribution<double> norm_d(0, sigma);
    std::uniform_real_distribution<double> uni_d;
    double capture_p = 1-exp(-lambda*dt);
    while(next_t_idx<Nt) {
        t+=dt;
        IdxT nFree=0;
        for(IdxT n=0;n<N;n++){
            if(capture_time(n)==0.){
                P(n,0)+= norm_d(rng);
                P(n,1)+= norm_d(rng);
                double r = sqrt(P(n,0)*P(n,0) + P(n,1)*P(n,1));
                if(r <= rho && uni_d(rng) <= capture_p) capture_time(n)=t;
                else nFree++;
            }
        }
        if(t>=next_t) {
            survivalProb(next_t_idx) = double(nFree)/N;
            next_t_idx++;
            if(next_t_idx<Nt) next_t = ts(next_t_idx);
        }
    }
}

#include "rdcapture.h"
