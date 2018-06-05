/** @file interaction_change_point.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 07-2015
 *  @brief The member definitions for InteractionChangePoint
 */
#include <cassert>
#include <iostream>
#include "interaction_change_point.h"

const InteractionChangePoint::FloatT InteractionChangePoint::log2pi = log(2*arma::Datum<FloatT>::pi);

InteractionChangePoint::InteractionChangePoint(FloatT _D_A, FloatT _D_B, FloatT _D_AB, FloatT _exposureT, FloatT _rho)
    : D_A(_D_A), D_B(_D_B), D_AB(_D_AB), exposureT(_exposureT), rho(_rho)
{
    assert(D_A>0 && D_B>0 && D_AB>0);
    assert(exposureT>0);
}


void InteractionChangePoint::freeSegmentLLH(const MatT &pair, MatT &LLH) const
{
    //Cols of LLH correspond to "start" locations. Rows correspond to "end" locations.
    int N = static_cast<int>(pair.n_rows);
    VecT T = pair.col(0); //observation times
    //Precompute the variance due to measurment and the variance due to diffusion (transition)
    MatT vD(N-1, 2), vM(N,4);
    vD.col(0) = computeDVariance(D_A, T); // variance T_a  (omega_a)
    vD.col(1) = computeDVariance(D_B, T); // variance T_b  (omega_b)
    vM.col(0) = computeMVariance(D_A, pair.col(9)); // variance M_a_x
    vM.col(1) = computeMVariance(D_A, pair.col(10)); // variance M_a_y
    vM.col(2) = computeMVariance(D_B, pair.col(13)); // variance M_b_x
    vM.col(3) = computeMVariance(D_B, pair.col(14)); // variance M_b_y

    LLH.zeros();
    #pragma omp parallel for schedule(static,1)
    for(int start=0; start<N-1; start++){
        LLH(arma::span(start+1,N-1), start) += free1DSegmentLLH(pair.col(1), vM.col(0), vD.col(0), start); // LLH_free_a_x
        LLH(arma::span(start+1,N-1), start) += free1DSegmentLLH(pair.col(2), vM.col(1), vD.col(0), start); // LLH_free_a_y
        LLH(arma::span(start+1,N-1), start) += free1DSegmentLLH(pair.col(5), vM.col(2), vD.col(1), start); // LLH_free_b_x
        LLH(arma::span(start+1,N-1), start) += free1DSegmentLLH(pair.col(6), vM.col(3), vD.col(1), start); // LLH_free_b_y
    }
}


void
InteractionChangePoint::boundSegmentLLH(const MatT &pair, MatT &LLH) const
{
    //Cols of LLH correspond to "start" locations. Rows correspond to "end" locations.
    int N = static_cast<int>(pair.n_rows);
    VecT T = pair.col(0); //observation times
    //Precompute the variance due to measurment and the variance due to diffusion (transition) Note: D_AB is the only diffusion constant used
    VecT vD(N-1);
    MatT vM(N,4);
    vD.col(0) = computeDVariance(D_AB, T); // variance T_a  (omega_a)
    vM.col(0) = computeMVariance(D_AB, pair.col(9)); // variance M_a_x
    vM.col(1) = computeMVariance(D_AB, pair.col(10)); // variance M_a_y
    vM.col(2) = computeMVariance(D_AB, pair.col(13)); // variance M_b_x
    vM.col(3) = computeMVariance(D_AB, pair.col(14)); // variance M_b_y
    LLH.zeros();
    #pragma omp parallel for schedule(static,1)
    for(int start=0; start<N-1; start++){
        LLH(arma::span(start+1,N-1), start) += bound1DSegmentLLH(pair.col(1), pair.col(5), vM.col(0), vM.col(2), vD, start); // LLH_bound_x
        LLH(arma::span(start+1,N-1), start) += bound1DSegmentLLH(pair.col(2), pair.col(6), vM.col(1), vM.col(3), vD, start); // LLH_bound_y
    }
}

InteractionChangePoint::FloatT
InteractionChangePoint::sequenceLLH(const MatT &pair, const Sequence &seq) const
{
    
    return full1DLLH(pair.col(0), pair.col(1), pair.col(5), pair.col(9)%pair.col(9), pair.col(13)%pair.col(13), seq)
         + full1DLLH(pair.col(0), pair.col(2), pair.col(6), pair.col(10)%pair.col(10), pair.col(14)%pair.col(14), seq);
}

InteractionChangePoint::FloatT
InteractionChangePoint::sequenceLLH_debug(const MatT &pair, const Sequence &seq, CubeT &LLH_components) const
{
    int N = seq.end-seq.begin+1; //maximum number of LLH LLH_dist_terms
    FloatT LLH_val=0;
    IndexT nRows = 5;
    MatT LLH_components_1D(nRows,2*N-1,arma::fill::zeros);
    IndexT Ncomponents=0;
    std::cout<<"N: "<<N<<" Ncomponents:"<<Ncomponents<<" size:["<<LLH_components_1D.n_rows<<","<<LLH_components_1D.n_cols<<"]"<<std::endl;
    LLH_val+=full1DLLH_debug(pair.col(0), pair.col(1), pair.col(5), pair.col(9)%pair.col(9), pair.col(13)%pair.col(13), seq, Ncomponents, LLH_components_1D);
    LLH_components.set_size(nRows,Ncomponents,2);
    std::cout<<"N: "<<N<<" Ncomponents:"<<Ncomponents<<" size:["<<LLH_components_1D.n_rows<<","<<LLH_components_1D.n_cols<<"]"<<std::endl;
    std::cout<<" cube size:["<<LLH_components.n_rows<<","<<LLH_components.n_cols<<","<<LLH_components.n_slices<<"]"<<std::endl;

    LLH_components.slice(0) = LLH_components_1D.cols(0,Ncomponents-1);

    LLH_components_1D.zeros();
    LLH_val+=full1DLLH_debug(pair.col(0), pair.col(2), pair.col(6), pair.col(10)%pair.col(10), pair.col(14)%pair.col(14), seq, Ncomponents, LLH_components_1D);
    std::cout<<"N: "<<N<<" Ncomponents:"<<Ncomponents<<" size:["<<LLH_components.n_rows<<","<<LLH_components.n_cols<<"]"<<std::endl;
    LLH_components.slice(1) = LLH_components_1D.cols(0,Ncomponents-1);
    return LLH_val;
}


InteractionChangePoint::VecT
InteractionChangePoint::compute1CPLLH(const MatT &pair, const Sequence &seq) const
{
    IndexT nCP = seq.end-seq.begin-1; // number change points to test seq.begin+1 .. seq.end-1
    VecT LLH(nCP,arma::fill::zeros);
    //This could be more efficient if we saved work from previous run
    VecT var_ax=pair.col(9)%pair.col(9);
    VecT var_ay=pair.col(10)%pair.col(10);
    VecT var_bx=pair.col(13)%pair.col(13);
    VecT var_by=pair.col(14)%pair.col(14);
    
    #pragma omp parallel for
    for(IndexT n=0; n<nCP; n++){
        IVecT cp = {seq.begin+1+n};
        auto can_seq = makeSequence(seq.begin, seq.end, seq.state0, cp);
        LLH(n)  = full1DLLH(pair.col(0), pair.col(1), pair.col(5), var_ax, var_bx, can_seq);
        LLH(n) += full1DLLH(pair.col(0), pair.col(2), pair.col(6), var_ay, var_by, can_seq);
    }
    return LLH;
}

InteractionChangePoint::MatT
InteractionChangePoint::compute2CPLLH(const MatT &pair, const Sequence &seq) const
{
    IndexT nCP = seq.end-seq.begin-1; // number change points to test seq.begin+1 .. seq.end-1
    MatT LLH(nCP,nCP,arma::fill::zeros);
    //This could be more efficient if we saved work from previous run
    VecT var_ax=pair.col(9)%pair.col(9);
    VecT var_ay=pair.col(10)%pair.col(10);
    VecT var_bx=pair.col(13)%pair.col(13);
    VecT var_by=pair.col(14)%pair.col(14);
    
    #pragma omp parallel for
    for(IndexT n=0; n<nCP-1; n++){
        for(IndexT m=n+1; m<nCP; m++){
            IVecT cp = {seq.begin+1+n, seq.begin+1+m};
            auto can_seq = makeSequence(seq.begin, seq.end, seq.state0, cp);
            LLH(m,n)  = full1DLLH(pair.col(0), pair.col(1), pair.col(5), var_ax, var_bx, can_seq);
            LLH(m,n) += full1DLLH(pair.col(0), pair.col(2), pair.col(6), var_ay, var_by, can_seq);
        }
    }
    return LLH;
}

/* Public Static Methods */

InteractionChangePoint::FloatT
InteractionChangePoint::segmentProductLLH(const MatT &freeLLH, const MatT &boundLLH,
                                          const Sequence &seq)
{
    //seq_cps - 0-based index of changepoints.  all(0<seq_cps<N-1), and they are in-order.
    //recall that freeLLH and boundLLH are lower-triangular column-major format, so syntax for liklihood of segment
    // i:j with i<j is freeLLH(j,i) (i.e., it is (end_idx,start_idx))
    // seq_start>=0;
    // seq_end<=N-1;
    bool state=seq.state0;
    int K = static_cast<IndexT>(seq.cps.n_elem); //number of change points
    FloatT llh=0;
    if(K==0) { //No change points --- 1 segment
        if(state) llh += boundLLH(seq.end, seq.begin);
        else      llh += freeLLH(seq.end, seq.begin);
    } else { //Has at least 2 segments
        //First segment
        if(state) llh += boundLLH(seq.cps(0), seq.begin);
        else      llh += freeLLH(seq.cps(0), seq.begin);
        state = !state;
        //Middle segments
        for(int k=0; k<K-1; k++) {
            if(state) llh += boundLLH(seq.cps(k+1), seq.cps(k));
            else      llh += freeLLH(seq.cps(k+1), seq.cps(k));
            state = !state;
        }
        //Final segment
        if(state) llh += boundLLH(seq.end, seq.cps(K-1));
        else      llh += freeLLH(seq.end, seq.cps(K-1));
    }
    return llh;
}


InteractionChangePoint::FloatT
InteractionChangePoint::segmentProductMLE(const MatT &freeLLH, const MatT &boundLLH, double beta, Sequence &seq)
{
    int N = static_cast<IndexT>(freeLLH.n_rows); //number of observations
    VecT maxLLH0(N-1,arma::fill::zeros), maxLLH1(N-1,arma::fill::zeros);
    IVecT mleCP0(N-1,arma::fill::zeros), mleCP1(N-1,arma::fill::zeros);
    computeProductMaxLikelihood(freeLLH, boundLLH, beta, maxLLH0, maxLLH1, mleCP0, mleCP1);
    seq.begin = 0;
    seq.end = N-1;
    seq.state0 = maxLLH1(0)>maxLLH0(0); //true if max likelihood happens when initial state is 1
    FloatT seq_llh = seq.state0 ? maxLLH1(0) : maxLLH0(0);
    IVecT tempCPs(N-1);
    int k=0;
    IndexT nextCP;
    bool state = seq.state0;
    if(state) nextCP = mleCP1(0);
    else      nextCP = mleCP0(0);
    while(nextCP!=N-1){
        tempCPs(k++) = nextCP;
        state = !state;
        if(state) nextCP = mleCP1(nextCP);
        else      nextCP = mleCP0(nextCP);
    }
    seq.cps.resize(k);
    seq.cps = tempCPs.head(k);
    return seq_llh;
}

void InteractionChangePoint::computeProductMaxLikelihood(const MatT &freeLLH, const MatT &boundLLH, double beta,
                                                  VecT &maxLLH0, VecT &maxLLH1, IVecT &mleCP0, IVecT &mleCP1)
{
    int N = static_cast<IndexT>(freeLLH.n_rows); //number of observations
    for(int n=N-2; n>=0; n--){
        //compute for init_state=0
        mleCP0(n) = N-1;
        maxLLH0(n) = freeLLH(N-1, n);
        for(int m=n+1; m<N-1; m++){
            FloatT temp = freeLLH(m, n) + maxLLH1(m) + beta;
            if(temp>maxLLH0(n)){
                maxLLH0(n) = temp;
                mleCP0(n) = m;
            }
        }
        //compute for init_state=1
        mleCP1(n) = N-1;
        maxLLH1(n) = boundLLH(N-1, n);
        for(int m=n+1; m<N-1; m++){
            FloatT temp = boundLLH(m, n) + maxLLH0(m) + beta;
            if(temp>maxLLH1(n)){
                maxLLH1(n) = temp;
                mleCP1(n) = m;
            }
        }
    }
}








/* Protected Methods  */

InteractionChangePoint::VecT
InteractionChangePoint::free1DSegmentLLH(const VecT &obs, const VecT &vM, const VecT &vD, int start) const
{
    int N = static_cast<FloatT>(obs.n_elem);  //Total number of observations
    VecT LLH(N-start-1, arma::fill::zeros);
    FloatT mu = obs(start); //mu_1
    FloatT eta = vD(start)+vM(start); //eta_1
    FloatT alpha = vM(start+1) + eta; //alpha_1
    LLH(0) = log2pi + log(alpha) + square(obs(start+1)-mu)/alpha;
    for(int n=start+1; n<N-1; n++){
        mu = (mu*vM(n) + obs(n)*eta)/alpha; //mu_n
        eta = vM(n)*eta/alpha + vD(n); //eta_n
        alpha = vM(n+1)+eta; //alpha_n
        LLH(n-start) = LLH(n-start-1) + log2pi + log(alpha) + square(obs(n+1)-mu)/alpha;
    }
    LLH *= -0.5; //Factored out -0.5 from everything above
    return LLH;
}

InteractionChangePoint::VecT
InteractionChangePoint::bound1DSegmentLLH(const VecT &obs_a, const VecT &obs_b, const VecT &vMa, const VecT &vMb, const VecT &vD, int start) const
{
    int N = static_cast<FloatT>(obs_a.n_elem);  //Total number of observations
    VecT LLHvec(N-start-1, arma::fill::zeros);
    FloatT rho_sq = rho*rho;
    FloatT beta = rho_sq + vMb(start); //beta(start)
    FloatT gamma = beta + vMa(start); //gamma(start)
    FloatT kappa = (obs_a(start)*beta+obs_b(start)*vMa(start))/gamma; //kappa(start)
    FloatT zeta = vMa(start)*beta/gamma; //zeta(start)
    FloatT LLH = log2pi+log(gamma)+square(obs_a(start)-obs_b(start))/gamma; //N(obs_a(start),obs_b(start),gamma(start))
    FloatT mu = kappa; //mu(start)
    FloatT eta = zeta + vD(start); //eta(start)
    //     std::cout<<"N:"<<N<<" beta:"<<beta<<" gamma:"<<gamma<<" kappa:"<<kappa<<" zeta:"<<zeta<<" mu:"<<mu<<" eta:"<<eta<<" LLH:"<<LLH<<std::endl;
    for(int n=start+1; n<N; n++){
        beta = rho_sq + vMb(n); //beta(n)
        gamma = beta+vMa(n); //gamma(n)
        kappa = (obs_a(n)*beta+obs_b(n)*vMa(n))/gamma; //kappa(n)
        zeta = vMa(n)*beta/gamma; //zeta(n)
        FloatT alpha = zeta + eta; //alpha(n) = zeta(n) + eta(n-1) [defined for n>start]
        if (alpha<0){
            std::cout<<"!!! Negative alpha. N="<<n<<"/"<<N<<" alpha:"<<alpha<<std::endl;
        }
        LLH += 2*log2pi + log(fabs(gamma*alpha)) + square(obs_a(n)-obs_b(n))/gamma + square(mu-kappa)/alpha; //N(mu(n-1),kappa(n),alpha(n))
        LLHvec(n-start-1) = LLH;
//         std::cout<<"N:"<<n<<"/"<<N<<" beta:"<<beta<<" gamma:"<<gamma<<" kappa:"<<kappa<<" zeta:"<<zeta<<" alpha:"<<alpha<<" mu:"<<mu<<" eta:"<<eta<<" LLH:"<<LLH<<std::endl;
        if(n<N-1) { //defined only for [n<N-1]
            mu = (mu*zeta+kappa*eta)/alpha; //mu(n)=(mu(n-1)*zeta(n)+kappa(n)*eta(n-1))/alpha(n)
            eta = zeta*eta/alpha+vD(n);// eta_n
        }
    }
    LLHvec *= -0.5; //Factored out -0.5 from everything above
    return LLHvec;
}


InteractionChangePoint::FloatT
InteractionChangePoint::full1DLLH(const VecT &T, const VecT &obs_a, const VecT &obs_b, const VecT &var_a, const VecT &var_b,
                                  const Sequence &seq) const
{
    FloatT LLH=0;
    FloatT epsilon_a, epsilon_b; //Variance due to measurement defined [0..N-1]
    FloatT omega_a, omega_b; //Variance due to diffusion jump defined [0..N-2]
    FloatT mu_a, mu_b; //Mean of carried Normal Distribution
    FloatT eta_a, eta_b; //Variance of carried Normal Distribution 
    FloatT alpha_a, alpha_b; //Varaince of factored Normal distribtions
    FloatT beta, b_pos; //mean and position for factored Normal distribution
    FloatT kappa, zeta, gamma;
    FloatT rho_sq = rho*rho; //variance of dimer distance constraint R
    VecT dT = arma::diff(T); 
    BVecT state_seq = seq.computeStateSeq();
    if(seq.state0){ //state0 = bound
        //Compute variances
        epsilon_a = boundVariance(var_a(seq.begin)-D_AB*exposureT/3.);
        epsilon_b = boundVariance(var_b(seq.begin)-D_AB*exposureT/3.);
        omega_a = boundVariance(2*D_AB*dT(seq.begin));
        
        beta = epsilon_b+rho_sq;
        b_pos = obs_b(seq.begin);
        gamma = epsilon_a + beta;
        kappa = (obs_a(seq.begin)*beta+b_pos*epsilon_a)/gamma;
        zeta = epsilon_a*beta/gamma;
        LLH += normalDistLLH(obs_a(seq.begin), b_pos, gamma);
        mu_a = kappa;
        eta_a = zeta + omega_a;
        mu_b = 0; //To prevent warning
        eta_b = 0; //To prevent warning
    }else{//state0 = free
        epsilon_a = boundVariance(var_a(seq.begin)-D_A*exposureT/3.);
        epsilon_b = boundVariance(var_b(seq.begin)-D_B*exposureT/3.);
        omega_a = boundVariance(2*D_A*dT(seq.begin));
        omega_b = boundVariance(2*D_B*dT(seq.begin));
        //mu_a and eta_a define L_Fa = N(a_i+1, mu_a, eta_a) - Carried to next level i+1
        mu_a = obs_a(seq.begin);
        eta_a = epsilon_a+omega_a;
        //mu_b and eta_b define L_Fb = N(b_i+1, mu_b, eta_b) - Carried to next level i+1
        mu_b = obs_b(seq.begin);
        eta_b = epsilon_b+omega_b;
    }
    for(int n=seq.begin+1; n<=seq.end; n++){
        bool state = state_seq[n-seq.begin];
        bool last_state = state_seq[n-seq.begin-1];
//         std::cout<<"n:"<<n<<" state:"<<state<<" laststate:"<<last_state<<std::endl;
        if(state){ //Current state is bound
            if(last_state!=state) { //F->B.  Carry in L_Fa(n-1) L_Fb(n-1)
                assert(n<seq.end); //cannot have a change point on last step
                epsilon_a = boundVariance(var_a(n)-D_A*exposureT/6.-D_AB*exposureT/6.);
                epsilon_b = boundVariance(var_b(n)-D_B*exposureT/6.-D_AB*exposureT/6.);
                //Now beta and b_pos work in the normal distirbution of L_Fb_n-1 which mentions b_n
                beta = epsilon_b*eta_b/(epsilon_b+eta_b)+rho_sq;
                b_pos = (mu_b*epsilon_b+obs_b(n)*eta_b)/(epsilon_b+eta_b);
                gamma = epsilon_a + beta;
                kappa = (obs_a(n)*beta+b_pos*epsilon_a)/gamma;
                zeta = epsilon_a*beta/gamma;
                alpha_b = epsilon_b+eta_b;
                alpha_a = zeta + eta_a;
                LLH += normalDistLLH(obs_b(n), mu_b, alpha_b); //drops out from integration over db_n
                LLH += normalDistLLH(obs_a(n), b_pos, gamma);
                LLH += normalDistLLH(mu_a, kappa, alpha_a);

                //Compute recusrive variables for next step
                omega_a = boundVariance(2*D_AB*dT(n));
                mu_a = (mu_a*zeta+kappa*eta_a)/alpha_a;
                eta_a = zeta*eta_a/alpha_a + omega_a;
            } else { //B->B
                epsilon_a = boundVariance(var_a(n)-D_AB*exposureT/3.);
                epsilon_b = boundVariance(var_b(n)-D_AB*exposureT/3.);
                beta = epsilon_b + rho_sq;
                b_pos = obs_b(n);
                gamma = epsilon_a + beta;
                kappa = (obs_a(n)*beta+b_pos*epsilon_a)/gamma;
                zeta = epsilon_a*beta/gamma;
                alpha_a = zeta + eta_a;
                LLH += normalDistLLH(obs_a(n), b_pos, gamma);
                LLH += normalDistLLH(mu_a, kappa, alpha_a);

                if(n==seq.end) break; //we are done recusing

                //Compute recusrive variables for next step
                omega_a = boundVariance(2*D_AB*dT(n));
                mu_a = (mu_a*zeta+kappa*eta_a)/alpha_a;
                eta_a = zeta*eta_a/alpha_a + omega_a;
//                 std::cout<<" omega_a:"<<omega_a<<" mu_a(n):"<<mu_a<<" eta_a(n):"<<eta_a<<std::endl;
            }
        } else { //Current state is free
            if(last_state!=state) { //B->F  Carry in L_B(n-1) using mu_a, eta_a
                assert(n<seq.end); //cannot have a change point on last step
                epsilon_a = boundVariance(var_a(n)-D_A*exposureT/6.-D_AB*exposureT/6.);
                epsilon_b = boundVariance(var_b(n)-D_B*exposureT/6.-D_AB*exposureT/6.);
                //Add the extracted normal distributions to the LLH sum
                //A carries in L_B(n-1) = N(a_n, mu_a_{n-1}, eta_a_{n-1});
                //B does not have any relation to previous integral so "start fresh"
                alpha_a = epsilon_a + eta_a;
                LLH += normalDistLLH(obs_a(n), mu_a, alpha_a);

                //Compute recusrsive variables for next step.
                mu_a = (mu_a*epsilon_a + obs_a(n)*eta_a)/alpha_a;
                omega_a = boundVariance(2*D_A*dT(n));
                eta_a = eta_a*epsilon_a/alpha_a + omega_a;

                mu_b = obs_b(n);
                omega_b = boundVariance(2*D_B*dT(n));
                eta_b = epsilon_b+omega_b;
            } else { //F->F
                epsilon_a = boundVariance(var_a(n)-D_A*exposureT/3.);
                epsilon_b = boundVariance(var_b(n)-D_B*exposureT/3.);

                //Add the extracted normal distributions to the LLH sum
                alpha_a = epsilon_a + eta_a;
                LLH += normalDistLLH(obs_a(n), mu_a, alpha_a);

                alpha_b = epsilon_b + eta_b;
                LLH += normalDistLLH(obs_b(n), mu_b, alpha_b);
                if(n==seq.end) break; //we are done recusing

                //Compute recusrive variables for next step
                mu_a = (mu_a*epsilon_a + obs_a(n)*eta_a)/alpha_a;
                omega_a = boundVariance(2*D_A*dT(n));
                eta_a = eta_a*epsilon_a/alpha_a + omega_a;

                mu_b = (mu_b*epsilon_b + obs_b(n)*eta_b)/alpha_b;
                omega_b = boundVariance(2*D_B*dT(n));
                eta_b = eta_b*epsilon_b/alpha_b + omega_b;
            }
        }
    }
    LLH *= -0.5;
    return LLH;
}

InteractionChangePoint::FloatT
InteractionChangePoint::full1DLLH_debug(const VecT &T, const VecT &obs_a, const VecT &obs_b, const VecT &var_a, const VecT &var_b,
                                        const Sequence &seq, IndexT &Ncomponents, MatT &LLH_components) const
{
    FloatT LLH=0;
    FloatT epsilon_a, epsilon_b; //Variance due to measurement defined [0..N-1]
    FloatT omega_a, omega_b; //Variance due to diffusion jump defined [0..N-2]
    FloatT mu_a, mu_b; //Mean of carried Normal Distribution
    FloatT eta_a, eta_b; //Variance of carried Normal Distribution
    FloatT alpha_a, alpha_b; //Varaince of factored Normal distribtions
    FloatT beta, b_pos; //mean and position for factored Normal distribution
    FloatT kappa, zeta, gamma;
    FloatT rho_sq = rho*rho; //variance of dimer distance constraint R
    VecT dT = arma::diff(T);
    BVecT state_seq = seq.computeStateSeq();

    //Debugging variables
    Ncomponents = 0;

    if(seq.state0){ //state0 = bound
        //Compute variances
        epsilon_a = boundVariance(var_a(seq.begin)-D_AB*exposureT/3.);
        epsilon_b = boundVariance(var_b(seq.begin)-D_AB*exposureT/3.);
        omega_a = boundVariance(2*D_AB*dT(seq.begin));

        beta = epsilon_b+rho_sq;
        b_pos = obs_b(seq.begin);
        gamma = epsilon_a + beta;
        kappa = (obs_a(seq.begin)*beta+b_pos*epsilon_a)/gamma;
        zeta = epsilon_a*beta/gamma;
        LLH_components(0,Ncomponents) = seq.begin;
        LLH_components(1,Ncomponents) = 0;
        LLH += normalDistLLH_debug(obs_a(seq.begin), b_pos, gamma, Ncomponents, LLH_components);
        mu_a = kappa;
        eta_a = zeta + omega_a;
        mu_b = 0; //To prevent warning
        eta_b = 0; //To prevent warning
    }else{//state0 = free
        epsilon_a = boundVariance(var_a(seq.begin)-D_A*exposureT/3.);
        epsilon_b = boundVariance(var_b(seq.begin)-D_B*exposureT/3.);
        omega_a = boundVariance(2*D_A*dT(seq.begin));
        omega_b = boundVariance(2*D_B*dT(seq.begin));
        //mu_a and eta_a define L_Fa = N(a_i+1, mu_a, eta_a) - Carried to next level i+1
        mu_a = obs_a(seq.begin);
        eta_a = epsilon_a+omega_a;
        //mu_b and eta_b define L_Fb = N(b_i+1, mu_b, eta_b) - Carried to next level i+1
        mu_b = obs_b(seq.begin);
        eta_b = epsilon_b+omega_b;
    }
    for(int n=seq.begin+1; n<=seq.end; n++){
        bool state = state_seq[n-seq.begin];
        bool last_state = state_seq[n-seq.begin-1];
        //         std::cout<<"n:"<<n<<" state:"<<state<<" laststate:"<<last_state<<std::endl;
        if(state){ //Current state is bound
            if(last_state!=state) { //F->B.  Carry in L_Fa(n-1) L_Fb(n-1)
                assert(n<seq.end); //cannot have a change point on last step
                epsilon_a = boundVariance(var_a(n)-D_A*exposureT/6.-D_AB*exposureT/6.);
                epsilon_b = boundVariance(var_b(n)-D_B*exposureT/6.-D_AB*exposureT/6.);
                //Now beta and b_pos work in the normal distirbution of L_Fb_n-1 which mentions b_n
                beta = epsilon_b*eta_b/(epsilon_b+eta_b)+rho_sq;
                b_pos = (mu_b*epsilon_b+obs_b(n)*eta_b)/(epsilon_b+eta_b);
                gamma = epsilon_a + beta;
                kappa = (obs_a(n)*beta+b_pos*epsilon_a)/gamma;
                zeta = epsilon_a*beta/gamma;
                alpha_b = epsilon_b+eta_b;
                alpha_a = zeta + eta_a;
                //                 std::cout<<"B->F: "<<"eps:"<<epsilon_a<<","<<epsilon_b<<" beta:"<<beta<<" b_pos:"<<b_pos
                //                          <<" gamma:"<<gamma<<" kappa:"<<kappa<<" zeta:"<<zeta<<" alpha_a:"<<alpha_a<<" alpha_b:"<<alpha_b
                //                          <<" mu_a(n-1):"<<mu_a<<" mu_b(n-1):"<<mu_b<<" eta_a(n-1):"<<eta_a<<" eta_b(n-1):"<<eta_b;
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 2;
                LLH += normalDistLLH_debug(obs_b(n), mu_b, alpha_b, Ncomponents, LLH_components);
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 2;
                LLH += normalDistLLH_debug(obs_a(n), b_pos, gamma, Ncomponents, LLH_components);
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 2;
                LLH += normalDistLLH_debug(mu_a, kappa, alpha_a, Ncomponents, LLH_components);

                //Compute recusrive variables for next step
                omega_a = boundVariance(2*D_AB*dT(n));
                mu_a = (mu_a*zeta+kappa*eta_a)/alpha_a;
                eta_a = zeta*eta_a/alpha_a + omega_a;
                //                 std::cout<<" omega_a"<<omega_a<<" mu_a(n):"<<mu_a<<" eta_a(n):"<<eta_a<<std::endl;
            } else { //B->B
                epsilon_a = boundVariance(var_a(n)-D_AB*exposureT/3.);
                epsilon_b = boundVariance(var_b(n)-D_AB*exposureT/3.);
                beta = epsilon_b + rho_sq;
                b_pos = obs_b(n);
                gamma = epsilon_a + beta;
                kappa = (obs_a(n)*beta+b_pos*epsilon_a)/gamma;
                zeta = epsilon_a*beta/gamma;
                alpha_a = zeta + eta_a;
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 3;
                LLH += normalDistLLH_debug(obs_a(n), b_pos, gamma, Ncomponents, LLH_components);
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 3;
                LLH += normalDistLLH_debug(mu_a, kappa, alpha_a, Ncomponents, LLH_components);
                /*
                 *                std::cout<<"B->B: "<<"eps:"<<epsilon_a<<","<<epsilon_b<<" beta:"<<beta<<" b_pos:"<<b_pos
                 *                         <<" gamma:"<<gamma<<" kappa:"<<kappa<<" zeta:"<<zeta<<" alpha_a:"<<alpha_a
                 *                         <<" mu_a(n-1):"<<mu_a<<" eta_a(n-1):"<<eta_a;*/

                if(n==seq.end) break; //we are done recusing

                //Compute recusrive variables for next step
                omega_a = boundVariance(2*D_AB*dT(n));
                mu_a = (mu_a*zeta+kappa*eta_a)/alpha_a;
                eta_a = zeta*eta_a/alpha_a + omega_a;
                //                 std::cout<<" omega_a:"<<omega_a<<" mu_a(n):"<<mu_a<<" eta_a(n):"<<eta_a<<std::endl;
            }
        } else { //Current state is free
            if(last_state!=state) { //B->F  Carry in L_B(n-1) using mu_a, eta_a
                assert(n<seq.end); //cannot have a change point on last step
                epsilon_a = boundVariance(var_a(n)-D_A*exposureT/6.-D_AB*exposureT/6.);
                epsilon_b = boundVariance(var_b(n)-D_B*exposureT/6.-D_AB*exposureT/6.);
                //Add the extracted normal distributions to the LLH sum
                //A carries in L_B(n-1) = N(a_n, mu_a_{n-1}, eta_a_{n-1});
                //B does not have any relation to previous integral so "start fresh"
                alpha_a = epsilon_a + eta_a;
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 4;
                LLH += normalDistLLH_debug(obs_a(n), mu_a, alpha_a, Ncomponents, LLH_components);

                //Compute recusrsive variables for next step.
                mu_a = (mu_a*epsilon_a + obs_a(n)*eta_a)/alpha_a;
                omega_a = boundVariance(2*D_A*dT(n));
                eta_a = eta_a*epsilon_a/alpha_a + omega_a;

                mu_b = obs_b(n);
                omega_b = boundVariance(2*D_B*dT(n));
                eta_b = epsilon_b+omega_b;
            } else { //F->F
                epsilon_a = boundVariance(var_a(n)-D_A*exposureT/3.);
                epsilon_b = boundVariance(var_b(n)-D_B*exposureT/3.);

                //Add the extracted normal distributions to the LLH sum
                alpha_a = epsilon_a + eta_a;
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 5;
                LLH += normalDistLLH_debug(obs_a(n), mu_a, alpha_a, Ncomponents, LLH_components);

                alpha_b = epsilon_b + eta_b;
                LLH_components(0,Ncomponents) = n;
                LLH_components(1,Ncomponents) = 5;
                LLH += normalDistLLH_debug(obs_b(n), mu_b, alpha_b, Ncomponents, LLH_components);
                if(n==seq.end) break; //we are done recusing

                //Compute recusrive variables for next step
                mu_a = (mu_a*epsilon_a + obs_a(n)*eta_a)/alpha_a;
                omega_a = boundVariance(2*D_A*dT(n));
                eta_a = eta_a*epsilon_a/alpha_a + omega_a;

                mu_b = (mu_b*epsilon_b + obs_b(n)*eta_b)/alpha_b;
                omega_b = boundVariance(2*D_B*dT(n));
                eta_b = eta_b*epsilon_b/alpha_b + omega_b;
            }
        }
    }
    LLH *= -0.5;
    return LLH;
}


InteractionChangePoint::VecT
InteractionChangePoint::computeDVariance(FloatT D, const VecT &T) const
{
    int N = static_cast<int>(T.n_elem);
    VecT vD(N-1);
    for(int n=0; n<N-1; n++){
        vD(n) = 2*D*(T(n+1)-T(n));
        if(vD(n)<min_variance) vD(n)=min_variance;
    }
    return vD;
}

InteractionChangePoint::FloatT
InteractionChangePoint::normalDistLLH_debug(FloatT a, FloatT b, FloatT var, IndexT &Ncomponents, MatT& LLH_components) const
{
    FloatT val = log2pi + log(var) + square(a-b)/var;
    LLH_components(2,Ncomponents) = a;
    LLH_components(3,Ncomponents) = b;
    LLH_components(4,Ncomponents) = var;
    Ncomponents++;
    return val;
}

InteractionChangePoint::VecT
InteractionChangePoint::computeMVariance(FloatT D, const VecT &SE) const{
    int N = static_cast<int>(SE.n_elem);
    VecT vM = (SE%SE) - D*exposureT/3; //Variance is square of SE (standard error).
    for(int n=0; n<N; n++)
        if(fabs(vM(n))<min_variance) vM(n) = sgn(vM(n))<0 ? -min_variance : min_variance;
    return vM;
}
