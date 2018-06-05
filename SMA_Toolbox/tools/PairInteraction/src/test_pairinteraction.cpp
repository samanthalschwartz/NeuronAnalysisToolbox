


#include <armadillo>
#include "pairinteractionsmc.h"
#include "unscentedtransform.h"
using namespace arma;
using namespace std;
using namespace pair_int;


void testData()
{
    PairInteractionSMC::ParamsT params;
    params.emplace("DA",0.1);
    params.emplace("DB",0.1);
    params.emplace("DAB",0.05);
    params.emplace("Dphi",0.2);
    params.emplace("rho",0.05);
    params.emplace("sigma_rho",0.015);
    params.emplace("gamma",1.0e3);
    params.emplace("rho_bind",0.05);
    params.emplace("lambda_bind",1);
    params.emplace("koff",1);
    
    PairInteractionSMC::MatListT d1,d2;
    d1.emplace_back(9,10,fill::zeros);
    d1.emplace_back(9,12,fill::zeros);
    d2.emplace_back(9,3,fill::zeros);
    d2.emplace_back(9,4,fill::zeros);
    
    PairInteractionSMC smc(params,d1);
    cout<<"Ndata: "<<smc.addData(d2)<<endl;
}


void testUT()
{
    vec mean{1,2};
    mat cov{ {1,0},{0,1}};
    UnscentedTransform ut(2);
    vec tmean;
    mat tcov;
    ut.transform([](const vec &m){return arma::exp(m);}, mean, cov, tmean, tcov);
    cout<<"Mean: "<<mean.t()<<"Cov: \n"<<cov<<endl;
    cout<<"TMean: "<<tmean.t()<<"TCov: \n"<<tcov<<endl;
    cout<<"LLH: "<<multiGaussianLLH(tmean,tcov)<<"\n";
    
}

void testSimulate()
{
    PairInteractionSMC::ParamsT params;
    params.emplace("DA",0.1);
    params.emplace("DB",0.1);
    params.emplace("DAB",0.05);
    params.emplace("Dphi",50);
    params.emplace("rho",0.05);
    params.emplace("sigma_rho",0.015);
    params.emplace("gamma",1.0e3);
    params.emplace("rho_bind",0.05);
    params.emplace("lambda_bind",1);
    params.emplace("koff",1);
    
    PairInteractionSMC smc(params);
    smc.setPriorState(vec{0.0,1});
    mat cov(2,2,fill::eye);
    smc.setPriorPositionGaussian(vec{0,0},cov);
    params = smc.getParams();
    cube obs,particles,llh;
    double tmax = 0.1;
    smc.simulate(1,10,tmax,"Ancestral",params,obs,particles,llh);
    cout<<"obs(1)\n"<<obs.slice(0)<<"\n";
    cout<<"particles(1)\n"<<particles.slice(0)<<"\n";
    cout<<"llh(1)\n"<<vec{arma::sum(arma::sum(llh,1),0)}.t()<<"\n";
    smc.addData(obs);
    VecT trueLLH = smc.computeLLH(particles);
    cout<<"TrueLLH: "<<trueLLH.t()<<"\n";
    
    smc.runParticleFilter(50,"Combined");
    VecT transLLH = smc.sampleParticleLLH();
    cout<<"TransLLH: "<<transLLH.t()<<"\n";
    smc.runParticleFilter(500,"Observation");
    VecT obsLLH = smc.sampleParticleLLH();
    cout<<"ObsLLH: "<<obsLLH.t()<<"\n";
//     VecT trueLLH = smc.computeLLH(particles);
    cout<<"TrueLLH: "<<trueLLH.t()<<"\n";
//     PairInteractionSMC::MatListT ps;
//     PairInteractionSMC::VecT ps_llh;
//     for(IdxT n=1; n<10;n++) {
//         smc.sampleParticle(ps, ps_llh);
//         cout<<"Sampled particle (LLH="<<ps_llh(0)<<"\n";
//         cout<<"All data sampled LLH: "<< smc.sampleParticleLLH().t()<<"\n";
//     }
//     PairInteractionSMC::CubeListT P;
//     PairInteractionSMC::MatT Pweights, Pllh;
//     smc.getAllParticles(P, Pweights, Pllh); //For debugging
//     cout<<"Particle 1:\n"<<P[0].slice(0)<<"\n";
//     cout<<"LLH: "<<Pllh(0,0)<<" weight:"<<Pweights(0,0)<<"\n";
}

int main()
{
    testSimulate();
//     testData();
//     testUT();
    return 0;
}
