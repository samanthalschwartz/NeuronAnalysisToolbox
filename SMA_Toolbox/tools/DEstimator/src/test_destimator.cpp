/**
 * @file test_destimator.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-06-2014
 */

#include "destimator.h"
// #include "Optimizer1D.h"
// #include <boost/bind.hpp>

using namespace arma;
using namespace std;

double func(double x)
{
    return  log(pow(x,log(x))+sqrt(1+pow(x-17,6)));
}
/*
void testOptimizer()
{
    Optimizer1D<double> opt(func,50);
    double xmin, Fmin;
    opt.minimize(0,1.1,xmin,Fmin);
    printf("xmin: %.16g Fmin: %.16g\n",xmin,Fmin);
    arma::vec X,F;
    opt.getStats(X,F);
    int N=X.n_elem;
    for(int i=0; i<N; i++){
        printf("%i: F(%.16g)=%.16g\n",i,X(i),F(i));
    }
}*/

void testDEstimator()
{
    typedef double FloatT;
    int N=70;
    FloatT dt=1.;
    FloatT D=1.;
    FloatT scalarV=10.0;
    FloatT exposureT=dt;
    Col<FloatT> Obs(N), T(N), V(N);
    int ND=20;
    Col<FloatT> Ds(ND);
    for(int n=0;n<N;n++) T(n)=n*dt;
    Obs.randn();
    Obs*=sqrt(2*D*dt);
    Obs=arma::cumsum(Obs);
    V.fill(scalarV);
    for(int n=0;n<ND;n++) Ds(n)=10*double(n)/ND;

    Col<FloatT> recursiveLLH(ND), laplaceLLH(ND), markovLLH(ND);
    cout<<"T: "<<T.t()<<endl;
    cout<<"Obs: "<<Obs.t()<<endl;
    cout<<"V: "<<V.t()<<endl;
    cout<<"exposureT: "<<exposureT<<endl;
    cout<<"Ds: "<<Ds.t()<<endl;

    DEstimator<FloatT>::LLH_recursive1D(Ds,Obs,T,V,exposureT,recursiveLLH);
    DEstimator<FloatT>::LLH_laplace1D(Ds,Obs,T,V,exposureT,laplaceLLH);
    DEstimator<FloatT>::LLH_markov1D(Ds,Obs,T,V,exposureT,markovLLH);
    cout<<"Recurisve LLH: "<<recursiveLLH.t()<<endl;
    cout<<"Laplace LLH: "<<laplaceLLH.t()<<endl;
    cout<<"Markov LLH: "<<markovLLH.t()<<endl;
    Mat<FloatT> Obs2(N,2), V2(N,2);
    Obs2.col(0)=Obs;
    Obs2.col(1)=Obs;
    V2.col(0)=V;
    V2.col(1)=V;
    DEstimator<FloatT> est(Obs2,T,V2,exposureT);
    cout<<"LLH(0): "<<est.LLH(0)<<"\n";
//     FloatT mleD, mllh;
//     int numEvals=est.MLEdim(2,mleD,mllh);
//     printf("MLE[#evals:%i]: LLH(%.16g)=%.16g\n",numEvals,mleD,mllh);
}


int main(){
//     testOptimizer();
    testDEstimator();
    return 0;
}
