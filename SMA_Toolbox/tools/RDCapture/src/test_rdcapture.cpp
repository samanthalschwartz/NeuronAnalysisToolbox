
#include <armadillo>
#include "bessel.h"
#include "rombergint.h"
#include "laplacetform.h"
#include "rdcapture.h"
// #include <random>

using namespace arma;
using namespace std;

double pi=3.1415926;
double xcircle(double x){ return sqrt(1-x*x);}
double xgauss(double x){ return 1/sqrt(2*pi)* exp(-x*x/2);}
double xsin(double x){ return sin(x);}
double xI0(double x) { return Bessel::AS::I0(x);}
double xK0(double x) { return Bessel::AS::K0(x);}
double xnr(double x) { return pow(x,4)*log(x+sqrt(x*x+1));}

double sexp(double s) {return 1./(s-(-1));}

void testIntr(double (f)(double),double a,double b,double eps,const char* str)
{
    double I;
    int Nevals;
    I=RombergInt::trapezoid(f,a,b,Nevals,eps);
    printf("[T] %s[%f,%f] Integral: %.9g  Nevals: %i\n",str,a,b,I,Nevals);
    I=RombergInt::simpsons(f,a,b,Nevals,eps);
    printf("[S] %s[%f,%f] Integral: %.9g  Nevals: %i\n",str,a,b,I,Nevals);
    I=RombergInt::romberg(f,a,b,Nevals,eps,5);
    printf("[R] %s[%f,%f] Integral: %.9g  Nevals: %i\n",str,a,b,I,Nevals);
}

void testOpenIntr(double (f)(double),double a,double b,double eps,const char* str)
{
    double I;
    int Nevals;
    I=RombergInt::trapezoid_open(f,a,b,Nevals,eps);
    printf("[T] %s[%f,%f] Integral: %.9g  Nevals: %i\n",str,a,b,I,Nevals);
    I=RombergInt::simpsons_open(f,a,b,Nevals,eps);
    printf("[S] %s[%f,%f] Integral: %.9g  Nevals: %i\n",str,a,b,I,Nevals);
    I=RombergInt::romberg_open(f,a,b,Nevals,eps,6);
    printf("[R] %s[%f,%f] Integral: %.9g  Nevals: %i\n",str,a,b,I,Nevals);
}


void testBessel()
{
    int nPoints = 150;
    double delta = 1.2;
    Bessel::gsl::setupErrorHandler();
    vec I0_AS(nPoints);
    vec K0_AS(nPoints);
    vec I0_NR(nPoints);
    vec K0_NR(nPoints);
    vec I0_boost(nPoints);
    vec K0_boost(nPoints);
    vec I0_gsl(nPoints);
    vec K0_gsl(nPoints);
    vec I1_AS(nPoints);
    vec K1_AS(nPoints);
    vec I1_NR(nPoints);
    vec K1_NR(nPoints);
    vec I1_boost(nPoints);
    vec K1_boost(nPoints);
    vec I1_gsl(nPoints);
    vec K1_gsl(nPoints);
    double x=1e-6;
    for(int n=0; n<nPoints; n++){
        if(n==nPoints-2) {
            x=0.;
        } else if(n==nPoints-1){
            x=INFINITY;
        } else {
            x*=delta;
        }
        I0_AS[n] = Bessel::AS::I0(x);
        K0_AS[n] = Bessel::AS::K0(x);
        I0_NR[n] = Bessel::NR::I0(x);
        K0_NR[n] = Bessel::NR::K0(x);
        I0_boost[n] = Bessel::boost::I0(x);
        K0_boost[n] = Bessel::boost::K0(x);
        I0_gsl[n] = Bessel::gsl::I0(x);
        K0_gsl[n] = Bessel::gsl::K0(x);
        
        I1_AS[n] = Bessel::AS::I1(x);
        K1_AS[n] = Bessel::AS::K1(x);
        I1_NR[n] = Bessel::NR::I1(x);
        K1_NR[n] = Bessel::NR::K1(x);
        I1_boost[n] = Bessel::boost::I1(x);
        K1_boost[n] = Bessel::boost::K1(x);
        I1_gsl[n] = Bessel::gsl::I1(x);
        K1_gsl[n] = Bessel::gsl::K1(x);
        printf("I0(%g): %g %g %g %g\n",x,I0_AS[n],I0_NR[n],I0_boost[n],I0_gsl[n]);
        printf("K0(%g): %g %g %g %g\n",x,K0_AS[n],K0_NR[n],K0_boost[n],K0_gsl[n]);
        printf("I1(%g): %g %g %g %g\n",x,I1_AS[n],I1_NR[n],I1_boost[n],I1_gsl[n]);
        printf("K1(%g): %g %g %g %g\n",x,K1_AS[n],K1_NR[n],K1_boost[n],K1_gsl[n]);
//         x*=delta;
    }
}
    
//     double D=0.1;
//     double rho=0.05;
//     double lambda=10;
//     double a0=0.1;
//     const int N = 1e6;
//     vec t = linspace(1e-3,10,N);
//     RDCapture rd(D,rho,lambda);
//     vec Q(N);
//     rd.survivalProb_parallel(Q,a0,t);
//     for(int i=0;i<static_cast<int>(t.n_elem);i+=static_cast<int>(t.n_elem)/10){
//         printf("Q(%.9f)=%.9f\n",t(i),Q(i));
//     }
//     
// }


void testRDCapture()
{
    double D=0.1;     // um^2/s
    double rho =0.05; // um
    double lambda=1e6; // 1/s
    
    RDCapture rd(D,rho,lambda);
    double r0=0.01; //um
    double t=1.0; //s
    double mu = rd.mu(r0,t);
    double nu = rd.nu(t);
    double Q = rd.survivalProb(r0,t);
    cout<<setprecision(16);
    cout<<"D:="<<D<<" rho:="<<rho<<" lambda:="<<lambda<<" r0:="<<r0<<" t:="<<t<<endl;
    cout<<"mu:="<<mu<<" nu:="<<nu<<" Q:="<<Q<<endl;
    int N=1e2;
    vec xs = logspace(-4,3,N);
    for(int n=0; n<N; n++){
        double I0=Bessel::AS::I0(xs(n));
        double I1=Bessel::AS::I1(xs(n));
        double K0=Bessel::AS::K0(xs(n));
        double K1=Bessel::AS::K1(xs(n));
        double logI0=Bessel::AS::logI0(xs(n));
        double logI1=Bessel::AS::logI1(xs(n));
        double logK0=Bessel::AS::logK0(xs(n));
        double logK1=Bessel::AS::logK1(xs(n));
        double I0K0=Bessel::AS::I0K0(xs(n));
        double I0K1=Bessel::AS::I0K1(xs(n));
        double I1K0=Bessel::AS::I1K0(xs(n));
        double I1K1=Bessel::AS::I1K1(xs(n));
        cout<<"x:"<<xs(n)<<" I0:"<<I0<<" K0:"<<K0<<" I1:"<<I1<<" K1:"<<K1<<"\n";
        cout<<"logI0:"<<logI0<<" log(I0):"<<log(I0)<<"\n";
        cout<<"logI1:"<<logI1<<" log(I1):"<<log(I1)<<"\n";
        cout<<"logK0:"<<logK0<<" log(K0):"<<log(K0)<<"\n";
        cout<<"logK1:"<<logK1<<" log(K1):"<<log(K1)<<"\n";
        cout<<"I0K0:"<<I0K0<<" log(I0K0):"<<I0*K0<<"\n";
        cout<<"I0K1:"<<I0K1<<" log(I0K1):"<<I0*K1<<"\n";
        cout<<"I1K0:"<<I1K0<<" log(I1K0):"<<I1*K0<<"\n";
        cout<<"I1K1:"<<I1K1<<" log(I1K1):"<<I1*K1<<"\n\n";
    }
    vec vals(N);
    rd.survivalProb_parallel(vals, xs,0.05);
    for(int i=0; i<N; i++) {
        cout<<"x:"<<xs(i)<<" survivalProb:"<<vals(i)<<"\n";
        if(std::isnan(vals(i))) throw std::runtime_error("Got a NAN");
    }
    std::cout<<vals.t();
}

int main(){
//     testBessel();
    testRDCapture();
    return 0;
}
