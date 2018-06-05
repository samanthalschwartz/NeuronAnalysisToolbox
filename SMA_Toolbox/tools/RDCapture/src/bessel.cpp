/** @file bessel.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 *
 *
 * Implementation of the NR modified bessel function approximations to I0 and K0.
 */

#include "bessel.h"
#include <armadillo>
#include "gsl/gsl_errno.h"

static const double ASSYMPTOTIC_THRESHOLD = 1e3;
static const double PROD_ASSYMPTOTIC_THRESHOLD = INFINITY;


void Bessel::gsl::error_handler(const char *reason, const char *file, int line, int gsl_errno)
{
    std::ostringstream msg;
    msg<<"Got GSL Error[#"<<gsl_errno<<"]: "<<reason<<" File:"<<file<<":"<<line<<"\n";
    std::cout<<msg.str();
    throw std::runtime_error(msg.str());
}

void Bessel::gsl::setupErrorHandler()
{
    gsl_set_error_handler(&Bessel::gsl::error_handler);
}

double Bessel::AS::I0(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return 1.;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<3.75 ) {
        return poly::I0_S(x); //[AS:9.8.1]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return exp(x)/sqrt(x) * poly::I0_L(x); //[AS:9.8.2]
    } else {
        return exp(x)/sqrt(2*arma::datum::pi*x) * poly::I0_A(x);//[AS:9.7.1]
    }
}

double Bessel::AS::I1(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return 0.;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<3.75 ) {
        return x * poly::I1_S(x); //[AS:9.8.3]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return exp(x)/sqrt(x) * poly::I1_L(x); //[AS:9.8.4]
    } else {
        return exp(x)/sqrt(2*arma::datum::pi*x) * poly::I1_A(x);//[AS:9.7.1]
    }
}

double Bessel::AS::K0(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<2) {
        return - log(x/2)*Bessel::AS::I0(x) + poly::K0_S(x); //[AS:9.8.5]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return exp(-x)/sqrt(x) * poly::K0_L(x); //[AS:9.8.6]
    } else {
        return exp(-x)/sqrt(arma::datum::pi/(2.*x)) * poly::K0_A(x);//[AS:9.7.2]
    }
}

double Bessel::AS::K1(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<2) {
        return log(x/2)*Bessel::AS::I1(x) + poly::K1_S(x) / x; //[AS:9.8.7]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return exp(-x)/sqrt(x) * poly::K1_L(x); //[AS:9.8.8]
    } else {
        return exp(-x)*sqrt(arma::datum::pi/(2.*x)) * poly::K0_A(x);//[AS:9.7.2]
    }
}

double Bessel::AS::logI0(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return 0.;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<3.75 ) {
        return log(poly::I0_S(x)); //[AS:9.8.1]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return x - .5*log(x) + log(poly::I0_L(x)); //[AS:9.8.2]
    } else {
        return x - .5*log(2*arma::datum::pi*x) + log(poly::I0_A(x));//[AS:9.7.1]
    }
}

double Bessel::AS::logI1(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return -INFINITY;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<3.75 ) {
        return log(x) + log(poly::I1_S(x)); //[AS:9.8.3]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return x -.5*log(x) + log(poly::I1_L(x)); //[AS:9.8.4]
    } else {
        return x - .5*log(2*arma::datum::pi*x) + log(poly::I1_A(x));//[AS:9.7.1]
    }
}

double Bessel::AS::logK0(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return -INFINITY;
    } else if (x<2) {
        return log(Bessel::K0(x)); //[AS:9.8.5]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return -x -.5*log(x) + log(poly::K0_L(x)); //[AS:9.8.6]
    } else {
        return -x +.5*log(arma::datum::pi/(2.*x)) +log(poly::K0_A(x));//[AS:9.7.2]
    }
}

double Bessel::AS::logK1(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return -INFINITY;
    } else if (x<2) {
        return log(Bessel::K1(x)); //[AS:9.8.7]
    } else if (x < ASSYMPTOTIC_THRESHOLD){
        return -x -.5*log(x) + log(poly::K1_L(x)); //[AS:9.8.8]
    } else {
        return -x + .5*log(arma::datum::pi/(2.*x)) + log(poly::K0_A(x));//[AS:9.7.2]
    }
}



double Bessel::AS::I0K0(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<3.75 ) {
        return Bessel::AS::I0(x) * Bessel::AS::K0(x); //Compute with the I0_S() polynomial
    } else if (x < PROD_ASSYMPTOTIC_THRESHOLD){
        return poly::I0_L(x) * poly::K0_L(x) / x; //[AS:9.8.2]
    } else {
        return poly::I0K0_A(x);//[AS:9.7.1]
    }
}


double Bessel::AS::I0K1(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<3.75 ) {
        return Bessel::AS::I0(x) * Bessel::AS::K1(x); //Compute with the I0_S() polynomial
    } else if (x < PROD_ASSYMPTOTIC_THRESHOLD){
        return poly::I0_L(x) * poly::K1_L(x) / x; //[AS:9.8.2]
    } else {
        return poly::I0K1_A(x);//[AS:9.7.1]
    }
}    

double Bessel::AS::I1K0(double x)
{
    if(x==0.){ 
        return 0;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<3.75 ) {
        return Bessel::AS::I1(x) * Bessel::AS::K0(x); //Compute with the I0_S() polynomial
    } else if (x < PROD_ASSYMPTOTIC_THRESHOLD){
        return poly::I1_L(x) * poly::K0_L(x) / x; //[AS:9.8.2]
    } else {
        return poly::I1K0_A(x);//[AS:9.7.1]
    }
}    

double Bessel::AS::I1K1(double x)
{
    if(x==0.){ 
        return .5;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<3.75 ) {
        return Bessel::AS::I1(x) * Bessel::AS::K1(x); //
    } else if (x < PROD_ASSYMPTOTIC_THRESHOLD){
        return poly::I1_L(x) * poly::K1_L(x) / x; //[AS:9.8.2]
    } else {
        return poly::I1K1_A(x);//[AS:9.7.1]
    }
}    

double Bessel::AS::I0K0(double x1, double x2)
{
    if(x1==0.) return Bessel::AS::K0(x2); //I0(0)==1
    else if(x2==0.) return INFINITY;  //K1(0)==INFINITY
    else if(std::isinf(x2)) return 0.; //K1(inf)==0
    else if(std::isinf(x1)) return INFINITY;
    else if (x1<3.75 || x2<2) return Bessel::AS::I0(x1) * Bessel::AS::K0(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::AS::poly::I0_L(x1)*Bessel::AS::poly::K0_L(x2); //[AS:9.8.2]
}


double Bessel::AS::I0K1(double x1, double x2)
{
    if(x1==0.) return Bessel::AS::K1(x2); //I0(0)==1
    else if(x2==0.) return INFINITY;  //K1(0)==INFINITY
    else if(std::isinf(x2)) return 0.; //K1(inf)==0
    else if(std::isinf(x1)) return INFINITY;
    else if (x1<3.75 || x2<2) return Bessel::AS::I0(x1) * Bessel::AS::K1(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::AS::poly::I0_L(x1)*Bessel::AS::poly::K1_L(x2); //[AS:9.8.2]
}

double Bessel::AS::I1K0(double x1, double x2)
{
    if(x1==0.) return 0.;
    else if (x2==0.) return INFINITY;
    else if (std::isinf(x2)) return 0.;
    else if (std::isinf(x1)) return INFINITY;
    else if (x1<3.75 || x2<2) return Bessel::AS::I1(x1) * Bessel::AS::K0(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::AS::poly::I1_L(x1)*Bessel::AS::poly::K0_L(x2); //[AS:9.8.2]
}    

double Bessel::AS::I1K1(double x1, double x2)
{
    if(x1==0.) return 0.;
    else if (x2==0.) return INFINITY;
    else if (std::isinf(x2)) return 0.;
    else if (std::isinf(x1)) return INFINITY;
    else if (x1<3.75 || x2<2) return Bessel::AS::I1(x1) * Bessel::AS::K1(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::AS::poly::I1_L(x1)*Bessel::AS::poly::K1_L(x2); //[AS:9.8.2]
}    

double Bessel::NR::I0(double x)
{
    double ax = fabs(x);
    double y,val;
    if(ax < 3.75) {
        y = x/3.75;
        y *= y; // y= (x/3.75)^2
        val = 1. + y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    } else  {
        y = 3.75/ax;
        val = 0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2)))))));
        val *= exp(ax)/sqrt(ax);
    }
    return val;
}

double Bessel::NR::K0(double x)
{
    assert(x>=0);
    double y,val;
    if(x==0) {
        val = INFINITY;
    } else if(x < 2) {
        y = x*x/4;
        val = -0.57721566 + y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5)))));
        val += -log(x/2)*Bessel::NR::I0(x);
    } else {
        y = 2/x;
        val = 1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3)))));
        val *= exp(-x)/sqrt(x);
    }
    return val;
}

double Bessel::NR::I1(double x)
{
    double ax = fabs(x);
    double y,val;
    if(ax < 3.75) {
        y = x/3.75;
        y *= y; // y= (x/3.75)^2
        val = ax*(0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934 + y*(0.2658733e-1
                    + y*(0.301532e-2 + y*0.32411e-3))))));
    } else {
        y = 3.75/ax;
        val = 0.2282967e-1 + y*(-0.2895312e-1 + y*(0.1787654e-1 - y*0.420059e-2));
        val = 0.39894228 + y*(-0.3988024e-1 + y*(-0.362018e-2 + y*(0.163801e-2 + y*(-0.103155e-1 + y*val))));
        val *= exp(ax)/sqrt(ax);
    }
    return std::fabs(val);
}


double Bessel::NR::K1(double x)
{
    assert(x>=0);
    double y,val;
    if(x==0) {
        val = INFINITY;
    } else if(x < 2) {
        y = x*x/4;
        val = log(x/2)*Bessel::NR::I1(x);
        val += (1 + y*(0.15443144 + y*(-0.67278579 + y*(-0.18156897 + y*(-0.1919402e-1 
                + y*(-0.110404e-2+y*-0.4686e-4))))))/x;
    } else {
        y = 2/x;
        val  = exp(-x)/sqrt(x);
        val *= 1.25331414 + y*(0.23498619 + y*(-0.3655620e-1 + y*(0.1504268e-1 + y*(-0.780353e-2 
                        + y*(0.325614e-2 + y*-0.68245e-3)))));
    }
    return val;
}



double Bessel::BE::I0(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return 1.;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<=15) {
        return Bessel::BE::poly::I0_1(x);
    } else {
        return exp(x)/sqrt(x) * Bessel::BE::poly::I0_2(x);
    }
}

double Bessel::BE::I1(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return 0.;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<=15) {
        return x * Bessel::BE::poly::I1_1(x);
    } else {
        return exp(x)/sqrt(x) * Bessel::BE::poly::I1_2(x);
    }
}

double Bessel::BE::K0(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=1) {
        return Bessel::BE::poly::K0_1(x) - log(x) * Bessel::BE::poly::K0_2(x);
    } else {
        return exp(-x)/sqrt(x) * Bessel::BE::poly::K0_3(x);
    }
}

double Bessel::BE::K1(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=1) {
        return (Bessel::BE::poly::K1_1(x) + log(x) * Bessel::BE::poly::K1_2(x))/x;
    } else {
        return exp(-x)/sqrt(x) * Bessel::BE::poly::K1_3(x);
    }
}


double Bessel::BE::logI0(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return 0.;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<=15) {
        return log(Bessel::BE::poly::I0_1(x));
    } else {
        return x-.5*log(x) + log(Bessel::BE::poly::I0_2(x));
    }
}

double Bessel::BE::logI1(double x)
{
    x=fabs(x);
    if(x==0.){ 
        return -INFINITY;
    } else if (std::isinf(x)) {
        return INFINITY;
    } else if (x<=15) {
        return log(x) + log(Bessel::BE::poly::I1_1(x)); 
    } else {
        return x-.5*log(x) + log(Bessel::BE::poly::I1_2(x));
    }
}

double Bessel::BE::logK0(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return -INFINITY;
    } else if (x<=1) {
        return log(Bessel::BE::poly::K0_1(x) - log(x) * Bessel::BE::poly::K0_2(x));
    } else {
        return -x-.5*log(x) + log(Bessel::BE::poly::K0_3(x));
    }
}

double Bessel::BE::logK1(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=1) {
        double logx = log(x);
        return log(Bessel::BE::poly::K1_1(x) - logx * Bessel::BE::poly::K1_2(x))-logx;
    } else {
        return -x-.5*log(x) + log(Bessel::BE::poly::K1_3(x));
    }
}


double Bessel::BE::I0K0(double x)
{
    if(x==0.){ 
        return INFINITY;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=15.) {
        return Bessel::BE::I0(x) * Bessel::BE::K0(x); //
    } else {
        return Bessel::BE::poly::I0_2(x) * Bessel::BE::poly::K0_3(x) / x; 
    }
}    

double Bessel::BE::I0K1(double x)
{
    if(x==0.){ 
        return .5;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=15.) {
        return Bessel::BE::I0(x) * Bessel::BE::K1(x); //
    } else {
        return Bessel::BE::poly::I0_2(x) * Bessel::BE::poly::K1_3(x) / x; 
    }
}    

double Bessel::BE::I1K0(double x)
{
    if(x==0.){ 
        return 0;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=15.) {
        return Bessel::BE::I1(x) * Bessel::BE::K0(x); //
    } else {
        return Bessel::BE::poly::I1_2(x) * Bessel::BE::poly::K0_3(x) / x; 
    }
}    


double Bessel::BE::I1K1(double x)
{
    if(x==0.){ 
        return .5;
    } else if (std::isinf(x)) {
        return 0;
    } else if (x<=15.) {
        return Bessel::BE::I1(x) * Bessel::BE::K1(x); //
    } else {
        return Bessel::BE::poly::I1_2(x) * Bessel::BE::poly::K1_3(x) / x; 
    }
}    

double Bessel::BE::I0K0(double x1, double x2)
{
    if(x1==0.) return Bessel::BE::K0(x2); //I0(0)==1
    else if(x2==0.) return INFINITY;  //K1(0)==INFINITY
    else if(std::isinf(x2)) return 0;
    else if(std::isinf(x1)) return INFINITY;
    else if (x1<=15 || x2<=1) return Bessel::BE::I0(x1) * Bessel::BE::K0(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::BE::poly::I0_2(x1)*Bessel::BE::poly::K0_3(x2); //[AS:9.8.2]
}


double Bessel::BE::I0K1(double x1, double x2)
{
    if(x1==0.) return Bessel::BE::K1(x2); //I0(0)==1
    else if(x2==0.) return INFINITY;  //K1(0)==INFINITY
    else if(std::isinf(x2)) return 0;
    else if(std::isinf(x1)) return INFINITY;
    else if (x1<=15 || x2<=1) return Bessel::BE::I0(x1) * Bessel::BE::K1(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::BE::poly::I0_2(x1)*Bessel::BE::poly::K1_3(x2); //[AS:9.8.2]
}

double Bessel::BE::I1K0(double x1, double x2)
{
    if(x1==0.) return 0.;
    else if (x2==0.) return INFINITY;
    else if (std::isinf(x2)) return 0.;
    else if (std::isinf(x1)) return INFINITY;
    else if (x1<=15 || x2<=1) return Bessel::BE::I1(x1) * Bessel::BE::K0(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::BE::poly::I1_2(x1)*Bessel::BE::poly::K0_3(x2); //[AS:9.8.2]
}    

double Bessel::BE::I1K1(double x1, double x2)
{
    if(x1==0. && x2==0.) return .5;
    else if (x2==0.) return INFINITY;
    else if (std::isinf(x2)) return 0.;
    else if (std::isinf(x1)) return INFINITY;
    else if (x1<=15 || x2<=1) return Bessel::BE::I1(x1) * Bessel::BE::K1(x2); //Compute with the I0_S() polynomial
    else return exp(x1-x2)/sqrt(x1*x2) * Bessel::BE::poly::I1_2(x1) * Bessel::BE::poly::K1_3(x2); //[AS:9.8.2]
}


// double Bessel::BE::I0K0(double x)
// 
// 
// 
// 
// 
// double Bessel::BE::I0K1(double x)
// {
//     
// }
// 
// double Bessel::BE::I1K0(double x)
// {
//     
// }
// 
// double Bessel::BE::I1K1(double x)
// {
//     
// }
// 
// double Bessel::BE::I0K1(double x1, double x2);
// double Bessel::BE::I1K0(double x1, double x2);
// 




/* Boost polynomial coefficients */
/* I0 */
const double Bessel::BE::poly::I0_P1[] = {
    -2.2335582639474375249e+15,
    -5.5050369673018427753e+14,
    -3.2940087627407749166e+13,
    -8.4925101247114157499e+11,
    -1.1912746104985237192e+10,
    -1.0313066708737980747e+08,
    -5.9545626019847898221e+05,
    -2.4125195876041896775e+03,
    -7.0935347449210549190e+00,
    -1.5453977791786851041e-02,
    -2.5172644670688975051e-05,
    -3.0517226450451067446e-08,
    -2.6843448573468483278e-11,
    -1.5982226675653184646e-14,
    -5.2487866627945699800e-18
};
const double Bessel::BE::poly::I0_Q1[] = {
    -2.2335582639474375245e+15,
    7.8858692566751002988e+12,
    -1.2207067397808979846e+10,
    1.0377081058062166144e+07,
    -4.8527560179962773045e+03,
    1.0
};
const double Bessel::BE::poly::I0_P2[] = {
    -2.2210262233306573296e-04,
    1.3067392038106924055e-02,
    -4.4700805721174453923e-01,
    5.5674518371240761397e+00,
    -2.3517945679239481621e+01,
    3.1611322818701131207e+01,
    -9.6090021968656180000e+00
};
const double Bessel::BE::poly::I0_Q2[] = {
    -5.5194330231005480228e-04,
    3.2547697594819615062e-02,
    -1.1151759188741312645e+00,
    1.3982595353892851542e+01,
    -6.0228002066743340583e+01,
    8.5539563258012929600e+01,
    -3.1446690275135491500e+01,
    1.0
};

/* I1 */
const double Bessel::BE::poly::I1_P1[] = {
    -1.4577180278143463643e+15,
    -1.7732037840791591320e+14,
    -6.9876779648010090070e+12,
    -1.3357437682275493024e+11,
    -1.4828267606612366099e+09,
    -1.0588550724769347106e+07,
    -5.1894091982308017540e+04,
    -1.8225946631657315931e+02,
    -4.7207090827310162436e-01,
    -9.1746443287817501309e-04,
    -1.3466829827635152875e-06,
    -1.4831904935994647675e-09,
    -1.1928788903603238754e-12,
    -6.5245515583151902910e-16,
    -1.9705291802535139930e-19,
};
const double Bessel::BE::poly::I1_Q1[] = {
    -2.9154360556286927285e+15,
    9.7887501377547640438e+12,
    -1.4386907088588283434e+10,
    1.1594225856856884006e+07,
    -5.1326864679904189920e+03,
    1.0,
};
const double Bessel::BE::poly::I1_P2[] = {
    1.4582087408985668208e-05,
    -8.9359825138577646443e-04,
    2.9204895411257790122e-02,
    -3.4198728018058047439e-01,
    1.3960118277609544334e+00,
    -1.9746376087200685843e+00,
    8.5591872901933459000e-01,
    -6.0437159056137599999e-02,
};
const double Bessel::BE::poly::I1_Q2[] = {
    3.7510433111922824643e-05,
    -2.2835624489492512649e-03,
    7.4212010813186530069e-02,
    -8.5017476463217924408e-01,
    3.2593714889036996297e+00,
    -3.8806586721556593450e+00,
    1.0,
};

/* K0 */

const double Bessel::BE::poly::K0_P1[] = {
    2.4708152720399552679e+03,
    5.9169059852270512312e+03,
    4.6850901201934832188e+02,
    1.1999463724910714109e+01,
    1.3166052564989571850e-01,
    5.8599221412826100000e-04
};
const double Bessel::BE::poly::K0_Q1[] = {
    2.1312714303849120380e+04,
    -2.4994418972832303646e+02,
    1.0
};
const double Bessel::BE::poly::K0_P2[] = {
    -1.6128136304458193998e+06,
    -3.7333769444840079748e+05,
    -1.7984434409411765813e+04,
    -2.9501657892958843865e+02,
    -1.6414452837299064100e+00
};
const double Bessel::BE::poly::K0_Q2[] = {
    -1.6128136304458193998e+06,
    2.9865713163054025489e+04,
    -2.5064972445877992730e+02,
    1.0
};
const double Bessel::BE::poly::K0_P3[] = {
    1.1600249425076035558e+02,
    2.3444738764199315021e+03,
    1.8321525870183537725e+04,
    7.1557062783764037541e+04,
    1.5097646353289914539e+05,
    1.7398867902565686251e+05,
    1.0577068948034021957e+05,
    3.1075408980684392399e+04,
    3.6832589957340267940e+03,
    1.1394980557384778174e+02
};
const double Bessel::BE::poly::K0_Q3[] = {
    9.2556599177304839811e+01,
    1.8821890840982713696e+03,
    1.4847228371802360957e+04,
    5.8824616785857027752e+04,
    1.2689839587977598727e+05,
    1.5144644673520157801e+05,
    9.7418829762268075784e+04,
    3.1474655750295278825e+04,
    4.4329628889746408858e+03,
    2.0013443064949242491e+02,
    1.0
};

/* K1 */
const double Bessel::BE::poly::K1_P1[] = {
    -2.2149374878243304548e+06,
    7.1938920065420586101e+05,
    1.7733324035147015630e+05,
    7.1885382604084798576e+03,
    9.9991373567429309922e+01,
    4.8127070456878442310e-01
};

const double Bessel::BE::poly::K1_Q1[] = {
    -2.2149374878243304548e+06,
    3.7264298672067697862e+04,
    -2.8143915754538725829e+02,
    1.0
};

const double Bessel::BE::poly::K1_P2[] = {
    0.0,
    -1.3531161492785421328e+06,
    -1.4758069205414222471e+05,
    -4.5051623763436087023e+03,
    -5.3103913335180275253e+01,
    -2.2795590826955002390e-01
};

const double Bessel::BE::poly::K1_Q2[] = {
    -2.7062322985570842656e+06,
    4.3117653211351080007e+04,
    -3.0507151578787595807e+02,
    1.0
};

const double Bessel::BE::poly::K1_P3[] = {
    2.2196792496874548962e+00,
    4.4137176114230414036e+01,
    3.4122953486801312910e+02,
    1.3319486433183221990e+03,
    2.8590657697910288226e+03,
    3.4540675585544584407e+03,
    2.3123742209168871550e+03,
    8.1094256146537402173e+02,
    1.3182609918569941308e+02,
    7.5584584631176030810e+00,
    6.4257745859173138767e-02
};

const double Bessel::BE::poly::K1_Q3[] = {
    1.7710478032601086579e+00,
    3.4552228452758912848e+01,
    2.5951223655579051357e+02,
    9.6929165726802648634e+02,
    1.9448440788918006154e+03,
    2.1181000487171943810e+03,
    1.2082692316002348638e+03,
    3.3031020088765390854e+02,
    3.6001069306861518855e+01,
    1.0
};
