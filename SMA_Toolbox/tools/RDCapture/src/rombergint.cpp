/** @file rombergint.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief Implementation of RombergInt numerical integration class
 */

#include "rombergint.h"
#include <cmath>
#include <cassert>


/**
 * 
 * @param f The function to integrate.  Should be double->double.
 * @param a The left bounds.
 * @param b The right bounds [assert: a < b]
 * @param order The order of the iteration.  Must be called in sequence [order=1,2,3,etc.]
 * @param old_sum The sum from order n-1.  If order==1, old_sum is ignored.  This should be exactly as
 *                returned from previous call
 * 
 * Uses the relation to recusively refine the estimate with 2^{order-2} new points.  Here order=n and h_n=(b-a)/2^{n-2} is the step size
 * $ I_n = h_n/h_{n-1} [ I_{n-1} + h_{n-1} \sum_{k=1}^{2^{n-2}} f_{n,k} ] $
 */
double RombergInt::extending_trapezoid_rule(RealFunc f, double a, double b,  int order, double old_sum)
{
    assert(a<b);
    assert(order>=1);
    if(order==1) return .5*(b-a)*(f(a)+f(b));
    long K= pow(2,order-2); //Numer of new points to evaluate;
    double h2 = (b-a)/K; // step size this order
    double h = h2/2; // step size for order-1
    double x = a+h; //start point
    double sum=0.;
    for(int i=0; i<K; i++) {
        sum += f(x); 
        x += h2;
    }
    double new_sum = h*(old_sum/h2+sum);
    double rel_error = fabs(old_sum - new_sum)/fabs(old_sum);
    printf("TRAP[order:%i, K:%li, old_sum:%.9g, new_sum:%.9g, rel_err:%.9g, stepsize:%.9g, h2:%.9g\n",order,K,old_sum,new_sum,rel_error,h,h2);
    return new_sum;
}

double RombergInt::extending_midpoint_rule(RealFunc f, double a, double b,  int order, double old_sum)
{
    assert(a<b);
    assert(order>=1);
    if(order==1) return (b-a)*f(.5*(a+b));
    long K= pow(3,order-2); //Numer of new points to evaluate;
    double h = (b-a)/(3*K); // step size this order
    double h2 = h*2; // step size for order-1
    double x = a+.5*h; //start point
    double sum=0.;
    for(int i=0; i<K; i++) {
        sum += f(x);
        x += h2;
        sum += f(x);
        x += h;
    }
    double new_sum = h*(old_sum/(3*h)+sum);
    double rel_error = fabs(old_sum - new_sum)/fabs(old_sum);
    printf("MIDP[order:%i, K:%li, old_sum:%.9g, new_sum:%.9g, rel_err:%.9g, stepsize:%.9g\n",order,K,old_sum,new_sum,rel_error,h);
    return new_sum;
}

double RombergInt::trapezoid(RealFunc f, double a, double b, int &Nevals, double eps, int maxIter)
{
    double sum, old_sum=0.;
    Nevals=0;
    for(int order=1; order<=maxIter; order++) {
        sum = extending_trapezoid_rule(f,a,b,order,old_sum);
        if(order>=5 && (fabs(sum-old_sum) < eps*fabs(old_sum) || (sum==0 && old_sum==0))) {
            Nevals = 1+pow(2,order-1);
            return sum;
        }
        old_sum=sum;
    }
    return NAN;//Failed convergence
}

double RombergInt::trapezoid_open(RealFunc f, double a, double b, int &Nevals, double eps, int maxIter)
{
    double sum, old_sum=0.;
    Nevals=0;
    for(int order=1; order<=maxIter; order++) {
        sum = extending_midpoint_rule(f,a,b,order,old_sum);
        if(order>=5 && (fabs(sum-old_sum) < eps*fabs(old_sum) || (sum==0 && old_sum==0))) {
            Nevals = pow(3,order-1);
            return sum;
        }
        old_sum=sum;
    }
    return NAN;//Failed convergence
}


double RombergInt::simpsons(RealFunc f, double a, double b, int &Nevals, double eps, int maxIter)
{
    double sum, old_sum=0.;
    double sum_t, old_sum_t=0.;
    Nevals=0;
    for(int order=1; order<=maxIter; order++) {
        sum_t = extending_trapezoid_rule(f,a,b,order,old_sum_t);
        sum = (4*sum_t - old_sum_t)/3;
        if(order>=5 && (fabs(sum-old_sum) < eps*fabs(old_sum) || (sum==0 && old_sum==0))) {
            Nevals = 1+pow(2,order-1); //2+2^(order-1)-1
            return sum;
        }
        old_sum = sum;
        old_sum_t = sum_t;
    }
    printf("Open Trapezoid Convergence Failed!");
    return NAN; //Failed convergence
}

double RombergInt::simpsons_open(RealFunc f, double a, double b, int &Nevals, double eps, int maxIter)
{
    double sum=0., old_sum=0.;
    double sum_t, old_sum_t=0.;
    Nevals=0;
    for(int order=1; order<=maxIter; order++) {
        sum_t = extending_midpoint_rule(f,a,b,order,old_sum_t);
        sum = (9*sum_t - old_sum_t)/8;
        if(order>=5 && (fabs(sum-old_sum) < eps*fabs(old_sum) || (sum==0 && old_sum==0))) {
            Nevals = pow(3,order-1);
            return sum;
        }
        old_sum = sum;
        old_sum_t = sum_t;
    }
    printf("Open Simpsons Convergence Failed!");
    return sum; //Failed convergence
}


double RombergInt::romberg(RealFunc f, double a, double b, int &Nevals, double eps, int Korder, int maxIter)
{
    VecT s(maxIter);
    VecT h(maxIter+1);
    double old_sum=0;
    double ss,dss;

    h(0)=1;
    Nevals=2;
    for(int j=1; j<=maxIter; j++) {
        s(j-1) = extending_trapezoid_rule(f,a,b,j,old_sum);
        old_sum=s(j-1);
        if(j>1) Nevals+=exp2(j-2);
        if (j>= Korder) {
            poly_interpolate(h(arma::span(j-Korder,j-1)),s(arma::span(j-Korder,j-1)),0.0,ss,dss);
            if (fabs(dss) <= eps*fabs(ss)) return ss;
        }
        h(j)=0.25*h(j-1);
    }
    printf("Open Romberg Convergence Failed!");
    return NAN; //Failed convergence
}

double RombergInt::romberg_open(RealFunc f, double a, double b, int &Nevals, double eps, int Korder, int maxIter)
{
    VecT s(maxIter);
    VecT h(maxIter+1);
    double old_sum=0;
    double ss,dss;
    
    h(0) = 1;
    Nevals = 2;
    for(int j=1; j<=maxIter; j++) {
        s(j-1) = extending_midpoint_rule(f,a,b,j,old_sum);
        old_sum = s(j-1);
        if (j >= Korder) {
            poly_interpolate(h(arma::span(j-Korder,j-1)), s(arma::span(j-Korder,j-1)), 0.0,ss,dss);
            if (fabs(dss) <= eps*fabs(ss)) {
                Nevals = pow(3,j-1);
                printf("sum:%.9f ss:%.9f dss:%.9f\n",s(j-1),ss,dss);
                return ss;
            }
        }
        h(j) = h(j-1)/9;
    }
    return NAN; //Failed convergence
}


void RombergInt::poly_interpolate(VecT x, VecT y, double x0, double &y0, double &y0_err)
{
    int N=static_cast<int>(x.n_elem);
    VecT c=y;
    VecT d=y;
    //Find the closest x
    int min_idx=arma::abs(x0-x).index_min();
    y0=y[min_idx--];
    for(int m=1; m<N; m++) {
        for(int i=0;i<N-m;i++) {
            double ho = x[i]-x0;
            double hp = x[i+m]-x0;
            double weight = c[i+1]-d[i];
            double den = ho-hp;
            assert(den!=0.0);
            double frac=weight/den;
            c[i]=ho*frac;
            d[i]=hp*frac;
        }
        if(2*(min_idx+1) < (N-m)){
            y0_err = c[min_idx+1];
        } else {
            y0_err = d[min_idx--];
        }
        y0 += y0_err;
    }
}
