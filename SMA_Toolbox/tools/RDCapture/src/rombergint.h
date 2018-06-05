/** @file rombergint.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief The class declaration and inline and templated functions for RombergInt.
 *
 * Integration with open and closed Newton-Coates adaptive quadrature rules using romberg polynomial interpolations.
 */

#ifndef _ROMBERGINT_H
#define _ROMBERGINT_H

#include <armadillo>
#include <limits> 

class RombergInt{
public:
    typedef double (RealFunc)(double);
    typedef arma::vec VecT;

    static const int DefaultMaxIter=24;
    static const int DefaultOpenMaxIter=20;
    
    inline static double trapezoid(RealFunc f, double a, double b);
    static double trapezoid(RealFunc f, double a, double b, int &Nevals, double eps=std::numeric_limits<double>::epsilon(), int maxIter=DefaultMaxIter);
    inline static double simpsons(RealFunc f, double a, double b);
    static double simpsons(RealFunc f, double a, double b, int &Nevals, double eps=std::numeric_limits<double>::epsilon(), int maxIter=DefaultMaxIter);
    inline static double romberg(RealFunc f, double a, double b);
    static double romberg(RealFunc f, double a, double b, int &Nevals, double eps=std::numeric_limits<double>::epsilon(), int Korder=5, int maxIter=DefaultMaxIter);
    
    inline static double trapezoid_open(RealFunc f, double a, double b);
    static double trapezoid_open(RealFunc f, double a, double b, int &Nevals, double eps=std::numeric_limits<double>::epsilon(), int maxIter=DefaultOpenMaxIter);
    inline static double simpsons_open(RealFunc f, double a, double b);
    static double simpsons_open(RealFunc f, double a, double b, int &Nevals, double eps=std::numeric_limits<double>::epsilon(), int maxIter=DefaultOpenMaxIter);
    inline static double romberg_open(RealFunc f, double a, double b);
    static double romberg_open(RealFunc f, double a, double b, int &Nevals, double eps=std::numeric_limits<double>::epsilon(), int Korder=5, int maxIter=DefaultOpenMaxIter);
    
    
    static void poly_interpolate(VecT xs, VecT ys, double x0, double &y0, double &y0_err);
    
private:
    static double extending_trapezoid_rule(RealFunc f, double a, double b,  int order, double old_sum);
    static double extending_midpoint_rule(RealFunc f, double a, double b,  int order, double old_sum);
};


/* Inlined Methods */
double RombergInt::trapezoid(RealFunc f, double a, double b)
{
    int Nevals;
    return trapezoid(f,a,b,Nevals);
}

double RombergInt::simpsons(RealFunc f, double a, double b)
{
    int Nevals;
    return simpsons(f,a,b,Nevals);
}

double RombergInt::romberg(RealFunc f, double a, double b)
{
    int Nevals;
    return romberg(f,a,b,Nevals);
}

double RombergInt::trapezoid_open(RealFunc f, double a, double b)
{
    int Nevals;
    return trapezoid_open(f,a,b,Nevals);
}

double RombergInt::simpsons_open(RealFunc f, double a, double b)
{
    int Nevals;
    return simpsons_open(f,a,b,Nevals);
}

double RombergInt::romberg_open(RealFunc f, double a, double b)
{
    int Nevals;
    return romberg_open(f,a,b,Nevals);
}


#endif /* ROMBERGINT_H */
