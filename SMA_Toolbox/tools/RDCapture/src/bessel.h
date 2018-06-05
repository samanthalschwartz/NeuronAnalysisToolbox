/** @file bessel.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief Bessel functions.
 *
 * This class computes the modified bessel functions of the first and second kind of order 0.
 * These are known as I0() and K0() functions.  
 * 
 * We provide a common interface to the boost and GSL implementations as well as the perhaps faster
 * Numerical Recipes implementation.  We will be substituting various definitions to compare the trade-off of
 * accuracy with speed.
 * 
 *  We use the method of Abromowitz & Stegun [1964].
 */

#ifndef _BESSEL_H
#define _BESSEL_H

#include <utility>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/tools/rational.hpp>
#include <gsl/gsl_sf_bessel.h>

namespace Bessel {
/* Helper methods */
inline double square(double x) {return x*x;}
inline double logSum(double log_a, double log_b) 
{
    if (log_a<log_b) std::swap(log_a,log_b);
    double delta = log_a-log_b;
    if (delta > 20) return log_a + log1p(exp(-delta));
    else return log_a + log(1+exp(-delta));
}



/* Boost implementations */
namespace boost{
    inline double I0(double x) { try {return ::boost::math::cyl_bessel_i(0,x);} catch (...) {return NAN;}}
    inline double K0(double x) { try {return ::boost::math::cyl_bessel_k(0,x);} catch (...) {return NAN;}}
    inline double I1(double x) { try {return ::boost::math::cyl_bessel_i(1,x);} catch (...) {return NAN;}}
    inline double K1(double x) { try {return ::boost::math::cyl_bessel_k(1,x);} catch (...) {return NAN;}}
    
    
}

/* GSL implementations */
namespace gsl {
    inline double I0(double x) { try {return gsl_sf_bessel_I0(x);} catch (const std::exception &e) {return NAN;} }
    inline double K0(double x) { try {return gsl_sf_bessel_K0(x);} catch (const std::exception &e) {return NAN;} }
    inline double I1(double x) { try {return gsl_sf_bessel_I1(x);} catch (const std::exception &e) {return NAN;} }
    inline double K1(double x) { try {return gsl_sf_bessel_K1(x);} catch (const std::exception &e) {return NAN;} }
    void error_handler (const char *reason, const char *file, int line, int gsl_errno);
    void setupErrorHandler();
}

/* NR implementations */
namespace NR{
    double I0(double x);
    double K0(double x);
    double I1(double x);
    double K1(double x);
}
/* Blair & Edwards inspired implementation */
namespace BE {
    double I0(double x);
    double K0(double x);
    double I1(double x);
    double K1(double x);
    
    double logI0(double x);
    double logK0(double x);
    double logI1(double x);
    double logK1(double x);
    
    double I0K0(double x);
    double I0K1(double x);
    double I1K0(double x);
    double I1K1(double x);
    
    double I0K0(double x1, double x2);
    double I0K1(double x1, double x2);
    double I1K0(double x1, double x2);
    double I1K1(double x1, double x2);
    
    namespace poly {
        double I0_1(double x);
        double I0_2(double x);
        
        double I1_1(double x);
        double I1_2(double x);
        
        double K0_1(double x);
        double K0_2(double x);
        double K0_3(double x);

        double K1_1(double x);
        double K1_2(double x);
        double K1_3(double x);
        
        extern const double I0_P1[];
        extern const double I0_Q1[];
        extern const double I0_P2[];
        extern const double I0_Q2[];
        
        extern const double I1_P1[];
        extern const double I1_Q1[];
        extern const double I1_P2[];
        extern const double I1_Q2[];
        
        extern const double K0_P1[];
        extern const double K0_Q1[];
        extern const double K0_P2[];
        extern const double K0_Q2[];
        extern const double K0_P3[];
        extern const double K0_Q3[];
        
        extern const double K1_P1[];
        extern const double K1_Q1[];
        extern const double K1_P2[];
        extern const double K1_Q2[];
        extern const double K1_P3[];
        extern const double K1_Q3[];
    }
}


/* Abramowitz & Stegun inspired implementation */
namespace AS {
    double I0(double x);
    double K0(double x);
    double I1(double x);
    double K1(double x);
    
    double logI0(double x);
    double logK0(double x);
    double logI1(double x);
    double logK1(double x);
    
    double I0K0(double x);
    double I0K1(double x);
    double I1K0(double x);
    double I1K1(double x);

    double I0K0(double x1, double x2);
    double I0K1(double x1, double x2);
    double I1K0(double x1, double x2);
    double I1K1(double x1, double x2);
    
    namespace poly {
        double I0_S(double x);
        double I0_L(double x);
        double I1_S(double x);
        double I1_L(double x);
        double K0_S(double x);
        double K0_L(double x);
        double K1_S(double x);
        double K1_L(double x);
        
        
        /* Asymptotic polys */
        double I0_A(double x);
        double I1_A(double x);
        double K0_A(double x);
        double K1_A(double x);

        /* Asymptotic product polys */
        double I0K0_A(double x);
        double I0K1_A(double x);
        double I1K0_A(double x);
        double I1K1_A(double x);
    }
}; /* namespace AS */

/* Toplevel: best availible options */
inline double I0(double x) { return gsl::I0(x);}
inline double K0(double x) { return gsl::K0(x);}
inline double I1(double x) { return gsl::I1(x);}
inline double K1(double x) { return gsl::K1(x);}


/* inlined methods */

inline double AS::poly::I0_S(double x)
{
    double t = square(x/3.75);
    return 1. + t*(3.5156229 + t*(3.0899424 + t*(1.2067492 + 
                t*(0.2659732 + t*(0.360768e-1 + t*0.45813e-2)))));
}

inline double AS::poly::I0_L(double x)
{
    double t = 3.75/x;
    return 0.39894228 + t*(0.1328592e-1 + t*(0.225319e-2 + t*(-0.157565e-2 + t*(0.916281e-2 + 
                        t*(-0.2057706e-1 + t*(0.2635537e-1 + t*(-0.1647633e-1 + t*0.392377e-2)))))));
}

inline double AS::poly::I1_S(double x)
{
    double t = square(x/3.75);
    return 0.5 + t*(0.87890594 + t*(0.51498869 + t*(0.15084934 + 
                 t*(0.2658733e-1 + t*(0.301532e-2 + t*0.32411e-3)))));
}

inline double AS::poly::I1_L(double x)
{
    double t = 3.75/x;
    return 0.39894228 + t*(-0.3988024e-1 + t*(-0.362018e-2 + t*(0.163801e-2 + t*(-0.1031555e-1 + 
                        t*(0.2282967e-1 + t*(-0.2895312e-1 + t*(0.1787654e-1 - t*0.420059e-2)))))));
    
}

inline double AS::poly::K0_S(double x)
{
    double t = square(x/2);
    return -0.57721566 + t*(0.42278420 + t*(0.23069756 + t*(0.3488590e-1 + 
                         t*(0.262698e-2 + t*(0.10750e-3 + t*0.74e-5)))));
}

inline double AS::poly::K0_L(double x)
{
    double t = 2./x;
    return 1.25331414 + t*(-0.7832358e-1 + t*(0.2189568e-1 + t*(-0.1062446e-1 + 
                        t*(0.587872e-2 + t*(-0.251540e-2 + t*0.53208e-3)))));
}

inline double AS::poly::K1_S(double x)
{
    double t = square(x/2.);
    return 1. + t*(0.15443144 + t*(-0.67278579 + t*(-0.18156897 + 
                t*(-0.1919402e-1 + t*(-0.110404e-2 + t*-0.4686e-4)))));
}

inline double AS::poly::K1_L(double x)
{
    double t = 2./x;
    return 1.25331414 + t*(0.23498619 + t*(-0.3655620e-1 + t*(0.1504268e-1 + 
                        t*(-0.780353e-2 + t*(0.325614e-2 + t*-0.68245e-3)))));
}


inline double AS::poly::I0_A(double x)
{
    //AS 9.7.1 u=0
    double q = 1./(8*x);
//     double q2 = square(q);
//     return 1. + q + 4.5*q2 + (75./2.)*q2*q;
    return 1. + q*(1. + q*(4.5 + q*37.5));
}

inline double AS::poly::I1_A(double x)
{
    //AS 9.7.1 u=3
    double q = 1./(8*x);
//     double q2 = square(q);
//     return 1. - (3./8.) * q - (15./2.)*q2 - (105./2.)*q2*q;
    return 1. + q*(-0.375 + q*(-7.5 + q*-52.5));
}

inline double AS::poly::K0_A(double x)
{
    //AS 9.7.2 u=0
    double q = 1./(8*x);
//     double q2 = square(q);
//     return 1. - q + 4.5*q2 - (75./2.)*q2*q;
    return 1. + q*(-1. + q*(4.5 + q*-37.5));
}

inline double AS::poly::K1_A(double x)
{
    //AS 9.7.2 u=3
    double q = 1./(8*x);
//     double q2 = square(q);
//     return 1. + (3./8.) * q - (15./2.)*q2 + (105./2.)*q2*q;
    return 1. + q*(0.375 + q*(-7.5 + q*52.5));
}


inline double AS::poly::I0K0_A(double x)
{
    double u = 1./(2*x);
    double u2 = square(u);
    return u*(1. + u2*(0.5 + u2*(27./8.)));
}

inline double AS::poly::I0K1_A(double x)
{
    double u = 1./(2*x);
    double u2 = square(u);
    return u + u2*(1. + u2*(1.5+ u2*(135./8.)));
}

inline double AS::poly::I1K0_A(double x)
{
    double u = 1./(2*x);
    double u2 = square(u);
    return u + u2*(-1. + u2*(-1.5 + u2*(135./8.)));
}

inline double AS::poly::I1K1_A(double x)
{
    double u = 1./(2*x);
    double u2 = square(u);
    return u*(1. + u2*(-1.5 + u2*(45./8.)));
}

inline double Bessel::BE::poly::I0_1(double x)
{
    double u = x*x;
    return  ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I0_P1, u, 15) / 
            ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I0_Q1, u, 6);
}

inline double Bessel::BE::poly::I0_2(double x)
{
    double u = 1./x - 1./15.;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I0_P2, u, 7) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I0_Q2, u, 8);
}

inline double Bessel::BE::poly::I1_1(double x)
{
    double u = x*x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I1_P1, u, 15) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I1_Q1, u, 6);
}

inline double Bessel::BE::poly::I1_2(double x)
{
    double u = 1./x - 1./15.;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I1_P2, u, 8) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::I1_Q2, u, 7);
}


inline double Bessel::BE::poly::K0_1(double x)
{
    double u = x*x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K0_P1, u, 6) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K0_Q1, u, 3);
}

inline double Bessel::BE::poly::K0_2(double x)
{
    double u = x*x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K0_P2, u, 5) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K0_Q2, u, 4);
}

inline double Bessel::BE::poly::K0_3(double x)
{
    double u = 1./x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K0_P3, u, 10) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K0_Q3, u, 11);
}

inline double Bessel::BE::poly::K1_1(double x)
{
    double u = x*x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K1_P1, u, 6) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K1_Q1, u, 4);
}

inline double Bessel::BE::poly::K1_2(double x)
{
    double u = x*x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K1_P2, u, 6) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K1_Q2, u, 4);
}

inline double Bessel::BE::poly::K1_3(double x)
{
    double u = 1./x;
    return ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K1_P3, u, 11) / 
           ::boost::math::tools::evaluate_polynomial(Bessel::BE::poly::K1_Q3, u, 10);
}


}; /* namespace Bessel */



#endif /* _BESSEL_H */
