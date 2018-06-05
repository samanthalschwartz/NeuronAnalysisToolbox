/** @file util.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2015
 * @brief Utility helper functions
 */

#ifndef _UTIL_H
#define _UTIL_H

#include <armadillo>
#include <limits>
#include <cmath>

// /** Fast computation of a logarithm of a product */
// template<class FloatT>
// FloatT logprod(const arma::Col<FloatT> &X)
// {
//     static const FloatT max=std::numeric_limits<FloatT>::max();
//     static const FloatT min=std::numeric_limits<FloatT>::min();
//     FloatT prod=1;
//     FloatT Lprod=0;
//     for(unsigned n=0; n<X.n_elem; n++) {
//         FloatT x=fabs(X(n));
//         FloatT temp_prod=prod*x;
//         if(temp_prod<min || temp_prod>max) {
//             Lprod+=log(prod);
//             prod=x;
//         } else {
//             prod=temp_prod;
//         }
//     }
//     Lprod+=log(prod);
//     return Lprod;
// }


/** Fast computation of sign [-1,0,+1] for any type */
template <class T>
inline
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

/** Fast computation of square of a number */
template <class T>
inline
T square(const T &val)
{
    return val*val;
}

#endif /* _UTIL_H */
