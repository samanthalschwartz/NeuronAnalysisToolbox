/**
 * @file DCompLib.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-21-2014
 * @brief Some fast routines for diffusion estimation computations
 */

#ifndef _DCOMPLIB_H
#define _DCOMPLIB_H

#include <armadillo>

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

/**
 * Finite difference of a vector (equivalent to Matlab diff)
 * @param[in] V - A vector of values
 * @returns diffV - A vector of length N-1, where diffV(i) = V(i+1)-V(i)
 */
template<class FloatT> 
inline
arma::Col<FloatT> vec_delta(const arma::Col<FloatT> &V)
{
    int N=V.n_elem-1;
    arma::Col<FloatT> D(N);
    for(int n=0;n<N;n++) D(n)=V(n+1)-V(n);
    return D;
}

/**
 * Implements log(prod(X)) quickly and accurately using a combination of sums of logs and logs of products.
 * This minimizes the number of logarithms necessary for the computation.  See supplementary material for
 * more information.
 * @param[in] X - vector of values
 * @returns log(prod(X))
 */
template<class FloatT>
FloatT logprod(const arma::Col<FloatT> &X);

/**
 * A templated class to represent a symmetric tri-diagonal matrix and allow efficient computation
 * of the log-determinant and a the solution of a linear equation M*x=y.
 * member variable b is the main diagonal and a is the off-diagonal.
 */
template<class FloatT>
class SymTriDiag {
public:
    typedef arma::Col<FloatT> VecT; //The type of a vector of FloatT's

    int N; //Size of the symmetric square matrix
    VecT a; //Off diagonal (index=-1,+1) length=N-1
    VecT b; //Main diagonal (index=0) length=N

    /** Make an empty matrix of size N */
    SymTriDiag(int _N) : N(_N), a(_N-1), b(_N) 
    {
        if(N<3) throw std::length_error("N<3.");
    }
    
    /** Make a new symmetric tri-diagonal matrix
     * @param[in] a - Vector length N-1 of off-diagonal elements
     * @param[in] b - Vector length N of on-diagonal elements
     */
    SymTriDiag(const VecT &_a, const VecT &_b) : N(_b.n_elem), a(_a), b(_b)
    {
        if(N<3) throw std::length_error("N<3.");
        if(static_cast<int>(a.n_elem) != N-1) throw std::length_error("Length inconsistency: len(a)!=len(b)-1");
    }

    /** Solve the linear system M*X=Y in time O(N).
     * Based on http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
     * @param[in] Y - vector size N - The RHS of the linear system to solve
     * @param[out] X - vector size N - The solution of the linear system
     */
    void solve(const VecT &Y, VecT &X) const;

    /** Return log(det(M)) in time O(N). */
    FloatT logdet() const;
};

#endif /* _DCOMPLIB_H */
