/**
 * @file Optimizer1D.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-21-2014
 * @brief The class method definitions for Optimizer1D.
 */

#ifndef _OPTIMIZER1D_H
#define _OPTIMIZER1D_H

#include <armadillo>
#include <boost/function.hpp>

template<class FloatT>
class Optimizer1D
{
public:
    typedef arma::Col<FloatT> VecT;
    typedef arma::Col<int> IVecT;
    typedef boost::function<FloatT (FloatT)> FuncT;
    
    const static FloatT phi;
    const static FloatT phi_inv;
    const static FloatT phi_conj;

    const static FloatT x_tolerance;
    const static FloatT eps;
    const static FloatT max_search_ratio;
    
    FuncT func;

    Optimizer1D(const FuncT &func, int max_eval=1000);
    void maximize(FloatT A, FloatT B, FloatT &xmax, FloatT &Fmax);
    void minimize(FloatT A, FloatT B, FloatT &xmax, FloatT &Fmax);
    void getStats(VecT &X, VecT &F) const;
    int getNFcalls() const;

    

private:
    bool maximize_mode=false; //false=minimize true=maximize
    int max_eval;
    int N; //Number of function evaluations
    VecT X; //Evaluated x location sequence
    VecT F; //Evaluated F value sequence cast as a minimization problem

    int eval(FloatT x);
    FloatT golden_step(int a, int b) const;
    FloatT parabolic_min(int a, int b, int c) const;

    IVecT bracket_min(FloatT xA, FloatT xB);
    int golden_min(const IVecT &bracket);
    int brent_min(const IVecT &bracket);
    
    FloatT max_interval_size(int a, int b, int x, FloatT step) const;
    static inline 
    void shift(int &x1, int &x2, const int x3)
    {
        x1=x2;
        x2=x3;
    }
    
    static inline 
    void shift(int &x1, int &x2, int &x3, const int x4)
    {
        x1=x2;
        x2=x3;
        x3=x4;
    }
    
};



#endif /* _OPTIMIZER1D_H */
