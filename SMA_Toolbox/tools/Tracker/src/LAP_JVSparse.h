/** @file LAP_JVSparse.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 05-2015
 * @brief The class declaration for the LAP Jonker Volgenant algorithm
 *
 * This is a modern dense/sparse C++ implementation of Jonker Volgenant algoirthm using armadillo
 * and presenting C++ and Matlab interface.
 * 
 * Adapted from text of Jonker and Volgenant. Computing 38, 324-340 (1986)
 * 
 */
#ifndef _LAP_JVSPARSE_H
#define _LAP_JVSPARSE_H

#include <armadillo>
#include <vector>

template<class FloatT>
class LAP_JVSparse {
    typedef int32_t IndexT; //The type for the indexes
    typedef arma::SpMat<FloatT> SpMatT;
    typedef arma::Col<FloatT> VecT;
    typedef arma::Col<IndexT> IVecT;
    typedef arma::Mat<IndexT> IMatT;

public:
    static IVecT solve(const SpMatT &C);
    static void solveLAP_orig(const SpMatT &C, IVecT &x, IVecT &y, VecT &u, VecT &v);
    static VecT computeCost(const SpMatT &C, const IVecT &row_sol);

    static bool checkCosts(const SpMatT &C);
    static bool checkSolution(const SpMatT &C,const IVecT &x, const IVecT &y, const VecT &u, const VecT &v);

private:
    /* The original sparse lapjv code which is outdated and should be updated. */
    static void lap_orig(IndexT n, const FloatT C_vals[], const IndexT C_cols[], const IndexT C_row_ptrs[],
                           IndexT x[], IndexT y[], FloatT u[], FloatT v[]);
    
};


#endif /* _LAP_JVSPARSE_H */
