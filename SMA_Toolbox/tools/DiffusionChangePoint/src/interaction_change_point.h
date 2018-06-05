/** @file interaction_change_point.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 02-2015
 * @brief The class declaration and inline and templated functions for InteractionChangePoint
 */

#ifndef _INTERACTIONCHANGEPOINT_H
#define _INTERACTIONCHANGEPOINT_H


#include <armadillo>
#include "util.h"

    /* 
     * States: 0=free; 1=bound.
     * 
     * Units:
     * This class is unit-flexible.  We assume inputs and outputs are all given
     * using the same fundamental distance and time units, [e.g. dist=micron, time=s].
     * 
     * Types:
     *  A state sequence is seq_state0 and IVecT seq_cp - the seqeuence of changepoints
     *
     * Space/time tradeoff:  This class generally
     */
class InteractionChangePoint {
public:
    typedef double FloatT; /* Set this to control float/double settings */
    typedef int32_t IndexT; /* Set this to control int32/int64 settings */
    typedef arma::Col<IndexT> IVecT; //Type for integer vecs is based on IndexT.
    typedef arma::Col<FloatT> VecT; //Type for floating point vecs is based on FloatT.
    typedef arma::Cube<FloatT> CubeT; //Type for floating point cubes is based on FloatT.
    typedef std::vector<bool> BVecT; //Type for vector of boolean
    typedef arma::Mat<FloatT> IMatT; //Type for integer mats is based on IndexT.
    typedef arma::Mat<FloatT> MatT; //Type for floating point mats is based on FloatT.

    class Sequence {
        public:
            typedef InteractionChangePoint::FloatT FloatT;
            typedef InteractionChangePoint::IndexT IndexT;
            typedef InteractionChangePoint::IVecT IVecT;
            typedef InteractionChangePoint::BVecT BVecT;
            IndexT begin;
            IndexT end;
            bool state0;
            IVecT cps;
            Sequence() {};
            Sequence(IndexT begin, IndexT end, bool state0, const IVecT &cps);
            BVecT computeStateSeq() const;
    };

    FloatT D_A, D_B, D_AB; //Diffusion constant units units:dist^2/time
    FloatT exposureT; // frame exposure time of microscope in units:time
    FloatT rho; //SE (standard error) of bound-pair-distance in units:dist.

    InteractionChangePoint(FloatT D_A, FloatT D_B, FloatT D_AB, FloatT exposureT, FloatT rho);

    void freeSegmentLLH(const MatT &pair, MatT &LLH) const;
    void boundSegmentLLH(const MatT &pair, MatT &LLH) const;

    FloatT sequenceLLH(const MatT &pair, const Sequence &seq) const;
    FloatT sequenceLLH_debug(const MatT &pair, const Sequence &seq, CubeT &LLH_components) const;

    
    VecT compute1CPLLH(const MatT &pair, const Sequence &seq) const;
    MatT compute2CPLLH(const MatT &pair, const Sequence &seq) const;
    
    
    static FloatT segmentProductLLH(const MatT &freeSegLLH, const MatT &boundSegLLH, const Sequence &seq);
    static FloatT segmentProductMLE(const MatT &freeLLH, const MatT &boundLLH, double beta, Sequence &seq);
    static void computeProductMaxLikelihood(const MatT &segLLH0, const MatT &segLLH1, double beta,
                                     VecT &maxLLH0, VecT &maxLLH1, IVecT &mleCP0, IVecT &mleCP1);

    Sequence makeSequence(IndexT begin, IndexT end, bool state0, const IVecT &cps) const;
    Sequence makeSequence(IndexT begin, IndexT end, bool state0) const;

private:
    FloatT min_variance = 1E-8;  //A limit on the minimum variance for the LLH computations
    static const FloatT log2pi;

    VecT free1DSegmentLLH(const VecT &obs, const VecT &vM, const VecT &vD, int startIdx) const;
    VecT bound1DSegmentLLH(const VecT &obs_a, const VecT &obs_b, const VecT &vMa, const VecT &vMb, const VecT &vD, int startIdx) const;
    FloatT full1DLLH(const VecT &T, const VecT &obs_a, const VecT &obs_b, const VecT &var_a, const VecT &var_b, const Sequence &seq) const;
    FloatT full1DLLH_debug(const VecT &T, const VecT &obs_a, const VecT &obs_b, const VecT &var_a, const VecT &var_b,
                           const Sequence &seq, IndexT &Ncomponents, MatT &LLH_components) const;

                           VecT computeDVariance(FloatT D, const VecT &T) const;
    VecT computeMVariance(FloatT D, const VecT &SE) const;
    FloatT boundVariance(FloatT var) const;
    FloatT normalDistLLH(FloatT a, FloatT b, FloatT var) const; //Quick method for computing the LLH of a single normal dist N(a,b,var);
    FloatT normalDistLLH_debug(FloatT a, FloatT b, FloatT v, IndexT &Ncomponents, MatT &LLH_components) const;

};

InteractionChangePoint::Sequence::Sequence(IndexT begin, IndexT end, bool state0, const IVecT &cps)
    : begin(begin), end(end), state0(state0), cps(cps)
{
    //TODO:check sequence
}

InteractionChangePoint::Sequence::BVecT
InteractionChangePoint::Sequence::computeStateSeq() const
{
    IndexT nstates = end-begin+1;
    BVecT state_seq(nstates);
    IndexT next_cp_idx = 0;
    IndexT next_cp = cps.is_empty() ? end : cps(next_cp_idx);
    bool state = state0;
    for(IndexT i=begin; i<=end; i++){
        if(i==next_cp && next_cp<end){
            state = !state;
            next_cp_idx++;
            next_cp = next_cp_idx == static_cast<IndexT>(cps.n_elem) ? end : cps(next_cp_idx);
        }
        state_seq[i-begin] = state;
    }
    return state_seq;
}

inline
InteractionChangePoint::FloatT
InteractionChangePoint::boundVariance(FloatT var) const
{
    if(fabs(var)<min_variance) var = sgn(var)<0 ? -min_variance : min_variance;
    return var;
}

inline
InteractionChangePoint::FloatT
InteractionChangePoint::normalDistLLH(FloatT a, FloatT b, FloatT var) const
{
    return log2pi + log(var) + square(a-b)/var;
}

inline
InteractionChangePoint::Sequence
InteractionChangePoint::makeSequence(IndexT begin, IndexT end, bool state0, const IVecT &cps) const
{
    return Sequence(begin,end,state0,cps);
}

inline
InteractionChangePoint::Sequence
InteractionChangePoint::makeSequence(IndexT begin, IndexT end, bool state0) const
{
    IVecT cps;
    return Sequence(begin,end,state0,cps);
}

#endif/* _INTERACTIONCHANGEPOINT_H */
