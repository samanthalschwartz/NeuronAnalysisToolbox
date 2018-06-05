/** @file InteractionChangePoint_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 04-2015
 * @brief The class declaration and inline and templated functions for InteractionChangePoint_Iface.
 */

#ifndef _INTERACTIONCHANGEPOINT_IFACE_H
#define _INTERACTIONCHANGEPOINT_IFACE_H

#include "Mex_Iface.h"
#include "interaction_change_point.h"

class InteractionChangePoint_Iface : public Mex_Iface
{
public:
    InteractionChangePoint_Iface();
private:
    InteractionChangePoint *obj;
    //Abstract member functions inherited from Mex_Iface
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);
    //Exposed method calls
    void objFreeSegmentLLH();
    void objBoundSegmentLLH();
    void objSegmentProductMLE();
    void objSequenceLLH();
    void objSequenceLLH_debug();
    void objCompute1CPLLH();
    void objCompute2CPLLH();
};


InteractionChangePoint_Iface::InteractionChangePoint_Iface()
    : Mex_Iface("InteractionChangePoint_Iface")
{
    methodmap["freeSegmentLLH"] = boost::bind(&InteractionChangePoint_Iface::objFreeSegmentLLH, this);
    methodmap["boundSegmentLLH"] = boost::bind(&InteractionChangePoint_Iface::objBoundSegmentLLH, this);
    methodmap["segmentProductMLE"] = boost::bind(&InteractionChangePoint_Iface::objSegmentProductMLE, this);
    methodmap["sequenceLLH"] = boost::bind(&InteractionChangePoint_Iface::objSequenceLLH, this);
    methodmap["sequenceLLH_debug"] = boost::bind(&InteractionChangePoint_Iface::objSequenceLLH_debug, this);
    methodmap["compute1CPLLH"] = boost::bind(&InteractionChangePoint_Iface::objCompute1CPLLH, this);
    methodmap["compute2CPLLH"] = boost::bind(&InteractionChangePoint_Iface::objCompute2CPLLH, this);
}

void InteractionChangePoint_Iface::objConstruct()
{
    checkNumArgs(1,5);
    double D_A = getDouble();
    double D_B = getDouble();
    double D_AB = getDouble();
    double exposureT = getDouble();
    double rho = getDouble();
    auto *obj = new InteractionChangePoint(D_A, D_B, D_AB, exposureT, rho);
    outputMXArray(Handle<InteractionChangePoint>::makeHandle(obj));
}

void InteractionChangePoint_Iface::objDestroy()
{
    checkNumArgs(0,1);
    Handle<InteractionChangePoint>::destroyObject(rhs[0]);
}

void InteractionChangePoint_Iface::getObjectFromHandle(const mxArray *mxhandle)
{
    obj = Handle<InteractionChangePoint>::getObject(mxhandle);
}

void InteractionChangePoint_Iface::objFreeSegmentLLH()
{
    // [in] pair: double mat with nObs rows and 18 columns as defined by the PairAnalysis matlab class
    // [out] LLH: log-likelihood matrix to be filled in with the liklihoods for [j,i] where i<j corresponding to observations i:j comming from
    //               particles in the free state.  The size should be [nObs,nObs].

    checkNumArgs(1,1);
    auto pair=getMat<InteractionChangePoint::FloatT>();
    auto LLH = makeMat<InteractionChangePoint::FloatT>(pair.n_rows,pair.n_rows); // create output LLH matrix
    obj->freeSegmentLLH(pair, LLH);
}

void InteractionChangePoint_Iface::objBoundSegmentLLH()
{
    // [in] pair: double mat with nObs rows and 18 columns as defined by the PairAnalysis matlab class
    // [out] LLH: log-likelihood matrix to be filled in with the liklihoods for [j,i] where i<j corresponding to observations i:j comming from
    //               particles in the bound state.  The size should be [nObs,nObs].

    checkNumArgs(1,1);
    auto pair=getMat<InteractionChangePoint::FloatT>();
    auto LLH = makeMat<InteractionChangePoint::FloatT>(pair.n_rows,pair.n_rows); // create output LLH matrix
    obj->boundSegmentLLH(pair, LLH);
}

void InteractionChangePoint_Iface::objSegmentProductMLE()
{
    // This is the liklihood for the rigid attachment model.  (i.e. there is a fixed displacment vector maintained in the interaction.)
    // [in] freeLLH
    // [in] boundLLH
    // [in] beta - cp penalty
    // [out] seq_state0
    // [out] seq_cps (0-based) 1<=cps<=N-2;
    // [out] seq_llh
    checkNumArgs(3,3);
    auto freeLLH = getMat<InteractionChangePoint::FloatT>();
    auto boundLLH = getMat<InteractionChangePoint::FloatT>();
    double beta = getDouble();
    InteractionChangePoint::Sequence seq;
    InteractionChangePoint::FloatT seq_llh = obj->segmentProductMLE(freeLLH, boundLLH, beta, seq);
    outputBool(seq.state0);
    outputVec(seq.cps);
    outputDouble(seq_llh);
}

void InteractionChangePoint_Iface::objCompute1CPLLH()
{
    // This is the liklihood for the rigid attachment model.  (i.e. there is a fixed displacment vector maintained in the interaction.)
    // [in] pairMat
    // [in] seq_beg
    // [in] seq_end
    // [in] seq_state0
    // [out] seq_llh
    checkNumArgs(1,4);
    auto pair = getMat<InteractionChangePoint::FloatT>();
    int seq_beg = getInt();
    int seq_end = getInt();
    bool seq_state0 = getBool();
    outputVec(obj->compute1CPLLH(pair, obj->makeSequence(seq_beg,seq_end,seq_state0)));
}

void InteractionChangePoint_Iface::objCompute2CPLLH()
{
    // This is the liklihood for the rigid attachment model.  (i.e. there is a fixed displacment vector maintained in the interaction.)
    // [in] pairMat
    // [in] seq_beg
    // [in] seq_end
    // [in] seq_state0
    // [out] seq_llh
    checkNumArgs(1,4);
    auto pair = getMat<InteractionChangePoint::FloatT>();
    int seq_beg = getInt();
    int seq_end = getInt();
    bool seq_state0 = getBool();
    outputMat(obj->compute2CPLLH(pair, obj->makeSequence(seq_beg,seq_end,seq_state0)));
}

void InteractionChangePoint_Iface::objSequenceLLH()
{
    // [in] pair: double mat with nObs rows and 18 columns as defined by the PairAnalysis matlab class
    // [in] seq_beg: integer 0<=seq_beg<N-2
    // [in] seq_end: integer 1<=seq_end<N-1
    // [in] seq_state0: bool
    // [in] seq_cps: integer vector all(1<=cps<=N-2)
    // [out] LLH: log-likelihood matrix to be filled in with the liklihoods for [j,i] where i<j corresponding to observations i:j comming from
    //               particles in the bound state.  The size should be [nObs,nObs].

    checkNumArgs(1,5);
    auto pair=getMat<InteractionChangePoint::FloatT>();
    int seq_beg = getInt();
    int seq_end = getInt();
    bool seq_state0 = getBool();
    auto seq_cps = getIVec();
    auto seq = obj->makeSequence(seq_beg, seq_end, seq_state0, seq_cps);
    auto LLH = obj->sequenceLLH(pair, seq);
    outputDouble(LLH);
}

void InteractionChangePoint_Iface::objSequenceLLH_debug()
{
    // [in] pair: double mat with nObs rows and 18 columns as defined by the PairAnalysis matlab class
    // [in] seq_beg: integer 0<=seq_beg<N-2
    // [in] seq_end: integer 1<=seq_end<N-1
    // [in] seq_state0: bool
    // [in] seq_cps: integer vector all(1<=cps<=N-2)
    // [out] LLH: log-likelihood matrix to be filled in with the liklihoods for [j,i] where i<j corresponding to observations i:j comming from
    //               particles in the bound state.  The size should be [nObs,nObs].

    checkNumArgs(2,5);
    auto pair=getMat<InteractionChangePoint::FloatT>();
    int seq_beg = getInt();
    int seq_end = getInt();
    bool seq_state0 = getBool();
    auto seq_cps = getIVec();
    auto seq = obj->makeSequence(seq_beg, seq_end, seq_state0, seq_cps);
    InteractionChangePoint::CubeT LLH_components;
    auto LLH = obj->sequenceLLH_debug(pair, seq, LLH_components);
    outputDouble(LLH);
    outputStack(LLH_components);
}


#endif /* _INTERACTIONCHANGEPOINT_IFACE_H */
