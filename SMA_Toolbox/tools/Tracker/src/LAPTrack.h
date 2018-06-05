/** @file LAPTrack.h
* @author Mark J. Olah (mjo\@cs.unm.edu)
* @date 02-2015
* @brief The class declaration and inline and templated functions for LAPTrack.
*
* A simple  LAP/Jaquman based tracker
*/

#ifndef _LAPTRACK_H
#define _LAPTRACK_H

#include "Tracker.h"

class LAPTrack : public Tracker {
public:
    typedef arma::SpMat<FloatT> SpMatT;
    typedef arma::Col<arma::uword> UVecT;

    typedef arma::umat UMatT;
    FloatT D; //  D - um^2/s
    FloatT kon;//  kon  - s^-1
    FloatT koff;//  koff - s^-1
    FloatT rho;
    VecT featureVar; //(nFeatures,1) The sigma^2 of error allowed in each feature dimension.  Should have same units as feature.
    FloatT maxSpeed = 0;  //Maximum speed
    FloatT maxPositionDisplacementSigma = 5.0; //Maximum standard deviations out to propose a connection
    VecT maxFeatureDisplacementSigma;  //Maximum standard deviations out to propose a connection
    IndexT maxGapCloseFrames = 20;
    IndexT minGapCloseTrackLength = 1;
    IndexT minFinalTrackLength = 1;
    
    
    const FloatT cost_epsilon = std::numeric_limits<FloatT>::epsilon();

    LAPTrack(const VecParamT &param);
    VecParamT getStats() const;
    void initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_);
    void initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_, const MatT &feature_, const MatT &SE_feature_);
    void linkF2F();
    void closeGaps();
    SpMatT computeF2FCostMat(int curFrame, int nextFrame) const;
    void debugF2F(int frameIdx, IVecT &cur_locs, IVecT &next_locs, SpMatT &cost, IMatT &connections, VecT &conn_costs) const;
    void debugCloseGaps(SpMatT &cost, IMatT &connections, VecT &conn_costs) const;
    
    SpMatT computeGapCloseMatrix() const;
    void generateTracks();
    void checkFrameIdxs();
protected:
    FloatT minCost = 1e-6; // The minimum cost to put in the matrix.  Should this be bigger than machine eps?
    FloatT log1mkoff; //log(1-koff);
    FloatT log1mkon; //log(1-kon);
    FloatT logrho; //log(rho);
    FloatT logkon; //log(kon);
    FloatT logkoff; //log(koff);

    enum StateT {UNTRACKED, F2F_LINKED, GAPS_CLOSED};
    StateT state;
    //These member variables make it easier for gap closing to get the information on the track beginning/ends
    //We can assemble the gap closing index without searching for tracks that are born at a particular frame.
    IndexVectorT birthFrameIdx; //frame index for the beginning of each track.
    IVecT frameBirthStartIdx; //for each frameIdx list the first track to start at that frame Idx or later.
    
};

#endif /* _LAPTRACK_H */
