/** @file Tracker.h
* @author Mark J. Olah (mjo\@cs.unm.edu)
* @date 02-2015
* @brief The class declaration and inline and templated functions for Tracker.
*
* The base class for all Tracking models
* 
* 
* Insted of templating on the FloatT type, which is problematic for inheritance hierarchies of templated
* base classes.  Instead wuse a typedef to allow configuration of use with either float/double.  Default
* is double.
* 
*/
#ifndef _TRACKER_H
#define _TRACKER_H

#include <cstdint>
#include <armadillo>
#include <map>
#include <string>
#include <list>
#include <vector>


class Tracker {
public:
    typedef double FloatT; /* Set this to control float/double settings */
    typedef int32_t IndexT; 
    typedef arma::Col<FloatT> VecT;
    typedef arma::Mat<FloatT> MatT;
    typedef arma::Col<IndexT> IVecT;
    typedef arma::Mat<IndexT> IMatT;
    typedef arma::field<IVecT> IVecFieldT;
    typedef std::vector<IndexT> IndexVectorT;
    typedef std::list<IndexT> TrackT; /**< A type for an individual track*/
    typedef std::vector<TrackT> TrackVecT;       /**< A type for a vector of tracks*/
    typedef std::map<std::string,FloatT> ParamT;  /**< A convenient form for reporting dictionaries of named FP data to matlab */
    typedef std::map<std::string,VecT> VecParamT;  /**< A convenient form for reporting dictionaries of named FP data to matlab */

    
    
    int N = 0; // Number of emitters
    int nDims = 0; //number of columns for postions 
    int nFeatures = 0;  //number of columns for features
    IVecT frameIdx; // length: N
    MatT position; // N x nDims;
    MatT SE_position; // N x nDims;
    MatT feature; // N x nFeatures;
    MatT SE_feature; // N x nFeatures;
    int firstFrame = 0; //index of first frame
    int lastFrame = 0; //index of last frame
    int nFrames = 0; //lastFrame-firstFrame+1 

    //Pre-computed on initialization
    IVecT nFrameLocs; //number of localizations for each frame, continuous indexing from firstFrame=0 to lastFrame=nFrames-1
    IVecFieldT frameLocIdx; //A field for each frame giving the indexes of localizations, continuous indexing from firstFrame=0 to lastFrame=nFrames-1

    //Computed by tracking
    TrackVecT tracks; //A vector of vectors.  Each vector reprsents a track by a sequence of localization indexes

    /**
     * param - A dictionary of floating point values to pass in.  This is a flexible interface to
     * the higher-level matlab code allowing each subclass to take in arbitrary floating point arguments!
     */
    Tracker(const VecParamT &param);
    virtual ~Tracker() {}
    virtual VecParamT getStats() const;
    virtual void initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_);
    virtual void initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_, const MatT &feature_, const MatT &SE_feature_);
    virtual void generateTracks()=0;
    
    void printTracks() const;
protected:
    static const FloatT log2pi;// = log(2*pi);
    IVecT trackAssignment; //A vector giving the track index of each loalizations
    
};

#endif /* _TRACKER_H */
