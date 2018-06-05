/** @file Tracker.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 04-2015
 *  @brief The member definitions for Tracker
 */
#include <cassert>
#include <cmath>
#include "Tracker.h"

const Tracker::FloatT Tracker::log2pi = log(2*arma::Datum<Tracker::FloatT>::pi);

Tracker::Tracker(const VecParamT &param)
{
}

Tracker::VecParamT Tracker::getStats() const
{
    VecParamT stats;
    stats["nLocalizations"] = N;
    stats["nDims"] = nDims;
    stats["nFeatures"] = nFeatures;
    stats["firstFrame"] =firstFrame;
    stats["lastFrame"] = lastFrame;
    stats["nFrames"] = nFrames;
    stats["nTracks"] = tracks.size();
    stats["nLocalizationsAssigned"] = arma::sum(trackAssignment>=0);
    return stats;
}

void Tracker::initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_)
{
    MatT feature_, SE_feature_;
    initializeTracks(frameIdx_, position_, SE_position_,feature_, SE_feature_); //Call with empty features
}

void Tracker::initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_, const MatT &feature_, const MatT &SE_feature_)
{
    assert(frameIdx_.n_elem == position_.n_rows);
    assert(frameIdx_.n_elem == SE_position_.n_rows);
    assert(position_.n_cols == SE_position_.n_cols);
    assert(feature_.n_cols  == SE_feature_.n_cols);
    if(!feature.is_empty()) assert(frameIdx.n_elem == feature_.n_rows);
    N = static_cast<int>(frameIdx_.n_elem);
    nDims = static_cast<int>(position_.n_cols)/2;
    nFeatures = static_cast<int>(feature_.n_cols)/2;
    frameIdx = frameIdx_;
    position = position_;
    SE_position = SE_position_;
    feature = feature_;
    SE_feature = SE_feature_;

    //Clear track data structures
    tracks.clear();
    tracks.reserve(static_cast<int>(ceil(sqrt(N))));
    trackAssignment.set_size(N);
    trackAssignment.fill(-1);

    //Initialize number of frames and range
    arma::uvec sFrameIdx = arma::stable_sort_index(frameIdx);//Ensure sort is stable.
    firstFrame = frameIdx(sFrameIdx(0));
    lastFrame = frameIdx(sFrameIdx(N-1));
    nFrames = lastFrame-firstFrame+1;
    //Initialize frameLocIdx - the list of localizations for each frame
    IndexVectorT buf;//buffer to store the current frames indexs
    buf.reserve(std::max(N,10*(N/nFrames)));//Pre-allocate approximate ammount of needed space for efficiency
    nFrameLocs.set_size(nFrames);
    frameLocIdx.set_size(nFrames);
    int cur_frame=firstFrame; // The current frame we are processing
    for(int n=0; n<N; n++){
        int next_loc_idx = sFrameIdx(n);
        int next_loc_frame = frameIdx(next_loc_idx);
        if (next_loc_frame>cur_frame) { //Next loc is from a future frame.
            frameLocIdx(cur_frame-firstFrame) = IVecT(buf);
            buf.clear(); //Reset for new cur_frame;
            while(++cur_frame<next_loc_frame) frameLocIdx(cur_frame-firstFrame).reset(); //zero-out blank frames
        }
        buf.push_back(next_loc_idx);//Add to loc to cur_frame
    }
    frameLocIdx(cur_frame-firstFrame) = IVecT(buf);
    assert(cur_frame==lastFrame);
    for(int n=0; n<nFrames; n++) nFrameLocs(n)= frameLocIdx(n).n_elem;
    
    int sum=0;
    for(int n=0; n<nFrames; n++) sum+=frameLocIdx(n).n_elem;
//     std::cout<<"frameLocIdx(0): "<<frameLocIdx(0).t()<<std::endl;
//     std::cout<<"frameLocIdx(1): "<<frameLocIdx(1).t()<<std::endl;
//     std::cout<<"nFrameLocs: "<<nFrameLocs.t()<<std::endl;
//     
//     std::cout<<"Sum: "<<sum<<std::endl;
    assert(sum==N);
}

void Tracker::printTracks() const
{
    int nTracks = static_cast<int>(tracks.size());
    std::cout<<"Number of tracks: "<<nTracks<<"\n";
    for(int n=0; n<nTracks; n++){
        int start = frameIdx(tracks[n].front());
        int end = frameIdx(tracks[n].back());
        int len = static_cast<int>(tracks[n].size());
        std::cout<<"Track["<<n<<"]: StartFrame: "<<start<<" EndFrame: "<<end<<" #Locs:"<<len<<"\n";
        std::cout<<"  Locs: ";
        for(auto k=tracks[n].cbegin(); k!=tracks[n].cend(); ++k)
            std::cout<<*k<<", ";
        std::cout<<std::endl;
    }
}

