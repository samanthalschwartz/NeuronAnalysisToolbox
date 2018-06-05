/** @file LAPTrack.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 04-2015
 *  @brief The member definitions for LAPTrack
 */
#include <cassert>
#include "LAPTrack.h"
#include "LAP_JVSparse.h"

LAPTrack::LAPTrack(const VecParamT &param) : Tracker(param)
{
//     for(auto &stat: param)
//         std::cout<<stat.first<<":"<<stat.second<<std::endl;
    //Read in parameters
    if (param.find("D") != param.end())
        D = static_cast<FloatT>(param.at("D")(0));
    if (param.find("kon") != param.end())
        kon = static_cast<FloatT>(param.at("kon")(0));
    if (param.find("koff") != param.end())
        koff = static_cast<FloatT>(param.at("koff")(0));
    if (param.find("rho") != param.end())
        rho = static_cast<FloatT>(param.at("rho")(0));
    if (param.find("maxSpeed") != param.end())
        maxSpeed = static_cast<FloatT>(param.at("maxSpeed")(0));
    if (param.find("maxPositionDisplacementSigma") != param.end())
        maxPositionDisplacementSigma = static_cast<FloatT>(param.at("maxPositionDisplacementSigma")(0));
    if (param.find("maxFeatureDisplacementSigma") != param.end())
        maxFeatureDisplacementSigma =  static_cast<FloatT>(param.at("maxFeatureDisplacementSigma")(0));
    if (param.find("maxGapCloseFrames") != param.end())
        maxGapCloseFrames =  static_cast<int>(param.at("maxGapCloseFrames")(0));
    if (param.find("minGapCloseTrackLength") != param.end())
        minGapCloseTrackLength =  static_cast<int>(param.at("minGapCloseTrackLength")(0));
    if (param.find("minFinalTrackLength") != param.end())
        minFinalTrackLength =  static_cast<int>(param.at("minFinalTrackLength")(0));
    if (param.find("featureVar") != param.end())
        featureVar = param.at("featureVar");
    //Pre-compute logarithms of commonly used values
    logkon = log(kon);
    log1mkon = log(1-kon);
    logkoff = log(koff);
    log1mkoff = log(1-koff);
    logrho = log(rho);
}

LAPTrack::VecParamT LAPTrack::getStats() const
{
    auto stats = Tracker::getStats();
    stats["D"] = D;
    stats["kon"] =kon;
    stats["koff"] = koff;
    stats["rho"] = rho;
    stats["maxSpeed"] = maxSpeed;
    stats["maxPositionDisplacementSigma"] = maxPositionDisplacementSigma;
    stats["maxFeatureDisplacementSigma"] = maxFeatureDisplacementSigma;
    stats["maxGapCloseFrames"] = maxGapCloseFrames;
    stats["minGapCloseTrackLength"] = minGapCloseTrackLength;
    stats["minFinalTrackLength"] = minFinalTrackLength;
    stats["featureVar"] = featureVar;
    return stats;
}

void LAPTrack::initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_)
{
    MatT feature_, SE_feature_;
    initializeTracks(frameIdx_, position_, SE_position_,feature_, SE_feature_);
}

void LAPTrack::initializeTracks(const IVecT &frameIdx_, const MatT &position_, const MatT &SE_position_, const MatT &feature_, const MatT &SE_feature_)
{
    Tracker::initializeTracks(frameIdx_, position_, SE_position_, feature_,SE_feature_);
    frameBirthStartIdx.clear();
    birthFrameIdx.clear();
    state = UNTRACKED;
}

void LAPTrack::generateTracks()
{
    //Do whatever is still needed to produce the tracks
    switch(state){
        case UNTRACKED:
            linkF2F();
        case F2F_LINKED:
            closeGaps();
        case GAPS_CLOSED:
            break;
    }
}

void LAPTrack::debugF2F(int curFrame, IVecT &cur_locs, IVecT &next_locs, SpMatT &cost, IMatT &connections, VecT &conn_costs) const
{
    assert(curFrame!=lastFrame);
    int nextFrame = curFrame+1;
    while(frameLocIdx(nextFrame-firstFrame).is_empty()) nextFrame++; //Find next frame with localizations
    cur_locs = frameLocIdx(curFrame-firstFrame);
    next_locs = frameLocIdx(nextFrame-firstFrame);
    cost = computeF2FCostMat(curFrame, nextFrame);
    IVecT frame_assignment = LAP_JVSparse<FloatT>::solve(cost);
    int nCurLocs = nFrameLocs(curFrame-firstFrame);
    int nNextLocs = nFrameLocs(nextFrame-firstFrame);
    //Comput nCons - number of connections.  This includes track->track, births, and deaths, but not phantom connections.
    int nCons = nCurLocs;
    for(int n=nCurLocs; n<nCurLocs+nNextLocs; n++) {
        if(frame_assignment(n)<nNextLocs) {
            nCons++;
        }//otherwise is a phantom connection so don't report it.
    }
    connections.set_size(nCons,2);
    int conIdx=0;
    for(int n=0;n<nCurLocs+nNextLocs;n++){
        if (n>=nCurLocs) {
            if (frame_assignment(n)>=nNextLocs) continue; //phantom
            connections(conIdx,0) = -1; //birth
        } else {
            connections(conIdx,0) = cur_locs(n); //connection
        }
        connections(conIdx,1) = (frame_assignment(n)>=nNextLocs) ? -1 : next_locs(frame_assignment(n)); //deaths
        conIdx++;
    }
    conn_costs = LAP_JVSparse<FloatT>::computeCost(cost, frame_assignment);
    conn_costs = conn_costs.elem(arma::find(conn_costs>cost_epsilon));
}

void LAPTrack::linkF2F()
{
    assert(state==UNTRACKED);
    //firstFrame and lastFrame are guaranteed to have the localizations others may not
    //we choose to connect with the next non-empty frame
    IndexT curFrame = firstFrame;
    //Initialize first frame of tracks
    IVecT &initLocs = frameLocIdx(0);
    frameBirthStartIdx.set_size(nFrames);
    for(int i=0; i< nFrameLocs(0); i++){
        int locIdx = initLocs(i);
        tracks.push_back(TrackT());
        tracks[i].push_back(locIdx);
        trackAssignment[locIdx]=i;
        frameBirthStartIdx(curFrame-firstFrame) = 0; // record births
        birthFrameIdx.push_back(curFrame); // record birth frame time
    }

    while(curFrame < lastFrame){  //When curFrame==lastFrame we have linked all frames
        IndexT nextFrame = curFrame+1;
//         std::cout<<"------------F2F------------"<<"\n";
        while(frameLocIdx(nextFrame-firstFrame).is_empty()){
            frameBirthStartIdx(nextFrame-firstFrame) = static_cast<int>(tracks.size());//Record absence of new births for frame nextFrame
            nextFrame++;
        }
//         std::cout<<"curFrame:"<<curFrame<<" nextFrame:"<<nextFrame<<"\n";
//         std::cout<<"trackAssignment: "<<trackAssignment.t()<<"\n";
//         std::cout<<"NTracks:"<<tracks.size()<<"\n";
        int nCur = nFrameLocs(curFrame-firstFrame);
        assert(nCur>0);
        int nNext=nFrameLocs(nextFrame-firstFrame);
//         std::cout<<"Ncur:"<<nCur<<" Nnext:"<<nNext<<"\n";
        
        SpMatT cost = computeF2FCostMat(curFrame, nextFrame); //Make the cost sparse matrix
        arma::mat dC(cost);
//         std::cout<<"Cost: ("<<cost.n_rows<<","<<cost.n_cols<<"):\n"<<dC<<"\n";
        IVecT frame_assignment = LAP_JVSparse<FloatT>::solve(cost); //Solve for the assignments.
//         std::cout<<"frameAssignment: "<<frame_assignment.t()<<"\n";
        IVecT &curFrameIdxs = frameLocIdx(curFrame-firstFrame);
        IVecT &nextFrameIdxs = frameLocIdx(nextFrame-firstFrame);
//         std::cout<<"curFrameIdxs: "<<curFrameIdxs.t()<<"\n";
//         std::cout<<"nextFrameIdxs: "<<nextFrameIdxs.t()<<"\n";

        
        
        for(int i=0; i<nCur; i++){
            //process frame_assignment for each of the current frame localizations
            int asgn = frame_assignment(i); //In terms of the possible next frames connections or deaths.
            int cur_id = curFrameIdxs(i);
            int track_id = trackAssignment(cur_id);
//              {//death
//                 std::cout<<"Death: Loc:"<<cur_id<<"Track:"<<track_id<<"\n";
//                 deathLocIdx[track_id] = cur_id; // record death index
//                 deathFrameIdx[track_id] = curFrame; // record deth frame
            if (asgn < nNext) { //connection - extend track corresponding to current frame localization
                assert(track_id>=0);
                int next_loc_idx = nextFrameIdxs(asgn);
                assert(trackAssignment(next_loc_idx)==-1);  //unassigned
                trackAssignment(next_loc_idx) = track_id;
                tracks[track_id].push_back(next_loc_idx);
//                 std::cout<<"Connect: Track:"<<track_id<<" "<<cur_id<<"->"<<next_loc_idx<<"\n";
            }
        }
        //Process births
        frameBirthStartIdx(nextFrame-firstFrame) = tracks.size(); //Births for next frame start with next track added.
//         std::cout<<"Frame Birth Start Idx: "<<frameBirthStartIdx.t()<<"\n";
        for(int i=nCur; i<nCur+nNext; i++) {
            int birth_id = i-nCur; //index into nextFrameIdx of this localization
            int asgn = frame_assignment(i);
            if (asgn < nNext) { //birth - make a track of size 1.
                int track_id = tracks.size(); //new track_id
//                 std::cout<<"Birthing NewTrackID: "<<track_id<<std::endl;
                
//                 assert(track_id>=1);
                int birth_loc_idx = nextFrameIdxs(birth_id); //Actual localization index of the birth localization
                assert(trackAssignment(birth_loc_idx)==-1);  //unassigned
                trackAssignment(birth_loc_idx) = track_id;
                tracks.push_back(TrackT());
                tracks[track_id].push_back(birth_loc_idx);
                birthFrameIdx.push_back(nextFrame); // record birth frame time as happening in next frame
//                 std::cout<<"birthId:"<<birth_id<<" birthLocIdx:"<<birth_loc_idx<<"\n";
//                 std::cout<<"trackAssignment:"<<trackAssignment.t()<<"\n";
//                 std::cout<<"recorded birthFrameIdx: "<<birthFrameIdx.back()<<"\n";
            }
        }
        curFrame=nextFrame;
    }
    //Finalize deaths of tracks still going in last frame
//     for(unsigned i=0; i<deathLocIdx.size();i++){
//         if(deathLocIdx[i]==-1) {
//             int lastIdx = tracks[i].back();
//             std::cout<<"Finalize Tracks: TrackID: "<<i<<" lastIdx:"<<lastIdx<<" frameIdx:"<<frameIdx(lastIdx)<<" deathLocIdx:"<<deathLocIdx[i]<<"\n";
//             assert(frameIdx(lastIdx)==lastFrame);
//         }
//     }
//     std::cout<<"NTracks: "<<tracks.size()<<"\n";
//     std::cout<<"TrackAssignment: "<<trackAssignment.t()<<"\n";
//     std::cout<<"frameBirthStartIdx: "<<IVecT(frameBirthStartIdx).t()<<"\n";
//     std::cout<<"BirthFrameIdx: "<<IVecT(birthFrameIdx).t()<<"\n";
    state = F2F_LINKED;
}

LAPTrack::SpMatT
LAPTrack::computeF2FCostMat(int curFrame, int nextFrame) const
{
    int nCur = nFrameLocs(curFrame-firstFrame);
    int nNext = nFrameLocs(nextFrame-firstFrame);
    int nTot = nCur+nNext;
    //These will be our sparse matrix format vectors
    std::vector<arma::uword> row_index;
    std::vector<arma::uword> col_index;
    std::vector<FloatT> values;
    int reserve_size = nTot+2*std::min(nCur*nNext, std::max(nCur,nNext)*10); //Guesstimate amount of entries used
    row_index.reserve(reserve_size);
    col_index.reserve(reserve_size);
    values.reserve(reserve_size);

    int deltaT = nextFrame - curFrame; //The number of frames spanned in the link
    FloatT DdT = 2*D*deltaT;
    FloatT position_gaussian_exponent_cuttoff = (maxPositionDisplacementSigma*maxPositionDisplacementSigma)/2.; //Only allow connections within 5 sigma
    VecT feature_gaussian_exponent_cuttoff = (maxFeatureDisplacementSigma%maxFeatureDisplacementSigma)/2.; //Only allow connections within 5 sigma
    FloatT norm_const = (nDims+nFeatures)*log2pi; //Pre-compute this

    //Fill in connection costs
    const IVecT &curFrameLocs = frameLocIdx(curFrame-firstFrame);
    const IVecT &nextFrameLocs = frameLocIdx(nextFrame-firstFrame);
//     std::cout<<"nCur:"<<nCur<<" nNext:"<<nNext<<"\n";
//     std::cout<<"nFrameLocs:"<<nFrameLocs.t()<<"\n";
    //connecting locI in current frame to locJ in next frame
    for(int j=0; j<nNext; j++) {
        int next_idx = nextFrameLocs(j);
        for(int i=0; i<nCur; i++){
            int cur_idx = curFrameLocs(i);
//             std::cout<<"i:"<<i<<" j:"<<j<<" curIdx:"<<cur_idx<<" nextIdx:"<<next_idx<<" DdT:"<<DdT<<"\n";
            FloatT C=0;
            bool feasible = true;
            FloatT total_dist_sq=0;
            for(int d=0; d<nDims; d++){
                FloatT dist_var = DdT + SE_position(cur_idx,d) + SE_position(next_idx,d);
                FloatT dist = position(cur_idx,d) - position(next_idx,d);
                FloatT dist_sq = dist*dist;
                total_dist_sq += dist_sq;
                FloatT cost_exponent = dist_sq/dist_var;
//                 std::cout<<"Dim:"<<d<<" dist:"<<dist<<" dist_var:"<<dist_var<<" costExp:"<<cost_exponent<<" ExpCuttoff:"<<position_gaussian_exponent_cuttoff<<"\n";
//                 std::cout<<"SE_position(cur_idx,d):"<<SE_position(cur_idx,d)<<"SE_position(next_idx,d):"<<SE_position(next_idx,d)<<"\n";
                if(cost_exponent > position_gaussian_exponent_cuttoff) { //Too far away to be connected
                    feasible=false;
                    break;
                }
                C+= cost_exponent + log(dist_var);
//                 std::cout<<"cost: "<<C<<"\n";
            }
            if(!feasible) continue; //gaussian sigma constraint violated: move to next pair.
            if(maxSpeed>0 && sqrt(total_dist_sq)/deltaT > maxSpeed) continue; //maxSpeed constraint violated
            for(int f=0; f<nFeatures; f++){
                FloatT feat_var = featureVar(f) + SE_feature(cur_idx,f)+ SE_feature(next_idx,f);
                FloatT feat_dist = feature(cur_idx,f) - feature(next_idx,f);
                FloatT cost_exponent = feat_dist*feat_dist/feat_var;
                if(cost_exponent > feature_gaussian_exponent_cuttoff(f)) { //Too far away to be connected
                    feasible=false;
                    break;
                }
                C+= cost_exponent + log(feat_var);
            }
            if(!feasible) continue; //move to next pair.
            //Otherwise we have a valid cost so normalize and record it
            C+= norm_const;
            C*= 0.5;
            C-= log1mkoff;
//             std::cout<<"C:"<<C<<" log1mkoff:"<<log1mkoff<<"\n";
            //Record cost
            row_index.push_back(i);
            col_index.push_back(j);
            values.push_back(C);
            //Record lower right block dummy cost
            row_index.push_back(nCur+j);
            col_index.push_back(nNext+i);
            values.push_back(cost_epsilon);
        }
    }
    //Fill in death costs
    FloatT deathC= -logkoff;
    for(int i=0; i<nCur; i++){
        row_index.push_back(i);
        col_index.push_back(nNext+i);
        values.push_back(deathC);
    }
    //Fill in birth costs
    FloatT birthC = -logrho-logkon;
    for(int j=0; j<nNext; j++){
        row_index.push_back(nCur+j);
        col_index.push_back(j);
        values.push_back(birthC);
    }
    //Assemble sparse matrix
    int nnz = values.size();
    UMatT locations(2,nnz);
    for(int n=0; n<nnz;n++){
        locations(0,n) = row_index[n];
        locations(1,n) = col_index[n];
    }
    VecT values_vec(values.data(), nnz, false);//Re-use the vector's memory directly.
    bool sort_them=true; //Make sure armadillo sorts the locations
    bool check_for_zeros=false; //Don't bother checking for zeros
    return SpMatT(locations, values_vec, nTot, nTot, sort_them, check_for_zeros);
}

void LAPTrack::checkFrameIdxs()
{
    if(state!=F2F_LINKED) throw std::runtime_error("state != F2F_LINKED");
    int trackIdx=0;
    for(int n=firstFrame; n<=lastFrame; n++){
        int stidx = frameBirthStartIdx(n-firstFrame);
//         std::cout<<"Frame:"<<n<<" TrackIdx:"<<trackIdx<<" StartIdx:"<<stidx<<"\n";
        if(stidx!=trackIdx) throw std::runtime_error("startidx != trackidx");
        if(!(trackIdx == static_cast<int>(tracks.size()) ||  frameIdx(tracks[trackIdx].front()) >=n)) throw std::runtime_error("Track indexing error");
        while(trackIdx < static_cast<int>(tracks.size()) && frameIdx(tracks[trackIdx].front())==n) {
//             std::cout<<"TrackIdx:"<<trackIdx<<" Correctly begins at frame:"<<n<<"\n";
            trackIdx++;
        }
    }
}

void LAPTrack::debugCloseGaps(SpMatT &cost, IMatT &connections, VecT &conn_costs) const
{
    assert(state==F2F_LINKED);
    cost = computeGapCloseMatrix();
    IVecT track_assignment = LAP_JVSparse<FloatT>::solve(cost);
    int nTracks = static_cast<int>(tracks.size());
    int nCons = nTracks;
    for(int n=nTracks; n<2*nTracks; n++) 
        if(track_assignment(n) < nTracks) nCons++;
    connections.set_size(nCons,2);
    int conIdx=0;
    for(int n=0;n<2*nTracks;n++){
        if (n>=nTracks) {
            if (track_assignment(n)>=nTracks) continue; //phantom
            connections(conIdx,0) = -1; //birth
        } else {
            connections(conIdx,0) = n; //connection
        }
        connections(conIdx,1) = (n>=nTracks) ? -1 : track_assignment(n); //deaths
        conIdx++;
    }
    conn_costs = LAP_JVSparse<FloatT>::computeCost(cost, track_assignment);
    conn_costs = conn_costs.elem(arma::find(conn_costs>cost_epsilon));
}

void LAPTrack::closeGaps()
{
//     assert(checkFrameIdxs());
    //Invariant: tracks are in birth order.  So when connecting trackM->trackN we have M<N;
    assert(state==F2F_LINKED);
    auto cost = computeGapCloseMatrix();
//     arma::mat dC(cost);
//     std::cout<<"Cost: ("<<cost.n_rows<<","<<cost.n_cols<<"):\n"<<dC<<"\n";
    IVecT track_assignment = LAP_JVSparse<FloatT>::solve(cost);
    int nTracks = tracks.size();
    int nNewTracks = nTracks;
    for(int m=nTracks-1; m>=0; m--){ //start at the end.  Last track cannot connect so skip it.
        //we are considerting the connection for trackM -> trackN
        int n = track_assignment(m);
        assert(m<n); //either we don't connect or we connect to an N born after M. so M<N.
        if(n < nTracks) {
            //Valid track join operation.
            assert(~tracks[m].empty());
//             std::cout<<"Tracks(m="<<m<<"):"<<tracks[m].size()<<" Tracks(n="<<n<<"):"<<tracks[n].size()<<std::endl;
            tracks[m].splice(tracks[m].end(),tracks[n]);
            nNewTracks--;
//             std::cout<<"Joined "<<m<<"->"<<n<<" New num tracks:"<<nNewTracks<<"\n";
        }
    }
    //Compress down tracks.
    TrackVecT new_tracks(nNewTracks);
    //Remove empty tracks and any not meeting minFinalTrackLength
    std::copy_if(tracks.cbegin(), tracks.cend(), new_tracks.begin(), 
                    [this](const TrackT &t) {return t.size()>0 && (minFinalTrackLength<=1 || static_cast<int>(t.size())>minFinalTrackLength);});
    tracks = new_tracks;
    //These are all now invalid since tracks variable is changed
    trackAssignment.clear();
    birthFrameIdx.clear();
    frameBirthStartIdx.clear();
    state = GAPS_CLOSED;
}

LAPTrack::SpMatT 
LAPTrack::computeGapCloseMatrix() const
{
    int nTracks = static_cast<int>(tracks.size());
    
    //These will be our sparse matrix format vectors
    std::vector<arma::uword> row_index;
    std::vector<arma::uword> col_index;
    std::vector<FloatT> values;
    int reserve_size = nTracks*10; //Guesstimate amount of entries used
    row_index.reserve(reserve_size);
    col_index.reserve(reserve_size);
    values.reserve(reserve_size);

    FloatT position_gaussian_exponent_cuttoff = (maxPositionDisplacementSigma*maxPositionDisplacementSigma)/2.; //Only allow connections within 5 sigma
    VecT feature_gaussian_exponent_cuttoff = (maxFeatureDisplacementSigma%maxFeatureDisplacementSigma)/2.; //Only allow connections within 5 sigma
    
    FloatT norm_const = (nDims+nFeatures)*log2pi; //Pre-compute this
    FloatT birthC = -logrho-logkon;
    FloatT deathC= -logkoff;
//     std::cout<<"frameBirthStartIdx: "<<frameBirthStartIdx.t()<<"\n";
    
    //connect trackI to trackJ so trackJ must start after trackI ends.
    for(int i=0; i<nTracks; i++){
        if(static_cast<int>(tracks[i].size()) < minGapCloseTrackLength) continue; //Don't connect tracks shorter than minGapCloseTrackLength
        int locI = tracks[i].back(); //last localization for track I.
        int trackIend = frameIdx(locI); //frame death
        if (trackIend >= lastFrame-1) continue; //Tracks ending on last 2 frames can be a "start" point since this would be connected by F2F
        for(int j=frameBirthStartIdx(trackIend+2-firstFrame); j<nTracks; j++){
            if(static_cast<int>(tracks[j].size()) < minGapCloseTrackLength) continue; //Don't connect tracks shorter than minGapCloseTrackLength
            int trackJstart = birthFrameIdx[j];
//             std::cout<<"Track:"<<j<<" Birth Loc:"<<tracks[j].front()<<std::endl;
//             std::cout<<" Birth Frame Idx:"<<frameIdx(tracks[j].front())<<std::endl;
//             std::cout<<" Recorded: "<<birthFrameIdx[j]<<std::endl;
            int deltaT = trackJstart - trackIend;
//             std::cout<<"i("<<i<<") -> j("<<j<<"): endI:"<<trackIend<<" startJ:"<<trackJstart<<" deltaT:"<<deltaT<<"\n";
            assert(deltaT>=1);
            if(deltaT>=maxGapCloseFrames) continue; //Gap must be at most maxGapCloseFrames
            int locJ = tracks[j].front();
            FloatT DdT = 2*D*deltaT;
            FloatT total_dist_sq=0;
            FloatT C=0;
            bool feasible = true;
            for(int d=0; d<nDims; d++){
                FloatT dist_var = DdT + SE_position(locI,d) + SE_position(locJ,d);
                FloatT dist = position(locI,d) - position(locJ,d);
                FloatT dist_sq = dist*dist;
                total_dist_sq += dist_sq;
                FloatT cost_exponent = dist_sq/dist_var;
//                 std::cout<<"Dim:"<<d<<" dist:"<<dist<<" dist_var:"<<dist_var<<" costExp:"<<cost_exponent<<" ExpCuttoff:"<<gaussian_exponent_cuttoff<<"\n";
//                 std::cout<<"SE_position(cur_idx,d):"<<SE_position(cur_idx,d)<<"SE_position(next_idx,d):"<<SE_position(next_idx,d)<<"\n";
                if(cost_exponent > position_gaussian_exponent_cuttoff) { //Too far away to be connected
                    feasible=false;
                    break;
                }
                C+= cost_exponent + log(dist_var);
            }
            if(!feasible) continue; //gaussian sigma constraint violated: move to next pair.
            if(maxSpeed>0 && sqrt(total_dist_sq)/deltaT > maxSpeed) continue; //maxSpeed constraint violated
            for(int f=0; f<nFeatures; f++){
                FloatT feat_var = featureVar(f) + SE_feature(locI,f)+ SE_feature(locJ,f);
                FloatT feat_dist = feature(locI,f) - feature(locJ,f);
                FloatT cost_exponent = feat_dist*feat_dist/feat_var;
                if(cost_exponent > feature_gaussian_exponent_cuttoff(f)) { //Too far away to be connected
                    feasible=false;
                    break;
                }
                C+= cost_exponent + log(feat_var);
            }
            if(!feasible) continue; //move to next pair.
            //Otherwise we have a valid cost so normalize and record it
            C+= norm_const;
            C*= 0.5;
            C-= logkon +logkoff*deltaT;
            //Record cost
            row_index.push_back(i);
            col_index.push_back(j);
            values.push_back(C);
            //Record lower right block dummy cost
            row_index.push_back(nTracks+j);
            col_index.push_back(nTracks+i);
            values.push_back(cost_epsilon);
        }
    }
    
    for(int i=0; i<nTracks; i++){
        //Fill in death costs
        row_index.push_back(i);
        col_index.push_back(nTracks+i);
        values.push_back(deathC);
        //Fill in birth costs
        row_index.push_back(nTracks+i);
        col_index.push_back(i);
        values.push_back(birthC);
    }
    int nnz = values.size();
    UMatT locations(2,nnz);
    for(int n=0; n<nnz;n++){
        locations(0,n) = row_index[n];
        locations(1,n) = col_index[n];
    }
//     VecT values_vec(values.data(), nnz, false);//Re-use the vector's memory directly.
    bool sort_them=true; //Make sure armadillo sorts the locations
    bool check_for_zeros=false; //Don't bother checking for zeros
    return SpMatT(locations, VecT(values), 2*nTracks, 2*nTracks, sort_them, check_for_zeros);
}
