
#include<iostream>
#include "LAPTrack.h"
#include "LAP_JVSparse.h"
// typedef double FloatT;
// typedef int32_t IndexT;
// typedef arma::SpMat<FloatT> SpMatT;
// typedef arma::Col<FloatT> VecT;
// typedef arma::Col<IndexT> IVecT;
// typedef arma::Mat<IndexT> IMatT;

using namespace arma;
using namespace std;

void testLAP()
{
    mat C;
    C<<11.1<<0<<5<<3<<9<<3<<endr
     <<5<<0<<0<<2<<1<<6<<endr
     <<0<<0<<1<<15<<10<<7<<endr
     <<7.1<<7.2<<7.3<<7.4<<7.5<<7.6<<endr
     <<3<<1<<1<<0<<0<<6<<endr
     <<0<<6<<3<<4<<0<<0<<endr;
    sp_mat Csp(C);    
    auto row_sol = LAP_JVSparse<double>::solve(Csp);
    std::cout<<"RowSol: "<<row_sol.t()<<"\n";
}

void testTracking()
{
    int Nframes = 20;
    auto frameIdx = randi<Tracker::IVecT>(Nframes,distr_param(1, 3*Nframes));
    frameIdx+=1000;
    frameIdx=unique(sort(frameIdx));
    int k = static_cast<int>(frameIdx.n_elem);
    if( k < Nframes) {
        frameIdx.resize(Nframes);
        for(int i=k; i<Nframes; i++) frameIdx(i)=frameIdx(i-1)+1;
    }
    std::cout<<"FrameIdx: "<<frameIdx.t()<<"\n";
    
    mat position = randn<mat>(Nframes,2)*3;
    mat SE_position = randn<mat>(Nframes,2)*0.3;
    std::cout<<"Position: "<<position<<"\n";
    std::cout<<"SEPosition: "<<SE_position<<"\n";
    Tracker::VecParamT params;
    params["D"]=0.3;
    params["kon"]=0.1;
    params["koff"]=0.1;
    params["rho"]=0.02;
    params["maxSpeed"] = -1;
    params["maxPositionDisplacementSigma"] = 5;
    params["maxFeatureDisplacementSigma"] = 5;
    params["maxGapCloseFrames"] = 5;
    params["minGapCloseTrackLength"] = 1;
    params["minFinalTrackLength"] = 1;
    params["featureVar"] = {2, 3};
    LAPTrack tracker(params);
    tracker.initializeTracks(frameIdx, position, SE_position);
    LAPTrack::SpMatT C;
    Tracker::IVecT cur_locs, next_locs;
    Tracker::VecT conn_costs;
    Tracker::IMatT connections;
    tracker.debugF2F(frameIdx(0),cur_locs, next_locs,C,connections, conn_costs);
    mat bigC(C);
    std::cout<<"CurLocs:\n"<<cur_locs.t()<<"\n";
    std::cout<<"NextLocs:\n"<<next_locs.t()<<"\n";
    std::cout<<"CostMat:\n"<<bigC<<"\n";
    std::cout<<"Connections: \n"<<connections<<"\n";
    std::cout<<"ConnectionCosts: \n"<<conn_costs<<"\n";
    tracker.printTracks();
    tracker.linkF2F();
    tracker.printTracks();
    tracker.closeGaps();
    tracker.printTracks();
// //     tracker.getTracks();
}

int main()
{
    testLAP();
    cout<<" =========== TRACKING ====================\n";
    testTracking();
    return 0;
}

/*

FloatT lap_arma(const SpMat &C, IVecT &x, IVecT &y, VecT &u VecT &v)
{

//     IndexT h,l,t,last,tel,td1=0,td2,i0,j0=0,j1=0,l0;
// 
//     IndexT *lab, *freeRow, *todo;
//     bool *ok;
//     FloatT min, v0, vj, dj, tmp;
//     FloatT *d;

    FloatT FLT_EPSILON = std::numeric_limits<FloatT>::epsilon();

    bool *ok = new bool[n + 1];
    IndexT *lab = new IndexT[n + 2];
    IndexT *freeRow = new IndexT[n + 2];
    IndexT *todo = new IndexT[n + 2];
    FloatT *d = new FloatT[n + 2];

    // Initialize
    int N = static_cast<IndexT>(C.n_rows);
    
    v.fill(INFINITY);
    x.zeros();
    u.zeros();
    const FloatT *cc=C.values;
    const IndexT *kk=C.row_indices;
    const IndexT *first=C.col_ptrs;
    
    for(IndexT i=0; i<n; i++) for(IndexT t=first[i]; t<first[i+1]; t++){
        IndexT j = kk[t]; //i=row, j=column, t=linear_index
        if (cc[t] < v[j]) {
            v[j] = cc[t];
            y[j] = i;
        }
    }
    
    // COLUMN REDUCTION 
    for(IndexT j=n-1; j>=0; j--) {
        IndexT i = y[j]; 
        if(x[i] == 0) {
            x[i] = j;
        } else {
            y[j] = 0;
            x[i] = -std::abs(x[i]);
        }
    }
    
    // REDUCTION TRANSFER
    IndexT freeRowCnt = 0;
    std::vector<IndexT> freeRows();
    for(IndexT i=0; i<n; i++) {
        if (x[i] < 0) x[i] = -x[i];
        else if (x[i] == 0) freeRows.push_back(i); 
        else {
            IndexT min_val = INFINITY;
            IndexT j1 = x[i];
            for (IndexT t=first[i]; t<first[i+1]; t++) {
                IndexT j = kk[t];
                if (j != j1) min_val = std::min(min_val, cc[t] - v[j])
            }
            u[i] = min_val;
            assert(C(i,j1));
            v[j1] = C(i,j1) - min_val;
        }
    }

    // AUGMENTING ROW REDUCTION 
    // Improve initial solution 
    for (int twice = 0; twice < 2; twice++) {
        IndexT ctr = 0;
        for(auto i: freeRows) {
            IndexT v0 = INFINITY;
            IndexT vj = INFINITY;
            for (IndexT t = first[i]; t < first[i+1]; t++) {
                IndexT j = kk[t];//i - row, j - col
                IndexT dj = cc[t] - v[j];
                if (dj < vj) {
                    if (dj >= v0) {
                        vj = dj;
                        j1 = j;
                    } else {
                        vj = v0;
                        v0 = dj;
                        j1 = j0;
                        j0 = j;
                    }
                }
            }

            IndexT i0 = y[j0];
            u[i] = vj;
            if (vj - v0 > FLT_EPSILON) {
                v[j0] = v[j0] - vj + v0;
            } else if (i0 > 0) {
                j0 = j1;
                i0 = y[j0];
            } 
            x[i] = j0;
            y[j0] = i;

            if (i0 > 0) {
                if (vj - v0 > FLT_EPSILON) {
                    freeRow[--h] = i0;
                } else {
                    freeRow[++l] = i0;
                }
            }
        } 
    }

    tmp = 0;
    for (i = 1; i <= n; i++) {
        tmp += u[i] + v[i];
    }

    // Augmentation part
    l0 = l;
    for (l = 1; l <= l0; l++) {

        for (j = 1; j <= n; j++) {
            d[j] = INFINITY;
            ok[j] = false;
        }

        min = INFINITY; i0 = freeRow[l];

        for (t = first[i0]; t < first[i0+1]; t++) {
            j = kk[t];
            dj = cc[t] - v[j];
            d[j] = dj;
            lab[j] = i0;

            if (dj <= min) {
                if (dj < min) {
                    td1 = 0;
                    min = dj;
                }
                todo[++td1] = j;
            }
        }

        for (h = 1; h <= td1; h++) {
            j = todo[h];
            if (y[j] == 0) {
                goto label2;
            } 
            ok[j] = true;
        }

        td2 = n;
        last = n + 1;

        // Repeat until a freeRow row found 
        while (true) {
            j0 = todo[td1--];
            i = y[j0];
            todo[td2--] = j0;
            t = first[i];

            for (t = first[i]; kk[t] != j0; t++) {
                // nothing 
            }

            tmp = cc[t] - v[j0] - min;

            for (t = first[i]; t < first[i+1]; t++) {
                j = kk[t];
                if (!ok[j]) {
                    vj = cc[t] - v[j] - tmp;
                    if (vj < d[j]) {
                        d[j] = vj;
                        lab[j] = i;
                        if (vj == min) {
                            if (y[j] == 0) {
                                goto label1;
                            }
                            td1++;
                            todo[td1] = j;
                            ok[j] = true;
                        }
                    }
                }
            }
            if (td1 == 0) {
                min = INFINITY - 1;
                last = td2 + 1;
                for (j = 1; j <= n; j++) {
                    if (d[j] <= min) {
                        if (!ok[j]) {
                            if (d[j] < min) {
                                td1 = 0;
                                min = d[j];
                            }
                            todo[++td1] = j;
                        }
                    }
                }
                for (h = 1; h <= td1; h++) {
                    j = todo[h];
                    if (y[j] == 0) {
                        goto label1;
                    }
                    ok[j] = true;
                }
            }
        } 
        label1:
        for (k = last; k <= n; k++) {
            j0 = todo[k];
            v[j0] += d[j0] - min;
        }

        label2:
        do {
            i = lab[j];
            y[j] = i;
            k = j;
            j = x[i];
            x[i] = k;
        } while (i != i0);
    }

    tmp = 0;
    for (i = 1; i <= n; i++) {
        j  = x[i];
        t = first[i];
        while (kk[t] != j) {
            t++;
        } 

        u[i] = cc[t] - v[j];
        tmp += cc[t];
    }


    delete [] ok;
    delete [] lab;
    delete [] freeRow;
    delete [] todo;
    delete [] d;

    return(tmp);
}*/
