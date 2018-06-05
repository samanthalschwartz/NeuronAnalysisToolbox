/**
 * @file Boxxer3D.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-24-2014
 * @brief The class method definitions for Boxxer3D.
 */

#include <cassert>
#include "GaussFilter.h"
#include "Maxima.h"
#include "Boxxer3D.h"
#include "omp.h"

/* Static member variables */
template<class FloatT>
const int Boxxer3D<FloatT>::dim = 3;
template<class FloatT>
const FloatT Boxxer3D<FloatT>::DefaultSigmaRatio = 1.1;

template<class FloatT>
Boxxer3D<FloatT>::Boxxer3D(const IVecT &imsize, const MatT &_sigma)
    : nScales(_sigma.n_cols),imsize(imsize), sigma(_sigma), sigma_ratio(DefaultSigmaRatio)
{
    assert(nScales>=1);
    assert(static_cast<int>(imsize.n_elem)==dim);
    assert(static_cast<int>(sigma.n_rows)==dim);
}

template<class FloatT>
void Boxxer3D<FloatT>::setDoGSigmaRatio(FloatT _sigma_ratio)
{
    assert(_sigma_ratio>1);
    sigma_ratio=_sigma_ratio;
}

template<class FloatT>
void Boxxer3D<FloatT>::filterScaledLoG(const ImageT &im, ScaledImageT &fim)
{
    #pragma omp parallel for
    for(int s=0; s<nScales; s++) {
        LoGFilter3D<FloatT> scale_filter(imsize,sigma.col(s));
        scale_filter.filter(im,fim.slice(s));
    }
}

template<class FloatT>
void Boxxer3D<FloatT>::filterScaledDoG(const ImageT &im, ScaledImageT &fim)
{
    #pragma omp parallel for
    for(int s=0; s<nScales; s++) {
        DoGFilter3D<FloatT> scale_filter(imsize,sigma.col(s),sigma_ratio);
        scale_filter.filter(im,fim.slice(s));
    }
}

/**
 * 
 * Get the maxima over all sacales and all frames.  Scale and maxfind on each frame individually to
 * cut down on memory size (otherwise it would be easier to decouple the filtering and maxfinding.
 */

template<class FloatT>
int Boxxer3D<FloatT>::scaleSpaceLoGMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, 
                                          int neighborhood_size, int scale_neighborhood_size)
{
    int nT=static_cast<int>(im.n_slices);
    arma::field<IMatT> frame_maxima(nT); //These will come back 3xN
    arma::field<VecT> frame_max_vals(nT);
    #pragma omp parallel
    {
        auto sim = make_scaled_image();
        std::vector<LoGFilter3D<FloatT>> scale_filters;
        for(int s=0; s<nScales; s++) scale_filters.push_back(LoGFilter3D<FloatT>(imsize,sigma.col(s)));
        #pragma omp for
        for(int n=0; n<nT; n++) {
            for(int s=0; s<nScales; s++) scale_filters[s].filter(im.slice(n),sim.slice(s));
            scaleSpaceFrameMaxima(sim, frame_maxima(n), frame_max_vals(n), neighborhood_size, scale_neighborhood_size);
        }
    }
    return combine_maxima(frame_maxima, frame_max_vals, maxima, max_vals);
}

template<class FloatT>
int Boxxer3D<FloatT>::scaleSpaceDoGMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, 
                                          int neighborhood_size, int scale_neighborhood_size)
{
    int nT=static_cast<int>(im.n_slices);
    arma::field<IMatT> frame_maxima(nT); //These will come back 3xN
    arma::field<VecT> frame_max_vals(nT);
    #pragma omp parallel
    {
        auto sim = make_scaled_image();
        std::vector<DoGFilter3D<FloatT>> scale_filters;
        for(int s=0; s<nScales; s++) scale_filters.push_back(DoGFilter3D<FloatT>(imsize,sigma.col(s),sigma_ratio));
        #pragma omp for
        for(int n=0; n<nT; n++) {
            for(int s=0; s<nScales; s++) scale_filters[s].filter(im.slice(n),sim.slice(s));
            scaleSpaceFrameMaxima(sim, frame_maxima(n), frame_max_vals(n), neighborhood_size, scale_neighborhood_size);
        }
    }
    return combine_maxima(frame_maxima, frame_max_vals, maxima, max_vals);
}

/**
 * Get the scale maxima for a single frame
 */
template<class FloatT>
int Boxxer3D<FloatT>::scaleSpaceFrameMaxima(const ScaledImageT &sim, IMatT &maxima, VecT &max_vals, 
                                   int neighborhood_size, int scale_neighborhood_size) const
{
    arma::field<IMatT> scale_maxima(nScales);
    arma::field<VecT> scale_max_vals(nScales);
    Maxima3D<FloatT> maxima3D(imsize, neighborhood_size);
    for(int s=0; s<nScales; s++) 
        maxima3D.find_maxima(sim.slice(s), scale_maxima(s), scale_max_vals(s));
    combine_maxima(scale_maxima, scale_max_vals, maxima, max_vals);
    return scaleSpaceFrameMaximaRefine(sim, maxima, max_vals, scale_neighborhood_size);
}

/**
 * Given a scaled image and scale maxima, refine to remove overlapping scale maxima 
 */
template<class FloatT>
int
Boxxer3D<FloatT>::scaleSpaceFrameMaximaRefine(const ScaledImageT &im, IMatT &maxima, VecT &max_vals, 
                                              int scale_neighborhood_size) const
{
    using std::min;
    using std::max;
    IMatT new_maxima(maxima.n_rows, maxima.n_cols);
    VecT new_max_vals(max_vals.n_elem);
    int nMaxima = static_cast<int>(maxima.n_cols);
    int nNewMaxima=0;
    int delta = static_cast<int>((scale_neighborhood_size-1)/2);
    for(int n=0; n<nMaxima; n++) {
        IVecT mx = maxima.col(n);
        double mxv = max_vals(n);
        bool ok=true;
        if ( (mx(0)-delta < 0 || mx(0)+delta>=imsize(0)) || 
             (mx(1)-delta < 0 || mx(1)+delta>=imsize(1)) || 
             (mx(2)-delta < 0 || mx(2)+delta>=imsize(2))) {
            for(int s=0; s<nScales; s++)
                for(int k=max(0,mx(2)-delta); k<=min(imsize(2)-1,mx(2)+delta); k++) 
                    for(int j=max(0,mx(1)-delta); j<=min(imsize(1)-1,mx(1)+delta); j++) 
                        for(int i=max(0,mx(0)-delta); i<=min(imsize(0)-1,mx(0)+delta); i++) {
                            if( im(i,j,k,s) > mxv) {ok=false; goto done;}
            }
        } else {
            for(int s=0; s<nScales; s++)
                for(int k=mx(2)-delta; k<=mx(2)+delta; k++) 
                    for(int j=mx(1)-delta; j<=mx(1)+delta; j++) 
                        for(int i=mx(0)-delta; i<=mx(0)+delta; i++) {
                            if( im(i,j,k,s) > mxv) {ok=false; goto done;}
            }
        }
        done:
        if(ok){
            new_maxima.col(nNewMaxima)=mx;
            new_max_vals(nNewMaxima)=mxv;
            nNewMaxima++;
        }
    }
    maxima = new_maxima(arma::span::all, arma::span(0,nNewMaxima-1));
    max_vals = new_max_vals(arma::span(0,nNewMaxima-1));
    return nNewMaxima;
}

/* Static Methods */


template<class FloatT>
void Boxxer3D<FloatT>::filterLoG(const ImageStackT &im, ImageStackT &fim, const VecT &sigma)
{
    int nT=static_cast<int>(fim.n_slices);
    IVecT imsize={im.sX,im.sY,im.sZ};
    #pragma omp parallel
    {
        LoGFilter3D<FloatT> filter(imsize,sigma);
        #pragma omp for
        for(int n=0; n<nT; n++) 
            filter.filter(im.slice(n),fim.slice(n));
    }
}

template<class FloatT>
void Boxxer3D<FloatT>::filterDoG(const ImageStackT &im, ImageStackT &fim, const VecT &sigma, FloatT sigma_ratio)
{
    int nT=static_cast<int>(fim.n_slices);
    IVecT imsize={im.sX,im.sY,im.sZ};
    #pragma omp parallel
    {
        DoGFilter3D<FloatT> filter(imsize,sigma,sigma_ratio);
        #pragma omp for
        for(int n=0; n<nT; n++) 
            filter.filter(im.slice(n),fim.slice(n));
    }
}

template<class FloatT>
void Boxxer3D<FloatT>::filterGauss(const ImageStackT &im, ImageStackT &fim, const VecT &sigma)
{
    int nT=static_cast<int>(fim.n_slices);
    IVecT imsize={im.sX,im.sY,im.sZ};
    #pragma omp parallel
    {
        GaussFilter3D<FloatT> filter(imsize,sigma);
        #pragma omp for
        for(int n=0; n<nT; n++) 
            filter.filter(im.slice(n),fim.slice(n));
    }
}

/**
 * This finds local maxima over an image stack in parallel.
 */
template<class FloatT>
int Boxxer3D<FloatT>::enumerateImageMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, 
                                           int neighborhood_size)
{
    int nT=static_cast<int>(im.n_slices);
    arma::field<IMatT> frame_maxima(nT);
    arma::field<VecT> frame_max_vals(nT);
    IVecT imsize={im.sX,im.sY,im.sZ};
    #pragma omp parallel
    {
        Maxima3D<FloatT> maxima3D(imsize, neighborhood_size);
        #pragma omp for
        for(int n=0; n<nT; n++) {
            maxima3D.find_maxima(im.slice(n), frame_maxima(n), frame_max_vals(n));
        }
    }
    return combine_maxima(frame_maxima, frame_max_vals, maxima, max_vals);
}

template<class FloatT>
int Boxxer3D<FloatT>::combine_maxima(const arma::field<IMatT> &frame_maxima, 
                                     const arma::field<VecT> &frame_max_vals,
                                     IMatT &maxima, VecT &max_vals)
{
    unsigned Nmaxima=0;
    for(unsigned n=0; n<frame_max_vals.n_elem; n++) Nmaxima += frame_max_vals(n).n_elem;
    int nrows = frame_maxima(0).n_rows;
    assert(nrows>=3);
    maxima.resize(nrows+1,Nmaxima);
    max_vals.resize(Nmaxima);
    unsigned Nsaved=0;
    for(unsigned n=0; n<frame_max_vals.n_elem; n++) {
        unsigned NFrameMaxima = frame_max_vals(n).n_elem;
        if(NFrameMaxima>0){
            maxima(arma::span(0,nrows-1), arma::span(Nsaved,Nsaved+NFrameMaxima-1))=frame_maxima(n);
            maxima(arma::span(nrows,nrows), arma::span(Nsaved,Nsaved+NFrameMaxima-1)).fill(n);
            max_vals.rows(Nsaved,Nsaved+NFrameMaxima-1)=frame_max_vals(n);
            Nsaved+=NFrameMaxima;
        }
    }
    return Nmaxima;
}


template<class FloatT>
void Boxxer3D<FloatT>::checkMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals)
{
    int Nmaxima=static_cast<int>(maxima.n_cols);
    for(int n=0; n<Nmaxima; n++){
        FloatT val=im(maxima(0,n), maxima(1,n), maxima(2,n), maxima(3,n));
        if (val!=max_vals(n)) {
            printf(" (%i,%i,%i,%i):%.9f!= %.9f\n",maxima(0,n), maxima(1,n), maxima(2,n), maxima(3,n), val, max_vals(n));
        }
    }
}



/* Explicit Template Instantiation */
template class Boxxer3D<float>;
template class Boxxer3D<double>;

