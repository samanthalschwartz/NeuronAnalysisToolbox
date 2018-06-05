/** @file SRRender.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-23-2014
 * @brief The class definition and template Specializations for SRRender.
 *
 * Rendering of SR emitter localizations
 */
#include "srrender.h"
#include <cassert>
#include "omp.h"

template<class FloatT>
const FloatT SRRender2D<FloatT>::DefaultSigmaAccuracy = 5.;

template<class FloatT>
const FloatT SRRender2D<FloatT>::normexp = 1/sqrt(2);


template<class FloatT>
void SRRender2D<FloatT>::renderHist(const EmitterVecT &points, const VecT &roi, ImageT &im)
{
//     #ifdef DEBUG
//     checkPoints(points);
//     #endif
    if (points.n_rows >= sqrt(im.n_elem)) {
        renderHistParallel(points,roi,im);
    } else {
        renderHistSingle(points,roi,im);
    }
}


template<class FloatT>
void SRRender2D<FloatT>::renderHistSingle(const EmitterVecT &points, const VecT &roi, ImageT &im)
{
    int pixelsX =  static_cast<int>(im.n_cols); //number of output pixels in the X direction (across rows)
    int pixelsY =  static_cast<int>(im.n_rows); //number of output pixels in the Y direction (down columns)
    FloatT xmin = roi(0);
    FloatT ymin = roi(2);
    FloatT sizeRatioX = static_cast<FloatT>(pixelsX) / (roi(1)-roi(0));
    FloatT sizeRatioY = static_cast<FloatT>(pixelsY) / (roi(3)-roi(2));
    for(unsigned n=0; n<points.n_rows; n++){
        int ix = static_cast<int>((points(n,1)-xmin)*sizeRatioX);
        int iy = static_cast<int>((points(n,2)-ymin)*sizeRatioY);
        if(0<=ix && ix<pixelsX && 0<=iy && iy<pixelsY) im(iy,ix) += points(n,0);
    }
}

template<class FloatT>
void SRRender2D<FloatT>::renderHistParallel(const EmitterVecT &points, const VecT &roi, ImageT &im)
{
    int pixelsX =  static_cast<int>(im.n_cols); //number of output pixels in the X direction (across rows)
    int pixelsY =  static_cast<int>(im.n_rows); //number of output pixels in the Y direction (down columns)
    FloatT xmin = roi(0);
    FloatT ymin = roi(2);
    FloatT sizeRatioX = static_cast<FloatT>(pixelsX) / (roi(1)-roi(0));
    FloatT sizeRatioY = static_cast<FloatT>(pixelsY) / (roi(3)-roi(2));

    int max_threads = omp_get_max_threads();
    arma::field<ImageT> histF(max_threads);
    int num_threads;
    #pragma omp parallel
    {
        ImageT hist(im.n_rows,im.n_cols);
        hist.zeros();
        num_threads=omp_get_num_threads(); //Save number of threads actually run
        #pragma omp for
        for(unsigned n=0; n<points.n_rows; n++){
            int ix = static_cast<int>((points(n,1)-xmin)*sizeRatioX);
            int iy = static_cast<int>((points(n,2)-ymin)*sizeRatioY);
            if(0<=ix && ix<pixelsX && 0<=iy && iy<pixelsY) hist(iy,ix) += points(n,0); //intensity
        }
        histF(omp_get_thread_num()) = hist;
    }
    //Parellelize sum of individual historgrams over columns
    #pragma omp parallel for
    for(int x=0; x<pixelsX; x++) for(int y=0; y<pixelsY; y++) {
        FloatT sum = 0.0;
        for(int n=0;n<num_threads; n++) sum += histF(n)(y,x);
        im(y,x) = sum;
    }
}


template<class FloatT>
void SRRender2D<FloatT>::renderHistMovie(const EmitterVecT &points, const VecT &roi, MovieT &im)
{
    int pixelsX =  static_cast<int>(im.n_cols); //number of output pixels in the X direction (across rows)
    int pixelsY =  static_cast<int>(im.n_rows); //number of output pixels in the Y direction (down columns)
    FloatT xmin = roi(0);
    FloatT ymin = roi(2);
    FloatT sizeRatioX = static_cast<FloatT>(pixelsX) / (roi(1)-roi(0));
    FloatT sizeRatioY = static_cast<FloatT>(pixelsY) / (roi(3)-roi(2));

    #pragma omp parallel
    {
        int num_threads = omp_get_num_threads(); //Number we actually created my be less than max
        int tid = omp_get_thread_num();
        for(unsigned n=0; n<points.n_rows; n++){
            int frame = points(n,5);
            if (frame%num_threads != tid) continue;
            int ix = static_cast<int>((points(n,1)-xmin)*sizeRatioX);
            int iy = static_cast<int>((points(n,2)-ymin)*sizeRatioY);
            if(0<=ix && ix<pixelsX && 0<=iy && iy<pixelsY) im(iy,ix,frame) += points(n,0);
        }
    }
}


template<class FloatT>
void SRRender2D<FloatT>::renderGauss(const EmitterVecT &points, const VecT &roi, ImageT &im, FloatT sigmaAccuracy)
{
//     #ifdef DEBUG
//     checkPoints(points);
//     #endif
    if (points.n_rows >= sqrt(im.n_elem)) {
        renderGaussParallel(points,roi,im,sigmaAccuracy);
    } else {
        renderGaussSingle(points,roi,im,sigmaAccuracy);
    }
}

template<class FloatT>
void SRRender2D<FloatT>::renderGaussSingle(const EmitterVecT &points, const VecT &roi, ImageT &im, FloatT sigmaAccuracy)
{
    int pixelsX =  static_cast<int>(im.n_cols); //number of output pixels in the X direction (across rows)
    int pixelsY =  static_cast<int>(im.n_rows); //number of output pixels in the Y direction (down columns)
    FloatT imageXmin = roi(0);
    FloatT imageYmin = roi(2);
    FloatT sizeRatioX = static_cast<FloatT>(pixelsX) / (roi(1)-roi(0));
    FloatT sizeRatioY = static_cast<FloatT>(pixelsY) / (roi(3)-roi(2));
    int N = static_cast<int>(points.n_rows);
    VecT xStencil(pixelsX), yStencil(pixelsY);
    im.zeros();
    for(int n=0; n<N; n++) {
        FloatT X = (points(n,1)-imageXmin)*sizeRatioX;
        FloatT sigmaX = points(n,3)*sizeRatioX;
        int xp = static_cast<int>(X); // 0 <= xp < pixelsX
        int xhw = static_cast<int>(0.5+sigmaAccuracy*sigmaX); //halfwidth for gaussian X
        int xmin = std::max(0,xp-xhw);
        int xmax = std::min(pixelsX-1,xp+xhw);
        int xspn = xmax-xmin+1;

        FloatT Y = (points(n,2)-imageYmin)*sizeRatioY;
        FloatT sigmaY = points(n,4)*sizeRatioY;
        int yp = static_cast<int>(Y); // 0 <= yp < pixelsY
        int yhw = static_cast<int>(0.5+sigmaAccuracy*sigmaY); //halfwidth for gauissian Y
        int ymin = std::max(0,yp-yhw);
        int ymax = std::min(pixelsY-1,yp+yhw);
        int yspn = ymax-ymin+1;

        if(xspn<=0 || yspn<=0) continue;
        fill_stencil(xspn, X-xmin, sigmaX, xStencil);
        fill_stencil(yspn, Y-ymin, sigmaY, yStencil);
        FloatT I = points(n,0);
        for(int x=0; x<xspn; x++) xStencil(x)*=I; //Pre-multiply by I;
        //Copy in image of new gaussian
        for(int x=xmin; x<=xmax; x++) for(int y=ymin; y<=ymax; y++)
            im(y,x) += xStencil(x-xmin) * yStencil(y-ymin);
    }
}


template<class FloatT>
void SRRender2D<FloatT>::renderGaussParallel(const EmitterVecT &points, const VecT &roi, ImageT &final_image, FloatT sigmaAccuracy)
{
    int pixelsX =  static_cast<int>(final_image.n_cols); //number of output pixels in the X direction (across rows)
    int pixelsY =  static_cast<int>(final_image.n_rows); //number of output pixels in the Y direction (down columns)
    FloatT imageXmin = roi(0);
    FloatT imageYmin = roi(2);
    FloatT sizeRatioX = static_cast<FloatT>(pixelsX) / (roi(1)-roi(0));
    FloatT sizeRatioY = static_cast<FloatT>(pixelsY) / (roi(3)-roi(2));
    int N = static_cast<int>(points.n_rows);
    int max_threads = omp_get_max_threads();
    int num_threads;
    arma::field<ImageT> imStack(max_threads);
    #pragma omp parallel
    {
        ImageT im(pixelsY,pixelsX,arma::fill::zeros);
        num_threads = omp_get_num_threads(); //Number we actually created my be less than max
        VecT xStencil(pixelsX), yStencil(pixelsY);
        num_threads = omp_get_num_threads(); //Number we actually created my be less than max
        #pragma omp for
        for(int n=0; n<N; n++) {
            FloatT X = (points(n,1)-imageXmin)*sizeRatioX;
            FloatT sigmaX = points(n,3)*sizeRatioX;
            int xp = static_cast<int>(X); // 0 <= xp < pixelsX
            int xhw = static_cast<int>(0.5+sigmaAccuracy*sigmaX); //halfwidth for gaussian X
            int xmin = std::max(0,xp-xhw);
            int xmax = std::min(pixelsX-1,xp+xhw);
            int xspn = xmax-xmin+1;

            FloatT Y = (points(n,2)-imageYmin)*sizeRatioY;
            FloatT sigmaY = points(n,4)*sizeRatioY;
            int yp = static_cast<int>(Y); // 0 <= yp < pixelsY
            int yhw = static_cast<int>(0.5+sigmaAccuracy*sigmaY); //halfwidth for gauissian Y
            int ymin = std::max(0,yp-yhw);
            int ymax = std::min(pixelsY-1,yp+yhw);
            int yspn = ymax-ymin+1;

            if(xspn<=0 || yspn<=0) continue;
            fill_stencil(xspn, X-xmin, sigmaX, xStencil);
            fill_stencil(yspn, Y-ymin, sigmaY, yStencil);
            FloatT I = points(n,0);
            for(int x=0; x<xspn; x++) xStencil(x)*=I; //Pre-multiply by I;
            //Copy in image of new gaussian
            for(int x=xmin; x<=xmax; x++) for(int y=ymin; y<=ymax; y++)
                im(y,x) += xStencil(x-xmin) * yStencil(y-ymin);
        }
        imStack(omp_get_thread_num()) = im;
    }
    //Parellelize sum of individual historgrams over columns
    #pragma omp parallel for
    for(int x=0; x<pixelsX; x++) for(int y=0; y<pixelsY; y++) {
        FloatT sum = 0.0;
        for(int n=0;n<num_threads; n++) sum += imStack(n)(y,x);
        final_image(y,x) = sum;
    }
}

template<class FloatT>
void SRRender2D<FloatT>::renderGaussMovie(const EmitterVecT &points, const VecT &roi, MovieT &im, FloatT sigmaAccuracy)
{
    int pixelsX =  static_cast<int>(im.n_cols); //number of output pixels in the X direction (across rows)
    int pixelsY =  static_cast<int>(im.n_rows); //number of output pixels in the Y direction (down columns)
    FloatT imageXmin = roi(0);
    FloatT imageYmin = roi(2);
    FloatT sizeRatioX = static_cast<FloatT>(pixelsX) / (roi(1)-roi(0));
    FloatT sizeRatioY = static_cast<FloatT>(pixelsY) / (roi(3)-roi(2));
    int N = static_cast<int>(points.n_rows);
    #pragma omp parallel
    {
        VecT xStencil(pixelsX), yStencil(pixelsY);
        int num_threads = omp_get_num_threads(); //Number we actually created my be less than max
        int tid = omp_get_thread_num();
        for(int n=0; n<N; n++) {
            int frame = points(n,5);
            if (frame%num_threads != tid) continue;
            FloatT X = (points(n,1)-imageXmin)*sizeRatioX;
            FloatT sigmaX = points(n,3)*sizeRatioX;
            int xp = static_cast<int>(X); // 0 <= xp < pixelsX
            int xhw = static_cast<int>(0.5+sigmaAccuracy*sigmaX); //halfwidth for gaussian X
            int xmin = std::max(0,xp-xhw);
            int xmax = std::min(pixelsX-1,xp+xhw);
            int xspn = xmax-xmin+1;

            FloatT Y = (points(n,2)-imageYmin)*sizeRatioY;
            FloatT sigmaY = points(n,4)*sizeRatioY;
            int yp = static_cast<int>(Y); // 0 <= yp < pixelsY
            int yhw = static_cast<int>(0.5+sigmaAccuracy*sigmaY); //halfwidth for gauissian Y
            int ymin = std::max(0,yp-yhw);
            int ymax = std::min(pixelsY-1,yp+yhw);
            int yspn = ymax-ymin+1;

            if(xspn<=0 || yspn<=0) continue;
            fill_stencil(xspn, X-xmin, sigmaX, xStencil);
            fill_stencil(yspn, Y-ymin, sigmaY, yStencil);
            FloatT I = points(n,0);
            for(int x=0; x<xspn; x++) xStencil(x)*=I; //Pre-multiply by I;
            //Copy in image of new gaussian
            for(int x=xmin; x<=xmax; x++) for(int y=ymin; y<=ymax; y++)
                im(y,x,frame) += xStencil(x-xmin) * yStencil(y-ymin);
        }
    }
}


template<class FloatT>
void SRRender2D<FloatT>::fill_stencil(int size, FloatT x, FloatT sigma, VecT& stencil)
{
    FloatT norm = normexp/sigma;
    FloatT derf = erf(-norm*x);
    for(int i=0;i<size;i++) {
        FloatT last_derf = derf;
        derf = erf(norm*((i+1)-x));
        stencil(i) = 0.5*(derf-last_derf);
    }
}
/*
template<class FloatT>
void SRRender2D<FloatT>::checkPoints(const EmitterVecT &points) const
{
    int nPoints = points.n_rows;
    for(int n=0; n<nPoints; n++){
        assert(points(n,0)>0);
        assert(points(n,1)>=0 && points(n,1)<size(0));
        assert(points(n,2)>=0 && points(n,2)<size(1));
        assert(points(n,3)>0);
        assert(points(n,4)>0);
    }
}*/

/* Explicit Template Instantiation */
template class SRRender2D<float>;
template class SRRender2D<double>;


/*
template<class FloatT>
SRRenderHS<FloatT>::SRRenderHS(const VecT &_size) : size(_size) {}



template<class FloatT>
void SRRenderHS<FloatT>::renderHist(const EmitterVecT &points, ImageT &im)
{
    for(unsigned n=0; n<points.n_cols; n++) im(points(n,2), points(n,3), points(n,4)) += points(n,1);
}*/
