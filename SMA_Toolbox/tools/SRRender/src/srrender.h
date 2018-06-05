/** @file SRRender.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-23-2014
 * @brief The class declaration and inline and templated functions for SRRender.
 *
 * Rendering of SR emitter localizations
 */

#ifndef _SRRENDER_H
#define _SRRENDER_H

#include <armadillo>

/**
 * 
 * Points format.  Row-oriented each row is a point, each column is a property
 * 2D renderHist Columns: [I X Y]
 * 2D renderGauss Columns:[I X Y sigmaX sigmaY]
 * 2D renderHistMovie Columns: [I X Y Frame]  - Frame is 0-indexed
 * 2D renderGaussMovie Columns:[I X Y sigmaX sigmaY Frame] - Frame is 0-indexed
 * 
 * The 'size' parameter gives the size of the entire field of view to be rendered in units 
 * correpsonding to the points format vectors. 
 * 
 */

template<class FloatT>
class SRRender2D{
public:
    typedef arma::Col<int32_t> IVecT;
    typedef arma::Col<FloatT> VecT;
    typedef arma::Mat<FloatT> ImageT;
    typedef arma::Cube<FloatT> MovieT;
    typedef arma::Mat<FloatT> EmitterVecT;
    static const FloatT DefaultSigmaAccuracy; //Default number of sigmas to render gaussian at

    static void renderHist(const EmitterVecT &points, const VecT &roi, ImageT &im);
    static void renderGauss(const EmitterVecT &points, const VecT &roi, ImageT &im, FloatT sigmaAccuracy=DefaultSigmaAccuracy);
    static void renderHistMovie(const EmitterVecT &points, const VecT &roi, MovieT &im);
    static void renderGaussMovie(const EmitterVecT &points, const VecT &roi, MovieT &im, FloatT sigmaAccuracy=DefaultSigmaAccuracy);
//     static void checkPoints(const EmitterVecT &points);
private:
    static const FloatT normexp; // 1/sqrt(2);

    static void fill_stencil(int size, FloatT x, FloatT sigma, VecT& stencil);
    static void renderHistSingle(const EmitterVecT &points, const VecT &roi, ImageT &im);
    static void renderHistParallel(const EmitterVecT &points, const VecT &roi, ImageT &im);
    static void renderGaussSingle(const EmitterVecT &points, const VecT &roi, ImageT &im, FloatT sigmaAccuracy);
    static void renderGaussParallel(const EmitterVecT &points, const VecT &roi, ImageT &im, FloatT sigmaAccuracy);
};


/* Inlined Methods */


#endif /* SRRENDER_H */
