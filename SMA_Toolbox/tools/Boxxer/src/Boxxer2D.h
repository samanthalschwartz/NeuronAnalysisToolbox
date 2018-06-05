/**
* @file Boxxer2D.h
* @author Mark J. Olah (mjo\@cs.unm.edu)
* @date 07-28-2014
* @brief The class declaration for Boxxer2D.
*/
#ifndef _BOXXER2D_H
#define _BOXXER2D_H

#include <armadillo>
#include "hypercube.h"
/**
 * @class Boxxer2D
 * 
 * In this class we make the assumtion that images are stored in column-major format
 * and that x=rows, y=cols, t=slices.  This relationship is important in the choice of imsize and sigma
 * parameters.  
 * 
 * imsize = [nrows, ncols, nframes];
 * sigma = [ sigma_rows (X scale=1), sigma_rows (X scale=2); 
 *           sigma_cols(Y scale=1), sigma_cols (Y scale=2)]
 *
 * This is contrary to normal image coordinates in matlab, but for this low level it is easier to think
 * about x as the first index into an image and understand that the meaning for "X" and "Y" will be reversed
 * from the matlab interpretation, but only internally within the Boxxer iface.
 * 
 */
template<class FloatT>
class Boxxer2D
{
public:
    typedef arma::Col<int> IVecT;
    typedef arma::Mat<int> IMatT;
    typedef arma::Col<FloatT> VecT;
    typedef arma::Mat<FloatT> MatT;
    typedef arma::Mat<FloatT> ImageT;
    typedef arma::Cube<FloatT> ImageStackT;
    typedef arma::Cube<FloatT> ScaledImageT;
    typedef Hypercube<FloatT>  ScaledImageStackT;
 
    static const FloatT DefaultSigmaRatio;
    static const int dim;
    
    int nScales;
    IVecT imsize; // [nrows x ncols] size of an individual frame
    MatT sigma; // size: [2 x nScales] row1=sigmaX (rows), row2=sigmaY (cols)
    FloatT sigma_ratio;
    Boxxer2D(const IVecT &imsize, const MatT &sigma);

    void setDoGSigmaRatio(FloatT sigma_ratio);

    void filterScaledLoG(const ImageStackT &im, ScaledImageStackT &fim) const;
    void filterScaledDoG(const ImageStackT &im, ScaledImageStackT &fim) const;
    int scaleSpaceLoGMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int neighborhood_size, int scale_neighborhood_size) const;
    int scaleSpaceDoGMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int neighborhood_size, int scale_neighborhood_size) const;

    ImageT make_image() const;
    ImageStackT make_image_stack(int nT) const;
    ScaledImageT make_scaled_image() const;
    ScaledImageStackT make_scaled_image_stack(int nT) const;

    /* Static Methods */
    static void filterLoG(const ImageStackT &im, ImageStackT &fim, const VecT &sigma);
    static void filterDoG(const ImageStackT &im, ImageStackT &fim, const VecT &sigma, FloatT sigma_ratio);
    static void filterGauss(const ImageStackT &im, ImageStackT &fim, const VecT &sigma);
    static void checkMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals);
    static int enumerateImageMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int neighborhood_size);

private:
    int scaleSpaceFrameMaximaRefine(const ScaledImageT &im, IMatT &maxima, VecT &max_vals, int scale_neighborhood_size) const;
    int scaleSpaceFrameMaxima(const ScaledImageT &im, IMatT &maxima, VecT &max_vals, 
                                   int neighborhood_size, int scale_neighborhood_size) const;
    static int combine_maxima(const arma::field<IMatT> &frame_maxima, const arma::field<VecT> &frame_max_vals,
                       IMatT &maxima, VecT &max_vals);
    static void computeDoGSigmas(const MatT &sigma, FloatT sigma_ratio, MatT &gauss_sigmaE, MatT &gauss_sigmaI);
};


template<class FloatT>
inline
typename Boxxer2D<FloatT>::ImageT
Boxxer2D<FloatT>::make_image() const
{
    return ImageT(imsize(0),imsize(1));
}


template<class FloatT>
inline
typename Boxxer2D<FloatT>::ImageStackT
Boxxer2D<FloatT>::make_image_stack(int nT) const
{
    return ImageStackT(imsize(0),imsize(1),nT);
}

template<class FloatT>
inline
typename Boxxer2D<FloatT>::ScaledImageT
Boxxer2D<FloatT>::make_scaled_image() const
{
    return ScaledImageT(imsize(0),imsize(1),nScales);
}

template<class FloatT>
inline
typename Boxxer2D<FloatT>::ScaledImageStackT
Boxxer2D<FloatT>::make_scaled_image_stack(int nT) const
{
    return ScaledImageStackT(imsize(0),imsize(1),nScales,nT);
}


#endif /* _BOXXER2D_H */
