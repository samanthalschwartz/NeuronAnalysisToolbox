/**
 * @file Boxxer3D.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-28-2014
 * @brief The class declaration for Boxxer3D.
 */
#ifndef _BOXXER3D_H
#define _BOXXER3D_H

#include <armadillo>
#include "hypercube.h"
#include "GaussFilter.h"

/**
 * @breif A box finding algroithm for 3D hyperspectral microscopy data.  Which estimates the center coordinates
 *  of gaussian blobs with anisotropic sigmas.
 *
 * All image data manipulated is stored as column-major FloatT arrays with dimension ordering [L Y X T].
 *
 * The Boxxer3D class makes uses of lower level class which are agnostic about the data source being hyperspectral,
 * they don't care what the coordinate dimensions represent scientfically, but this class is associated with the
 * Matlab Boxxer3D class and so maintains the knowledge that the actual coordinates are [L Y X T].
 */
template<class FloatT>
class Boxxer3D
{
public:
    typedef arma::Col<int> IVecT;
    typedef arma::Mat<int> IMatT;
    typedef arma::Col<FloatT> VecT;
    typedef arma::Mat<FloatT> MatT;
    typedef arma::Cube<FloatT> ImageT;
    typedef Hypercube<FloatT> ImageStackT;
    typedef Hypercube<FloatT> ScaledImageT;

    static const FloatT DefaultSigmaRatio;
    static const int dim;

    int nScales;
    IVecT imsize; // Size of each dimension for the column-major data. HSData is [L Y X] this is [row, col, slice]
    MatT sigma; // sized: [2 x nScales].  Rows are [psf_L, psf_y, psf_x] cols are the different scales 
                //CRITICAL: the order of sigma rows must match the order of dimension in imsize.
    FloatT sigma_ratio;
    Boxxer3D(const IVecT &size, const MatT &sigma);
    
    void setDoGSigmaRatio(FloatT sigma_ratio);

    void filterScaledLoG(const ImageT &im, ScaledImageT &fim);
    void filterScaledDoG(const ImageT &im, ScaledImageT &fim);
    int scaleSpaceLoGMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int neighborhood_size, int scale_neighborhood_size);
    int scaleSpaceDoGMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int neighborhood_size, int scale_neighborhood_size);

    ImageT make_image() const;
    ImageStackT make_image_stack(int nT) const;
    ScaledImageT make_scaled_image() const;

    /* Static Methods */
    static void filterLoG(const ImageStackT &im, ImageStackT &fim, const VecT &sigma);
    static void filterDoG(const ImageStackT &im, ImageStackT &fim, const VecT &sigma, FloatT sigma_ratio);
    static void filterGauss(const ImageStackT &im, ImageStackT &fim, const VecT &sigma);
    static void checkMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals);
    static int enumerateImageMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int neighborhood_size);

private:
//     std::vector<std::vector<LoGFilter3D<FloatT>>> log_scale_filters; //rows are scales columns are thread#
//     std::vector<std::vector<DoGFilter3D<FloatT>>> dog_scale_filters; //rows are scales columns are thread#
    int scaleSpaceFrameMaximaRefine(const ScaledImageT &im, IMatT &maxima, VecT &max_vals, int scale_neighborhood_size) const;
    int scaleSpaceFrameMaxima(const ScaledImageT &im, IMatT &maxima, VecT &max_vals, 
                              int neighborhood_size, int scale_neighborhood_size) const;
    static int combine_maxima(const arma::field<IMatT> &frame_maxima, const arma::field<VecT> &frame_max_vals,
                              IMatT &maxima, VecT &max_vals);
    void initialize_log_scale_filters();
    void initialize_dog_scale_filters();
};

template<class FloatT>
typename Boxxer3D<FloatT>::ImageT
Boxxer3D<FloatT>::make_image() const
{
    return ImageT(imsize(0),imsize(1),imsize(2));
}

template<class FloatT>
typename Boxxer3D<FloatT>::ImageStackT
Boxxer3D<FloatT>::make_image_stack(int nT) const
{
    return ImageStackT(imsize(0),imsize(1),imsize(2),nT);
}

template<class FloatT>
typename Boxxer3D<FloatT>::ScaledImageT
Boxxer3D<FloatT>::make_scaled_image() const
{
    return ScaledImageT(imsize(0),imsize(1),imsize(2),nScales);
}

#endif /* _BOXXER3D_H */
