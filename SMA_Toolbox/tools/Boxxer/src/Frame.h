/**
* @file Frame.h
* @author Mark J. Olah (mjo\@cs.unm.edu)
* @date 07-28-2014
* @brief The class declaration for Frames.
*/
#ifndef _BOXXER2D_H
#define _BOXXER2D_H

#include <armadillo>

class Frame
{
    
    int sequence;
    int dim;
    int Npoints;
    arma::umat points;
    arma::vec psfsigma;
    arma::uvec imsize;
    arma::uvec min_boxsize;
    arma::uvec max_boxsize;
    
    
    Frame(arma::umat points);
    
    void getBoxesSize(arma::vec size);
    
}





/*
template<class FloatT>
class Boxxer2D
{
public:
    typedef arma::Col<int> IVecT;
    typedef arma::Mat<int> IMatT;
    typedef arma::Col<FloatT> VecT;
    typedef arma::Mat<FloatT> ImageT;
    typedef arma::Cube<FloatT> ImageStackT;

    const int dim=2;
    IVecT imsize; // [X Y]
    VecT sigma; //[psf_x psf_y]
    Boxxer2D(const IVecT &imsize, const VecT &sigma);

    void boxIt(const ImageStackT &im, IMatT &box_coords, ImageStackT &boxes) const;
    void filterLoG(const ImageStackT &im, ImageStackT &fim) const;
    void filterDoG(const ImageStackT &im, ImageStackT &fim, FloatT sigma_ratio) const;
    void enumerateMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals, int boxsize) const;
    void checkMaxima(const ImageStackT &im, IMatT &maxima, VecT &max_vals) const;
    ImageStackT make_image_stack(int nT) const;
};*/


template<class FloatT>
typename Boxxer2D<FloatT>::ImageStackT
Boxxer2D<FloatT>::make_image_stack(int nT) const
{
    return ImageStackT(imsize(0),imsize(1),nT);
}

#endif /* _BOXXER2D_H */
