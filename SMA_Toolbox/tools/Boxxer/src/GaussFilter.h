/** @file GaussFilter.h
 * @author Mark J. Olah (mjo@@cs.unm.edu)
 * @date 07-23-2014
 * @brief 
 * 
 * These classes are meant to be a per-thread worker class or a direct interface
 * for single threaded processes.  Each object has its own local storage of 1 or
 * 2 frames for the transforms.
 */
#ifndef _GAUSSFILTER_H
#define _GAUSSFILTER_H

#include <vector>
#include <memory>
#include <armadillo>

template<class FloatT>
class GaussFIRFilter
{
public:
    typedef arma::Col<int>    IVecT;
    typedef arma::Col<FloatT> FVecT;
    typedef arma::Mat<FloatT> FMatT;
    
    int dim;
    IVecT size; //[nrows, ncols]
    FVecT sigma; //sigma to filter at [sigma_rows, sigma_cols]
    IVecT hw; //Half width Full kernel width is 2*hw+1. Recommended 3sigma

    GaussFIRFilter(int dim, const IVecT &size, const FVecT &sigma);

    virtual void set_kernel_hw(const IVecT &kernel_half_width)=0;

    static FVecT compute_Gauss_FIR_kernel(FloatT sigma, int hw);
    static FVecT compute_LoG_FIR_kernel(FloatT sigma, int hw);
protected:
    const int max_kernel_hw=30;
    const FloatT default_sigma_hw_ratio=3.; //hw will be 3sigma if not specified.
};

template<class FloatT>
class GaussFilter2D : public GaussFIRFilter<FloatT>
{
public:
    typedef typename GaussFIRFilter<FloatT>::IVecT IVecT;
    typedef typename GaussFIRFilter<FloatT>::FVecT FVecT;
    typedef typename GaussFIRFilter<FloatT>::FMatT FMatT;
    typedef arma::Mat<FloatT> ImageT;

    GaussFilter2D(const IVecT &size, const FVecT &sigma);
    GaussFilter2D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw);
    void set_kernel_hw(const IVecT &kernel_half_width);
    ImageT make_image() const;

    void filter(const ImageT &im, ImageT &out);
    void test_filter(const ImageT &im);

    template <class T>
    friend std::ostream& operator<<(std::ostream &out, const GaussFilter2D<T> &filt);
private:
    ImageT temp_im;
    arma::field<FVecT> kernels;
};

template<class FloatT>
class DoGFilter2D : public GaussFIRFilter<FloatT>
{
public:
    typedef typename GaussFIRFilter<FloatT>::IVecT IVecT;
    typedef typename GaussFIRFilter<FloatT>::FVecT FVecT;
    typedef typename GaussFIRFilter<FloatT>::FMatT FMatT;
    typedef arma::Mat<FloatT> ImageT;

    FloatT sigma_ratio;
    
    DoGFilter2D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio);
    DoGFilter2D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio, const IVecT &kernel_hw);
    void set_kernel_hw(const IVecT &kernel_half_width);
    void set_sigma_ratio(FloatT sigma_ratio);
    ImageT make_image() const;

    void filter(const ImageT &im, ImageT &out);
    void test_filter(const ImageT &im);

    template <class T>
    friend std::ostream& operator<<(std::ostream &out, const DoGFilter2D<T> &filt);
private:
    ImageT temp_im0;
    ImageT temp_im1;
    arma::field<FVecT> excite_kernels;
    arma::field<FVecT> inhibit_kernels;
};


template<class FloatT>
class GaussFilter3D : public GaussFIRFilter<FloatT>
{
public:
    typedef typename GaussFIRFilter<FloatT>::IVecT IVecT;
    typedef typename GaussFIRFilter<FloatT>::FVecT FVecT;
    typedef typename GaussFIRFilter<FloatT>::FMatT FMatT;
    typedef arma::Cube<FloatT> ImageT;


    GaussFilter3D(const IVecT &size, const FVecT &sigma);
    GaussFilter3D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw);
    void set_kernel_hw(const IVecT &kernel_half_width);
    ImageT make_image() const;

    void filter(const ImageT &im, ImageT &out);
    void test_filter(const ImageT &im);

    template <class T>
    friend std::ostream& operator<<(std::ostream &out, const GaussFilter3D<T> &filt);
private:
    ImageT temp_im0;
    ImageT temp_im1;
    arma::field<FVecT> kernels;
};

template<class FloatT>
class DoGFilter3D : public GaussFIRFilter<FloatT>
{
public:
    typedef typename GaussFIRFilter<FloatT>::IVecT IVecT;
    typedef typename GaussFIRFilter<FloatT>::FVecT FVecT;
    typedef typename GaussFIRFilter<FloatT>::FMatT FMatT;
    typedef arma::Cube<FloatT> ImageT;

    FloatT sigma_ratio;

    DoGFilter3D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio);
    DoGFilter3D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio, const IVecT &kernel_hw);
    void set_kernel_hw(const IVecT &kernel_half_width);
    void set_sigma_ratio(FloatT sigma_ratio);
    ImageT make_image() const;

    void filter(const ImageT &im, ImageT &out);
    void test_filter(const ImageT &im);

    template <class T>
    friend std::ostream& operator<<(std::ostream &out, const GaussFilter3D<T> &filt);
private:
    ImageT temp_im0;
    ImageT temp_im1;
    arma::field<FVecT> excite_kernels;
    arma::field<FVecT> inhibit_kernels;
};

template<class FloatT>
class LoGFilter2D : public GaussFIRFilter<FloatT>
{
public:
    typedef typename GaussFIRFilter<FloatT>::IVecT IVecT;
    typedef typename GaussFIRFilter<FloatT>::FVecT FVecT;
    typedef typename GaussFIRFilter<FloatT>::FMatT FMatT;
    typedef arma::Mat<FloatT> ImageT;

    LoGFilter2D(const IVecT &size, const FVecT &sigma);
    LoGFilter2D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw);
    void set_kernel_hw(const IVecT &kernel_half_width);
    ImageT make_image() const;

    void filter(const ImageT &im, ImageT &out);
    void test_filter(const ImageT &im);

    template <class T>
    friend std::ostream& operator<<(std::ostream &out, const LoGFilter2D<T> &filt);
private:
    ImageT temp_im0;
    ImageT temp_im1;
    arma::field<FVecT>  gauss_kernels; 
    arma::field<FVecT>  LoG_kernels; 
};

template<class FloatT>
class LoGFilter3D : public GaussFIRFilter<FloatT>
{
public:
    typedef typename GaussFIRFilter<FloatT>::IVecT IVecT;
    typedef typename GaussFIRFilter<FloatT>::FVecT FVecT;
    typedef typename GaussFIRFilter<FloatT>::FMatT FMatT;
    typedef arma::Cube<FloatT> ImageT;

    LoGFilter3D(const IVecT &size, const FVecT &sigma);
    LoGFilter3D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw);
    void set_kernel_hw(const IVecT &kernel_half_width);
    ImageT make_image() const;

    void filter(const ImageT &im, ImageT &out);
    void test_filter(const ImageT &im);

    template <class T>
    friend std::ostream& operator<<(std::ostream &out, const LoGFilter3D<T> &filt);
private:
    ImageT temp_im0, temp_im1;
    arma::field<FVecT> gauss_kernels; 
    arma::field<FVecT> LoG_kernels; 
};

/* Inlined Methods */

template<class FloatT>
inline
typename GaussFilter2D<FloatT>::ImageT
GaussFilter2D<FloatT>::make_image() const
{
    return ImageT(this->size(0),this->size(1));
}


template<class FloatT>
inline
typename GaussFilter3D<FloatT>::ImageT
GaussFilter3D<FloatT>::make_image() const
{
    return ImageT(this->size(0),this->size(1),this->size(2));
}

template<class FloatT>
inline
typename DoGFilter2D<FloatT>::ImageT
DoGFilter2D<FloatT>::make_image() const
{
    return ImageT(this->size(0),this->size(1));
}

template<class FloatT>
inline
typename DoGFilter3D<FloatT>::ImageT
DoGFilter3D<FloatT>::make_image() const
{
    return ImageT(this->size(0),this->size(1),this->size(2));
}


template<class FloatT>
inline
typename LoGFilter2D<FloatT>::ImageT
LoGFilter2D<FloatT>::make_image() const
{
    return ImageT(this->size(0),this->size(1));
}

template<class FloatT>
inline
typename LoGFilter3D<FloatT>::ImageT
LoGFilter3D<FloatT>::make_image() const
{
    return ImageT(this->size(0),this->size(1),this->size(2));
}

#endif /* _GAUSSFILTER_H */
