/** @file GaussFilter.cpp
 * @author Mark J. Olah (mjo@@cs.unm.edu)
 * @date 07-23-2014
 * @brief
 *
 */

#include <limits>
#include <cassert>
#include <iomanip>
#include <armadillo>
// #include "util.h"
#include "GaussFilter.h"
#include "filter_kernels.h"


/* GaussFIRFilter */

template<class FloatT>
GaussFIRFilter<FloatT>::GaussFIRFilter(int dim_, const IVecT &size_, const FVecT &sigma_)
    : dim(dim_), size(size_), sigma(sigma_), hw(dim_)
{
    assert(1<=dim && dim<=3);
    assert(static_cast<int>(size.n_elem)==dim);
    assert(arma::all(size>0));
    assert(static_cast<int>(sigma.n_elem)==dim);
    assert(arma::all(sigma>0));
}


template<class FloatT>
typename GaussFIRFilter<FloatT>::FVecT 
GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(FloatT sigma, int hw)
{
    int W=hw+1;  //full kernel size is hw*2+1, but we only compute the right half + center pixel
    arma::vec kernel(W);
    double exp_norm=-0.5/(sigma*sigma); // -1/(2*sigma^2)
    kernel(0)=1.;
    double sum=1.;
    for(int r=1; r<W; r++) {
        double val=exp(r*r*exp_norm);
        kernel(r)=val;
        sum+=2*val;
    }
    kernel/=sum;
    return arma::conv_to<FVecT>::from(kernel);
}


template<class FloatT>
typename GaussFIRFilter<FloatT>::FVecT 
GaussFIRFilter<FloatT>::compute_LoG_FIR_kernel(FloatT sigma, int hw)
{
    int W=hw+1;  //full kernel size is hw*2+1, but we only compute the right half + center pixel
    arma::vec kernel(W);
    double sigmanorm=1./(sigma*sigma);
    double norm=sigmanorm/(sqrt(2*arma::Datum<double>::pi)); //1/(sqrt(2*pi)*sigma^3) * sigma <-normalization term
    double exp_norm=-0.5*sigmanorm;
    kernel(0)=norm;
//     double sum=norm;
    for(int r=1; r<W; r++) { 
        double rsq=r*r;
        double val=norm*(1-rsq*sigmanorm)*exp(rsq*exp_norm);
        kernel(r)=val;
//         sum+=2*val;
    }
//     kernel(0)-=sum;
    return arma::conv_to<FVecT>::from(kernel);
}




/* GaussFilter2D */

template<class FloatT>
GaussFilter2D<FloatT>::GaussFilter2D(const IVecT &size, const FVecT &sigma)
    : GaussFIRFilter<FloatT>(2, size, sigma), kernels(2)
{
    auto hw=arma::conv_to<IVecT>::from(arma::ceil(this->default_sigma_hw_ratio * sigma));
    set_kernel_hw(hw);
    temp_im.set_size(size(0),size(1));
}

template<class FloatT>
GaussFilter2D<FloatT>::GaussFilter2D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw)
    : GaussFIRFilter<FloatT>(2, size, sigma), kernels(2)
{
    set_kernel_hw(kernel_hw);
    temp_im.set_size(size(0),size(1));
}


template<class FloatT>
void GaussFilter2D<FloatT>::set_kernel_hw(const IVecT &kernel_half_width)
{
    assert(arma::all(kernel_half_width>0));
    this->hw=kernel_half_width;
    for(int d=0; d<this->dim; d++) 
        kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d), this->hw(d));
}


template<class FloatT>
void GaussFilter2D<FloatT>::filter(const ImageT &im, ImageT &out)
{
    gaussFIR_2Dx<FloatT>(im, temp_im, kernels(0));
    gaussFIR_2Dy<FloatT>(temp_im, out, kernels(1));
}

template<class FloatT>
void GaussFilter2D<FloatT>::test_filter(const ImageT &im)
{
    ImageT fast_out=make_image();
    ImageT slow_out=make_image();
    filter(im, fast_out);
    gaussFIR_2Dx_small<FloatT>(im, temp_im, kernels(0));
    gaussFIR_2Dy_small<FloatT>(temp_im, slow_out, kernels(1));
    FloatT eps=4.*std::numeric_limits<FloatT>::epsilon();
    for(int y=0; y<this->size(1); y++) for(int x=0; x<this->size(0); x++)
        if( fabs(fast_out(x,y)-slow_out(x,y))>eps )
            printf("Fast (%i,%i):%.17f  != Slow (%i,%i):%.17f\n",x,y,fast_out(x,y),x,y,slow_out(x,y));
}

template<class FloatT>
std::ostream& operator<< (std::ostream &out, const GaussFilter2D<FloatT> &filt)
{
    out<<std::setprecision(15);
    auto k0=filt.kernels(0);
    auto k1=filt.kernels(1);
    
    out<<"GaussFilter2D:[size=["<<filt.size(0)<<","<<filt.size(1)<<"]"
        <<" sigma=["<<filt.sigma(0)<<","<<filt.sigma(1)<<"]"
       <<" hw=["<<filt.hw(0)<<","<<filt.hw(1)<<"]"
       <<"\n >>KernelX:(sum:="<<2*arma::sum(k0)-k0(0)<<")\n"<<k0
       <<"\n >>KernelY:(sum:="<<2*arma::sum(k1)-k1(0)<<")\n"<<k1<<"\n";
    out<<std::setprecision(9);
    return out;
}


/* GaussFilter3D */

template<class FloatT>
GaussFilter3D<FloatT>::GaussFilter3D(const IVecT &size, const FVecT &sigma)
    : GaussFIRFilter<FloatT>(3, size, sigma), kernels(3)
{
    auto hw=arma::conv_to<IVecT>::from(arma::ceil(this->default_sigma_hw_ratio * sigma));
    set_kernel_hw(hw);
    temp_im0.set_size(size(0),size(1),size(2));
    temp_im1.set_size(size(0),size(1),size(2));
}

template<class FloatT>
GaussFilter3D<FloatT>::GaussFilter3D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw)
    : GaussFIRFilter<FloatT>(3, size, sigma), kernels(3)
{
    set_kernel_hw(kernel_hw);
    temp_im0.set_size(size(0),size(1),size(2));
    temp_im1.set_size(size(0),size(1),size(2));
}


template<class FloatT>
void GaussFilter3D<FloatT>::set_kernel_hw(const IVecT &kernel_half_width)
{
    assert(arma::all(kernel_half_width>0));
    this->hw=kernel_half_width;
    for(int d=0; d<this->dim; d++) 
        kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d), this->hw(d));
}


template<class FloatT>
void GaussFilter3D<FloatT>::filter(const ImageT &im, ImageT &out)
{
    gaussFIR_3Dx<FloatT>(im, temp_im0, kernels(0));
    gaussFIR_3Dy<FloatT>(temp_im0, temp_im1, kernels(1));
    gaussFIR_3Dz<FloatT>(temp_im1, out, kernels(2));
}

template<class FloatT>
void GaussFilter3D<FloatT>::test_filter(const ImageT &im)
{
    ImageT fast_out=make_image();
    ImageT slow_out=make_image();
    filter(im, fast_out);
    gaussFIR_3Dx_small<FloatT>(im, temp_im0, kernels(0));
    gaussFIR_3Dy_small<FloatT>(temp_im0, temp_im1, kernels(1));
    gaussFIR_3Dz_small<FloatT>(temp_im1, slow_out, kernels(2));
    FloatT eps=4.*std::numeric_limits<FloatT>::epsilon();
    for(int z=0; z<this->size(2); z++) for(int y=0; y<this->size(1); y++) for(int x=0; x<this->size(0); x++)
        if( fabs(fast_out(x,y,z)-slow_out(x,y,z))>eps )
            printf("Fast (%i,%i,%i):%.17f  != Slow (%i,%i,.%i):%.17f\n",x,y,z,fast_out(x,y,z),x,y,z,slow_out(x,y,z));
}

template<class FloatT>
std::ostream& operator<< (std::ostream &out, const GaussFilter3D<FloatT> &filt)
{
    out<<std::setprecision(15);
    auto k0=filt.kernels(0);
    auto k1=filt.kernels(1);
    auto k2=filt.kernels(2);
    
    out<<"GaussFilter3D:[size=["<<filt.size(0)<<","<<filt.size(1)<<","<<filt.size(2)<<"]"
       <<" sigma=["<<filt.sigma(0)<<","<<filt.sigma(1)<<","<<filt.sigma(2)<<"]"
       <<" hw=["<<filt.hw(0)<<","<<filt.hw(1)<<","<<filt.hw(2)<<"]"
       <<"\n >>KernelX:(sum:="<<2*arma::sum(k0)-k0(0)<<")\n"<<k0
       <<"\n >>KernelY:(sum:="<<2*arma::sum(k1)-k1(0)<<")\n"<<k1<<"\n"
       <<"\n >>KernelZ:(sum:="<<2*arma::sum(k2)-k2(0)<<")\n"<<k2<<"\n";
    out<<std::setprecision(9);
    return out;
}

/* DoGFilter2D */
template<class FloatT>
DoGFilter2D<FloatT>::DoGFilter2D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio)
    : GaussFIRFilter<FloatT>(2, size, sigma), sigma_ratio(sigma_ratio), excite_kernels(2), inhibit_kernels(2)
{
    assert(sigma_ratio>1);
    auto hw=arma::conv_to<IVecT>::from(arma::ceil(this->default_sigma_hw_ratio * sigma));
    set_kernel_hw(hw);
    temp_im0.set_size(size(0),size(1));
    temp_im1.set_size(size(0),size(1));
}

template<class FloatT>
DoGFilter2D<FloatT>::DoGFilter2D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio, const IVecT &kernel_hw)
    : GaussFIRFilter<FloatT>(2, size, sigma), sigma_ratio(sigma_ratio), excite_kernels(2), inhibit_kernels(2)
{
    assert(sigma_ratio>1);
    set_kernel_hw(kernel_hw);
    temp_im0.set_size(size(0),size(1));
    temp_im1.set_size(size(0),size(1));
}


template<class FloatT>
void DoGFilter2D<FloatT>::set_kernel_hw(const IVecT &kernel_half_width)
{
    assert(arma::all(kernel_half_width>0));
    this->hw=kernel_half_width;
    for(int d=0; d<this->dim; d++) {
        excite_kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d), this->hw(d));
        inhibit_kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d)*sigma_ratio, this->hw(d));
    }
}

template<class FloatT>
void DoGFilter2D<FloatT>::set_sigma_ratio(FloatT _sigma_ratio)
{
    assert(_sigma_ratio>1);
    sigma_ratio = _sigma_ratio;
    set_kernel_hw(this->hw);
}

template<class FloatT>
void DoGFilter2D<FloatT>::filter(const ImageT &im, ImageT &out)
{
    gaussFIR_2Dx<FloatT>(im, temp_im0, excite_kernels(0));
    gaussFIR_2Dy<FloatT>(temp_im0, out, excite_kernels(1));
    
    gaussFIR_2Dx<FloatT>(im, temp_im1, inhibit_kernels(0));
    gaussFIR_2Dy<FloatT>(temp_im1, temp_im0, inhibit_kernels(1));
    out-=temp_im0;
}

template<class FloatT>
void DoGFilter2D<FloatT>::test_filter(const ImageT &im)
{
    ImageT fast_out=make_image();
    ImageT slow_out=make_image();
    filter(im, fast_out);
    gaussFIR_2Dx_small<FloatT>(im, temp_im0, excite_kernels(0));
    gaussFIR_2Dy_small<FloatT>(temp_im0, slow_out, excite_kernels(1));
    
    gaussFIR_2Dx_small<FloatT>(im, temp_im1, inhibit_kernels(0));
    gaussFIR_2Dy_small<FloatT>(temp_im1, temp_im0, inhibit_kernels(1));
    slow_out-=temp_im0;
    FloatT eps=4.*std::numeric_limits<FloatT>::epsilon();
    for(int y=0; y<this->size(1); y++) for(int x=0; x<this->size(0); x++)
        if( fabs(fast_out(x,y)-slow_out(x,y))>eps )
            printf("Fast (%i,%i):%.17f  != Slow (%i,%i):%.17f\n",x,y,fast_out(x,y),x,y,slow_out(x,y));
}

/* DoGFilter3D */
template<class FloatT>
DoGFilter3D<FloatT>::DoGFilter3D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio)
    : GaussFIRFilter<FloatT>(3, size, sigma), sigma_ratio(sigma_ratio), excite_kernels(3), inhibit_kernels(3)
{
    assert(sigma_ratio>1);
    auto hw=arma::conv_to<IVecT>::from(arma::ceil(this->default_sigma_hw_ratio * sigma));
    set_kernel_hw(hw);
    temp_im0.set_size(size(0),size(1),size(2));
    temp_im1.set_size(size(0),size(1),size(2));
}

template<class FloatT>
DoGFilter3D<FloatT>::DoGFilter3D(const IVecT &size, const FVecT &sigma, FloatT sigma_ratio, const IVecT &kernel_hw)
    : GaussFIRFilter<FloatT>(3, size, sigma), sigma_ratio(sigma_ratio), excite_kernels(3), inhibit_kernels(3)
{
    assert(sigma_ratio>1);
    set_kernel_hw(kernel_hw);
    temp_im0.set_size(size(0),size(1),size(2));
    temp_im1.set_size(size(0),size(1),size(2));
}


template<class FloatT>
void DoGFilter3D<FloatT>::set_kernel_hw(const IVecT &kernel_half_width)
{
    assert(arma::all(kernel_half_width>0));
    this->hw=kernel_half_width;
    for(int d=0; d<this->dim; d++) {
        excite_kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d), this->hw(d));
        inhibit_kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d)*sigma_ratio, this->hw(d));
    }
}

template<class FloatT>
void DoGFilter3D<FloatT>::set_sigma_ratio(FloatT _sigma_ratio)
{
    assert(_sigma_ratio>1);
    sigma_ratio = _sigma_ratio;
    set_kernel_hw(this->hw);
}

template<class FloatT>
void DoGFilter3D<FloatT>::filter(const ImageT &im, ImageT &out)
{
    gaussFIR_3Dx<FloatT>(im, temp_im0, excite_kernels(0));
    gaussFIR_3Dy<FloatT>(temp_im0, temp_im1, excite_kernels(1));
    gaussFIR_3Dz<FloatT>(temp_im1, out, excite_kernels(2));
    
    gaussFIR_3Dx<FloatT>(im, temp_im0, inhibit_kernels(0));
    gaussFIR_3Dy<FloatT>(temp_im0, temp_im1, inhibit_kernels(1));
    gaussFIR_3Dz<FloatT>(temp_im1, temp_im0, inhibit_kernels(2));
    out-=temp_im0;
}

template<class FloatT>
void DoGFilter3D<FloatT>::test_filter(const ImageT &im)
{
    ImageT fast_out=make_image();
    ImageT slow_out=make_image();
    filter(im, fast_out);    
    gaussFIR_3Dx_small<FloatT>(im, temp_im0, excite_kernels(0));
    gaussFIR_3Dy_small<FloatT>(temp_im0, temp_im1, excite_kernels(1));
    gaussFIR_3Dz_small<FloatT>(temp_im1, slow_out, excite_kernels(2));
    
    gaussFIR_3Dx_small<FloatT>(im, temp_im0, inhibit_kernels(0));
    gaussFIR_3Dy_small<FloatT>(temp_im0, temp_im1, inhibit_kernels(1));
    gaussFIR_3Dz_small<FloatT>(temp_im1, temp_im0, inhibit_kernels(2));
    slow_out-=temp_im0;
    FloatT eps=4.*std::numeric_limits<FloatT>::epsilon();
    for(int z=0; z<this->size(2); z++) for(int y=0; y<this->size(1); y++) for(int x=0; x<this->size(0); x++)
        if( fabs(fast_out(x,y,z)-slow_out(x,y,z))>eps )
            printf("Fast (%i,%i,%i):%.17f  != Slow (%i,%i,%i):%.17f\n",x,y,z,fast_out(x,y,z),x,y,z,slow_out(x,y,z));
}


/* LoGFilter2D */

template<class FloatT>
LoGFilter2D<FloatT>::LoGFilter2D(const IVecT &size, const FVecT &sigma)
    : GaussFIRFilter<FloatT>(2, size, sigma), gauss_kernels(2), LoG_kernels(2)
{
    auto hw=arma::conv_to<IVecT>::from(arma::ceil(this->default_sigma_hw_ratio * sigma));
    set_kernel_hw(hw);
    temp_im0.set_size(size(0),size(1));
    temp_im1.set_size(size(0),size(1));
}

template<class FloatT>
LoGFilter2D<FloatT>::LoGFilter2D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw)
    : GaussFIRFilter<FloatT>(2, size, sigma), gauss_kernels(2), LoG_kernels(2)
{
    set_kernel_hw(kernel_hw);
    temp_im0.set_size(size(0),size(1));
    temp_im1.set_size(size(0),size(1));
}


template<class FloatT>
void LoGFilter2D<FloatT>::set_kernel_hw(const IVecT &kernel_half_width)
{
    assert(arma::all(kernel_half_width>0));
    this->hw=kernel_half_width;
    for(int d=0; d<this->dim; d++) {
        gauss_kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d), this->hw(d));
        LoG_kernels(d)=GaussFIRFilter<FloatT>::compute_LoG_FIR_kernel(this->sigma(d), this->hw(d));
    }
}



template<class FloatT>
void LoGFilter2D<FloatT>::filter(const ImageT &im, ImageT &out)
{
    gaussFIR_2Dy<FloatT>(im, temp_im0, LoG_kernels(1)); //G''(y)fc
    gaussFIR_2Dx<FloatT>(temp_im0, out, gauss_kernels(0)); //G(x)

    gaussFIR_2Dy<FloatT>(im, temp_im0, gauss_kernels(1)); //G(y)
    gaussFIR_2Dx<FloatT>(temp_im0, temp_im1, LoG_kernels(0)); //G''(x)
    out+=temp_im1;

//     gaussFIR_2Dy<FloatT>(im, temp_im0, LoG_kernels(1)); //G''(y)fc
//     gaussFIR_2Dx<FloatT>(im, temp_im1, gauss_kernels(0)); //G(x)
//     temp_im0 = temp_im0 % temp_im1; //straight product
//
//     gaussFIR_2Dy<FloatT>(im, temp_im1, gauss_kernels(1)); //G(y)
//     gaussFIR_2Dx<FloatT>(im, out, LoG_kernels(0)); //G''(x)
//     out = (out%temp_im1)+temp_im0;
}

template<class FloatT>
void LoGFilter2D<FloatT>::test_filter(const ImageT &im)
{
    ImageT fast_out=make_image();
    ImageT slow_out=make_image();
    filter(im, fast_out);
    gaussFIR_2Dy_small<FloatT>(im, temp_im0, LoG_kernels(1)); //G''(y)fc
    gaussFIR_2Dx_small<FloatT>(temp_im0, slow_out, gauss_kernels(0)); //G(x)

    gaussFIR_2Dy_small<FloatT>(im, temp_im0, gauss_kernels(1)); //G(y)
    gaussFIR_2Dx_small<FloatT>(temp_im0, temp_im1, LoG_kernels(0)); //G''(x)
    slow_out+=temp_im1;

    FloatT eps=4.*std::numeric_limits<FloatT>::epsilon();
    for(int y=0; y<this->size(1); y++) for(int x=0; x<this->size(0); x++)
        if( fabs(fast_out(x,y)-slow_out(x,y))>eps )
            printf("Fast (%i,%i):%.17f  != Slow (%i,%i):%.17f\n",x,y,fast_out(x,y),x,y,slow_out(x,y));
}

template<class FloatT>
std::ostream& operator<< (std::ostream &out, const LoGFilter2D<FloatT> &filt)
{
    out<<std::setprecision(15);
    auto gk0=filt.gauss_kernels(0);
    auto gk1=filt.gauss_kernels(1);
    auto logk0=filt.LoG_kernels(0);
    auto logk1=filt.LoG_kernels(1);
    
    out<<"LoGFilter2D:[size=["<<filt.size(0)<<","<<filt.size(1)<<"]"
        <<" sigma=["<<filt.sigma(0)<<","<<filt.sigma(1)<<"]"
       <<" hw=["<<filt.hw(0)<<","<<filt.hw(1)<<"]"
       <<"\n >>GaussKernelX:(sum:="<<2*arma::sum(gk0)-gk0(0)<<")\n"<<gk0
       <<"\n >>GaussKernelY:(sum:="<<2*arma::sum(gk1)-gk1(0)<<")\n"<<gk1
       <<"\n >>LoGKernelX:(sum:="<<2*arma::sum(logk0)-logk0(0)<<")\n"<<logk0
       <<"\n >>LoGKernelY:(sum:="<<2*arma::sum(logk1)-logk1(0)<<")\n"<<logk1<<"\n";
    out<<std::setprecision(9);
    return out;
}

/* LoGFilter3D */

template<class FloatT>
LoGFilter3D<FloatT>::LoGFilter3D(const IVecT &size, const FVecT &sigma)
    : GaussFIRFilter<FloatT>(3, size, sigma), gauss_kernels(3), LoG_kernels(3)
{
    auto hw=arma::conv_to<IVecT>::from(arma::ceil(this->default_sigma_hw_ratio * sigma));
    set_kernel_hw(hw);
    temp_im0.set_size(size(0),size(1),size(2));
    temp_im1.set_size(size(0),size(1),size(2));
//     temp_im2.set_size(size(0),size(1),size(2));
}

template<class FloatT>
LoGFilter3D<FloatT>::LoGFilter3D(const IVecT &size, const FVecT &sigma, const IVecT &kernel_hw)
    : GaussFIRFilter<FloatT>(3, size, sigma), gauss_kernels(3), LoG_kernels(3)
{
    set_kernel_hw(kernel_hw);
    temp_im0.set_size(size(0),size(1),size(2));
    temp_im1.set_size(size(0),size(1),size(2));
//     temp_im2.set_size(size(0),size(1),size(2));
}


template<class FloatT>
void LoGFilter3D<FloatT>::set_kernel_hw(const IVecT &kernel_half_width)
{
    assert(arma::all(kernel_half_width>0));
    this->hw=kernel_half_width;
    for(int d=0; d<this->dim; d++) {
        gauss_kernels(d)=GaussFIRFilter<FloatT>::compute_Gauss_FIR_kernel(this->sigma(d), this->hw(d));
        LoG_kernels(d)=GaussFIRFilter<FloatT>::compute_LoG_FIR_kernel(this->sigma(d), this->hw(d));
    }
}


template<class FloatT>
void LoGFilter3D<FloatT>::filter(const ImageT &im, ImageT &out)
{
    gaussFIR_3Dz<FloatT>(im, temp_im0, gauss_kernels(2));
    gaussFIR_3Dy<FloatT>(temp_im0, temp_im1, gauss_kernels(1));
    gaussFIR_3Dx<FloatT>(temp_im1, out, LoG_kernels(0));

    gaussFIR_3Dz<FloatT>(im, temp_im0, gauss_kernels(2));
    gaussFIR_3Dy<FloatT>(temp_im0, temp_im1, LoG_kernels(1));
    gaussFIR_3Dx<FloatT>(temp_im1, temp_im0, gauss_kernels(0));
    out+=temp_im0;

    gaussFIR_3Dz<FloatT>(im, temp_im0, LoG_kernels(2));
    gaussFIR_3Dy<FloatT>(temp_im0, temp_im1, gauss_kernels(1));
    gaussFIR_3Dx<FloatT>(temp_im1, temp_im0, gauss_kernels(0));
    out+=temp_im0;
//     gaussFIR_3Dz<FloatT>(im, temp_im0, gauss_kernels(2));
//     gaussFIR_3Dy<FloatT>(im, temp_im1, gauss_kernels(1));
//     gaussFIR_3Dx<FloatT>(im, out, LoG_kernels(0));
//     out = out%temp_im0%temp_im1;
//
//     gaussFIR_3Dz<FloatT>(im, temp_im0, gauss_kernels(2));
//     gaussFIR_3Dy<FloatT>(im, temp_im1, LoG_kernels(1));
//     temp_im0 = temp_im0%temp_im1;
//     gaussFIR_3Dx<FloatT>(im, temp_im1, gauss_kernels(0));
//     out += temp_im0%temp_im1;
//
//     gaussFIR_3Dz<FloatT>(im, temp_im0, LoG_kernels(2));
//     gaussFIR_3Dy<FloatT>(im, temp_im1, gauss_kernels(1));
//     temp_im0 = temp_im0%temp_im1;
//     gaussFIR_3Dx<FloatT>(im, temp_im0, gauss_kernels(0));
//     out += temp_im0%temp_im1;
}

template<class FloatT>
void LoGFilter3D<FloatT>::test_filter(const ImageT &im)
{
    ImageT fast_out=make_image();
    ImageT slow_out=make_image();
    filter(im, fast_out);
    gaussFIR_3Dz_small<FloatT>(im, temp_im0, gauss_kernels(2));
    gaussFIR_3Dy_small<FloatT>(temp_im0, temp_im1, gauss_kernels(1));
    gaussFIR_3Dx_small<FloatT>(temp_im1, slow_out, LoG_kernels(0));

    gaussFIR_3Dz_small<FloatT>(im, temp_im0, gauss_kernels(2));
    gaussFIR_3Dy_small<FloatT>(temp_im0, temp_im1, LoG_kernels(1));
    gaussFIR_3Dx_small<FloatT>(temp_im1, temp_im0, gauss_kernels(0));
    slow_out+=temp_im0;

    gaussFIR_3Dz_small<FloatT>(im, temp_im0, LoG_kernels(2));
    gaussFIR_3Dy_small<FloatT>(temp_im0, temp_im1, gauss_kernels(1));
    gaussFIR_3Dx_small<FloatT>(temp_im1, temp_im0, gauss_kernels(0));
    slow_out+=temp_im0;
    FloatT eps=4.*std::numeric_limits<FloatT>::epsilon();
    for(int z=0; z<this->size(2); z++) for(int y=0; y<this->size(1); y++) for(int x=0; x<this->size(0); x++)
        if( fabs(fast_out(x,y,z)-slow_out(x,y,z))>eps )
            printf("Fast (%i,%i,%i):%.17f  != Slow (%i,%i,.%i):%.17f\n",x,y,z,fast_out(x,y,z),x,y,z,slow_out(x,y,z));
}

template<class FloatT>
std::ostream& operator<< (std::ostream &out, const LoGFilter3D<FloatT> &filt)
{
    out<<std::setprecision(15);
    auto gk0=filt.gauss_kernels(0);
    auto gk1=filt.gauss_kernels(1);
    auto gk2=filt.gauss_kernels(2);
    auto logk0=filt.LoG_kernels(0);
    auto logk1=filt.LoG_kernels(1);
    auto logk2=filt.LoG_kernels(2);
    
    out<<"LoGFilter3D:[size=["<<filt.size(0)<<","<<filt.size(1)<<","<<filt.size(2)<<"]"
       <<" sigma=["<<filt.sigma(0)<<","<<filt.sigma(1)<<","<<filt.sigma(2)<<"]"
       <<" hw=["<<filt.hw(0)<<","<<filt.hw(1)<<","<<filt.hw(2)<<"]"
       <<"\n >>Gauss KernelX:(sum:="<<2*arma::sum(gk0)-gk0(0)<<")\n"<<gk0
       <<"\n >>Gauss KernelY:(sum:="<<2*arma::sum(gk1)-gk1(0)<<")\n"<<gk1<<"\n"
       <<"\n >>Gauss KernelZ:(sum:="<<2*arma::sum(gk2)-gk2(0)<<")\n"<<gk2<<"\n"
       <<"\n >>LoG KernelX:(sum:="<<2*arma::sum(logk0)-logk0(0)<<")\n"<<logk0<<"\n"
       <<"\n >>LoG KernelY:(sum:="<<2*arma::sum(logk1)-logk1(0)<<")\n"<<logk1<<"\n"
       <<"\n >>LoG KernelZ:(sum:="<<2*arma::sum(logk2)-logk2(0)<<")\n"<<logk2<<"\n";
    out<<std::setprecision(9);
    return out;
}



/* Explicit Template Instantiation */
template class GaussFIRFilter<float>;
template class GaussFIRFilter<double>;

template class GaussFilter2D<float>;
template class GaussFilter2D<double>;

template class GaussFilter3D<float>;
template class GaussFilter3D<double>;

template class DoGFilter2D<float>;
template class DoGFilter2D<double>;

template class DoGFilter3D<float>;
template class DoGFilter3D<double>;

template class LoGFilter2D<float>;
template class LoGFilter2D<double>;

template class LoGFilter3D<float>;
template class LoGFilter3D<double>;


template std::ostream& operator<< <float>(std::ostream &out, const GaussFilter2D<float> &filt);
template std::ostream& operator<< <double>(std::ostream &out, const GaussFilter2D<double> &filt);

template std::ostream& operator<< <float>(std::ostream &out, const GaussFilter3D<float> &filt);
template std::ostream& operator<< <double>(std::ostream &out, const GaussFilter3D<double> &filt);

template std::ostream& operator<< <float>(std::ostream &out, const LoGFilter2D<float> &filt);
template std::ostream& operator<< <double>(std::ostream &out, const LoGFilter2D<double> &filt);

template std::ostream& operator<< <float>(std::ostream &out, const LoGFilter3D<float> &filt);
template std::ostream& operator<< <double>(std::ostream &out, const LoGFilter3D<double> &filt);
