/** @file Gauss2DModel.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-13-2014
 * @brief The class declaration and inline and templated functions for Gauss2DModel.
 */

#ifndef _GAUSS2DMODEL_H
#define _GAUSS2DMODEL_H

#include "PointEmitter2DModel.h"

/** @brief A base class for likelihood models for point emitters imaged in 2D with symmmetric PSF.
 *
 *
 *
 */
class Gauss2DModel : public PointEmitter2DModel {
public:
    /* Model matrix and vector types */
    typedef arma::vec::fixed<4> ParamT; /**< A type for the set of parameters estimated by the model */
    typedef arma::mat::fixed<4,4> ParamMatT; /**< A matrix type for the Hessian used by the CRLB estimation */
    static const std::vector<std::string> param_names;

    class Stencil {
    public:
        bool derivatives_computed=false;
        typedef Gauss2DModel::ParamT ParamT;
        Gauss2DModel const *model;
        ParamT theta;
        VecT dx, dy;
        VecT Gx, Gy;
        VecT X, Y;
        VecT DX, DY;
        VecT DXS, DYS;
        Stencil() {}
        Stencil(const Gauss2DModel &model, const ParamT &theta, bool compute_derivatives=true);
        void compute_derivatives();
        inline double x() const {return theta(0);}
        inline double y() const {return theta(1);}
        inline double I() const {return theta(2);}
        inline double bg() const {return theta(3);}
        friend std::ostream& operator<<(std::ostream &out, const Gauss2DModel::Stencil &s);
    };

    Gauss2DModel(const IVecT &size, const VecT &psf_sigma);

    /* Make arrays for working with model data */
    ParamT make_param() const;
    ParamT make_param(double x, double y, double I, double bg) const;
    ParamT make_param(const ParamT &theta) const;
    ParamMatT make_param_mat() const;
    Stencil make_stencil(const ParamT &theta, bool compute_derivatives=true) const;
    Stencil make_stencil(double x, double y, double I=1.0, double bg=0.0, bool compute_derivatives=true) const;


    /* Model Pixel Value And Derivatives */
    double model_value(int i, int j, const Stencil &s) const;
    void pixel_grad(int i, int j, const Stencil &s, ParamT &pgrad) const;
    void pixel_grad2(int i, int j, const Stencil &s, ParamT &pgrad2) const;
    void pixel_hess(int i, int j, const Stencil &s, ParamMatT &hess) const;
    void pixel_hess_update(int i, int j, const Stencil &s, double dm_ratio_m1, 
                           double dmm_ratio, ParamT &grad, ParamMatT &hess) const;
 
    virtual void bound_theta(ParamT &theta) const=0;
                           
    /* Compute the Log likelihood of an image at theta */
    Stencil initial_theta_estimate(const ImageT &im, const ParamT &theta_init) const;

    /* Posterior Sampling */
    void sample_candidate_theta(int sample_index, RNG &rng, ParamT &canidate_theta, double scale=1.0) const;

protected:
    VecT gaussian_Xstencil; /**< A stencil for gaussian filters with this size and psf*/
    VecT gaussian_Ystencil; /**< A stencil for gaussian filters with this size and psf*/
};

/* Function Declarations */



/* Inline Method Definitions */

inline
Gauss2DModel::ParamT
Gauss2DModel::make_param() const
{
    return ParamT();
}

inline
Gauss2DModel::ParamT
Gauss2DModel::make_param(double x, double y, double I, double bg) const
{
    ParamT theta;
    theta<<x<<y<<I<<bg;
    bound_theta(theta);
    return theta;
}

inline
Gauss2DModel::ParamT
Gauss2DModel::make_param(const ParamT &theta) const
{
    ParamT ntheta(theta);
    bound_theta(ntheta);
    return ntheta;
}


inline
Gauss2DModel::ParamMatT
Gauss2DModel::make_param_mat() const
{
    return ParamMatT();
}

inline
Gauss2DModel::Stencil
Gauss2DModel::make_stencil(const ParamT &theta, bool compute_derivatives) const
{
    return Stencil(*this,make_param(theta),compute_derivatives);
}

inline
Gauss2DModel::Stencil
Gauss2DModel::make_stencil(double x, double y, double I, double bg, bool compute_derivatives) const
{
    return Stencil(*this,make_param(x,y,I,bg),compute_derivatives);
}


inline
double Gauss2DModel::model_value(int i, int j, const Stencil &s) const
{
    return s.bg()+s.I()*s.X(i)*s.Y(j);
}

inline
void
Gauss2DModel::pixel_grad(int i, int j, const Stencil &s, ParamT &pgrad) const
{
    double I=s.I();
    pgrad(0) = I * s.DX(i) * s.Y(j);
    pgrad(1) = I * s.DY(j) * s.X(i);
    pgrad(2) = s.X(i) * s.Y(j);
    pgrad(3) = 1.;
}

inline
void
Gauss2DModel::pixel_grad2(int i, int j, const Stencil &s, ParamT &pgrad2) const
{
    double I=s.I();
    pgrad2(0)= I/psf_sigma(0) * s.DXS(i) * s.Y(j);
    pgrad2(1)= I/psf_sigma(1) * s.DYS(j) * s.X(i);
    pgrad2(2)= 0;
    pgrad2(3)= 0;
}

inline
void
Gauss2DModel::pixel_hess(int i, int j, const Stencil &s, ParamMatT &hess) const
{
    hess.zeros();
    double I=s.I();
    hess(0,0)= I/psf_sigma(0) * s.DXS(i) * s.Y(j);
    hess(0,1)= I * s.DX(i) * s.DY(j);
    hess(1,1)= I/psf_sigma(1) * s.DYS(j) * s.X(i);
    hess(0,2)= s.DX(i) * s.Y(j); 
    hess(1,2)= s.DY(j) * s.X(i); 
}


#endif /* _GAUSS2DMODEL_H */
