/** @file Gauss2DModel.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 03-13-2014
 * @brief The class definition and template Specializations for Gauss2DModel
 */

#include "Gauss2DModel.h"
#include "stencil.h"

const std::vector<std::string> Gauss2DModel::param_names({ "x", "y", "I", "bg" });

Gauss2DModel::Gauss2DModel(const IVecT &size, const VecT &psf_sigma)
    : PointEmitter2DModel(4, size, psf_sigma),
      gaussian_Xstencil(make_gaussian_stencil(size(0),psf_sigma(0))),
      gaussian_Ystencil(make_gaussian_stencil(size(1),psf_sigma(1)))
{
    /* Initialize MCMC step sizes */
    num_candidate_sampling_phases=2;
}


Gauss2DModel::Stencil::Stencil(const Gauss2DModel &model_,
                               const Gauss2DModel::ParamT &theta,
                               bool _compute_derivatives)
      : model(&model_),theta(theta)
{
    int szX=model->size(0);
    int szY=model->size(1);
    dx=make_d_stencil(szX, x());
    dy=make_d_stencil(szY, y());
    X=make_X_stencil(szX, dx,model->psf_sigma(0));
    Y=make_X_stencil(szY, dy,model->psf_sigma(1));
    if(_compute_derivatives) compute_derivatives();
}

void Gauss2DModel::Stencil::compute_derivatives()
{
    if(derivatives_computed) return;
    derivatives_computed=true;
    int szX=model->size(0);
    int szY=model->size(1);
    double sigmaX=model->psf_sigma(0);
    double sigmaY=model->psf_sigma(1);
    Gx=make_G_stencil(szX, dx,sigmaX);
    Gy=make_G_stencil(szY, dy,sigmaY);
    DX=make_DX_stencil(szX, Gx,sigmaX);
    DY=make_DX_stencil(szY, Gy,sigmaY);
    DXS=make_DXS_stencil(szX, dx, Gx,sigmaX);
    DYS=make_DXS_stencil(szY, dy, Gy,sigmaY);
}

std::ostream& operator<<(std::ostream &out, const Gauss2DModel::Stencil &s)
{
    int w=8;
    print_vec_row(out,s.theta,"Theta:",w,TERM_WHITE);
    print_vec_row(out,s.dx,"dx:",w,TERM_CYAN);
    print_vec_row(out,s.dy,"dy:",w,TERM_CYAN);
    print_vec_row(out,s.X,"X:",w,TERM_CYAN);
    print_vec_row(out,s.Y,"Y:",w,TERM_CYAN);
    if(s.derivatives_computed) {
        print_vec_row(out,s.Gx,"Gx:",w,TERM_BLUE);
        print_vec_row(out,s.Gy,"Gy:",w,TERM_BLUE);
        print_vec_row(out,s.DX,"DX:",w,TERM_BLUE);
        print_vec_row(out,s.DY,"DY:",w,TERM_BLUE);
        print_vec_row(out,s.DXS,"DXS:",w,TERM_BLUE);
        print_vec_row(out,s.DYS,"DYS:",w,TERM_BLUE);
    }
    return out;
}




void
Gauss2DModel::pixel_hess_update(int i, int j, const Stencil &s, double dm_ratio_m1, double dmm_ratio, ParamT &grad, ParamMatT &hess) const
{
    /* Caclulate pixel derivative */
    auto pgrad=make_param();
    pixel_grad(i,j,s,pgrad);
    double I=s.I();
    /* Update grad */
    grad+=dm_ratio_m1*pgrad;
    /* Update hess */
    hess(0,0)+=dm_ratio_m1 * I/psf_sigma(0) * s.DXS(i) * s.Y(j);
    hess(0,1)+=dm_ratio_m1 * I * s.DX(i) * s.DY(j);
    hess(1,1)+=dm_ratio_m1 * I/psf_sigma(1) * s.DYS(j) * s.X(i);
    hess(0,2)+=dm_ratio_m1 * pgrad(0) / I; 
    hess(1,2)+=dm_ratio_m1 * pgrad(1) / I; 
    //This is the pixel-gradient dependenent part of the hessian
    for(int c=0; c<(int)hess.n_cols; c++) for(int r=0; r<=c; r++)
        hess(r,c) -= dmm_ratio * pgrad(r) * pgrad(c);
}



Gauss2DModel::Stencil
Gauss2DModel::initial_theta_estimate(const ImageT &im, const ParamT &theta_init) const
{
    double x_pos=0, y_pos=0, I=0, bg=0;
    double min_bg=1; //default minimum background.  Will be updated only if estimate_gaussian_2Dmax is called.
//     std::cout<<"Theta_init: "<<theta_init.t()<<" -->";
    if (!theta_init.is_empty()) {
        x_pos = theta_init(0);
        y_pos = theta_init(1);
        I = theta_init(2);
        bg = theta_init(3);
    }
    if(x_pos<=0 || x_pos>size(0) || y_pos<=0 || y_pos>size(1)){ //Invlaid positions, estimate them
//         std::cout<<"Full init\n";
        int px_pos[2];
        estimate_gaussian_2Dmax(im, gaussian_Xstencil, gaussian_Ystencil, px_pos, min_bg);
        refine_gaussian_2Dmax(im, gaussian_Xstencil, gaussian_Ystencil, px_pos);
        x_pos = static_cast<double>(px_pos[0])+0.5;
        y_pos = static_cast<double>(px_pos[1])+0.5;
        auto unit_im = unit_model_image(size,px_pos,psf_sigma);
        bg = estimate_background(im, unit_im, min_bg);
        I = estimate_intensity(im, unit_im, bg);
    } else if(I<=0 || bg<=0) {
//         std::cout<<"Intenisty init\n";
        int px_pos[2];
        px_pos[0] = static_cast<int>(floor(x_pos));
        px_pos[1] = static_cast<int>(floor(y_pos));
        auto unit_im = unit_model_image(size,px_pos,psf_sigma);
        bg = estimate_background(im, unit_im, min_bg);
        I = estimate_intensity(im, unit_im, bg);
    } /*else {
        std::cout<<"Null init\n";
    }*/
    auto theta= make_stencil(x_pos, y_pos, I, bg);
//     std::cout<<"ThetaFinal: "<<theta.theta.t()<<"\n";
    return theta;
}


void Gauss2DModel::sample_candidate_theta(int sample_index, RNG &rng, ParamT &candidate_theta, double scale) const
{
    int phase=sample_index%num_candidate_sampling_phases;
    switch(phase) {
        case 0:  //change x,y
            candidate_theta(0)+=generate_normal(rng,0.0,candidate_eta_x*scale);
            candidate_theta(1)+=generate_normal(rng,0.0,candidate_eta_y*scale);
            break;
        case 1: //change I, bg
            candidate_theta(2)+=generate_normal(rng,0.0,candidate_eta_I*scale);
            candidate_theta(3)+=generate_normal(rng,0.0,candidate_eta_bg*scale);
    }
}

