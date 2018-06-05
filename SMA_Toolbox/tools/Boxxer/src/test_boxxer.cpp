
#include "filter_kernels.h"
#include "GaussFilter.h"
#include "Maxima.h"
#include "Boxxer2D.h"
#include "Boxxer3D.h"

using std::cout;
using std::endl;
using namespace arma;
// #include "maxima.h"

// void time1D(int N, double sigma, double tol)
// {
//     VecT data(N), fdata(N);
//     data.zeros();
//     data(0)=1;
//     data(3*N/4)=1;
//     int hw=compute_FIR_width(sigma, tol);
//     VecT kernel=compute_FIR_kernel(sigma, tol);
//     for(int k=1; k<N; k++){
//         gaussFIR_1D(data, fdata, kernel);
//         gaussFIR_1D_arma(data, fdata, kernel);
//         gaussFIR_1D_small(data.n_elem, data.memptr(), fdata.memptr(), kernel.n_elem-1, kernel.memptr());
//         gaussFIR_1D_inplace_arma(fdata, kernel);
//         gaussFIR_1D_inplace(data.n_elem, data.memptr(),kernel.n_elem-1, kernel.memptr());
//     }
// }


// void time2D(int N, double sigma, double tol)
// {
//     MatT data(N,N), fdata(N,N);
//     data.randu();
//     VecT kernel=compute_FIR_kernel(sigma, tol);
//     gaussFIR_2Dx(data, fdata, kernel);
//     gaussFIR_2Dx_small(data, fdata, kernel);
//     gaussFIR_2Dx_arma(data, fdata, kernel);
//     gaussFIR_2Dy_rowmajor(data, fdata, kernel);
//     gaussFIR_2Dy_colmajor(data, fdata, kernel);
//     gaussFIR_2Dy_small(data, fdata, kernel);
//     gaussFIR_2Dy(data, fdata, kernel);
// }
/*
void time3D(int N, double sigma, double tol)
{
    CubeT data(N,N,N), fdata(N,N,N);
    data.randu();
    VecT kernel=compute_FIR_kernel(sigma, tol);
    gaussFIR_3Dx_small(data, fdata, kernel);
    gaussFIR_3Dy_small(data, fdata, kernel);
    gaussFIR_3Dz_small(data, fdata, kernel);
}

void test3D(int N, double sigma, double tol)
{
    CubeT data(N,N,N), ref(N,N,N), fdata(N,N,N);
    data.randu();
    VecT kernel=compute_FIR_kernel(sigma, tol);
    gaussFIR_3Dx_small(data, ref, kernel);
    gaussFIR_3Dx(data, fdata, kernel);
    for(int z=0; z<N; z++) for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y,z)-fdata(x,y,z))>1e-9) printf("(%i,%i,%i): ref=%.9g != %.9g\n",x,y,z,ref(x,y,z),fdata(x,y,z));
    gaussFIR_3Dy_small(data, ref, kernel);
    gaussFIR_3Dy(data, fdata, kernel);
    for(int z=0; z<N; z++) for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y,z)-fdata(x,y,z))>1e-9) printf("(%i,%i,%i): ref=%.9g != %.9g\n",x,y,z,ref(x,y,z),fdata(x,y,z));
    gaussFIR_3Dz_small(data, ref, kernel);
    gaussFIR_3Dz(data, fdata, kernel);
    for(int z=0; z<N; z++) for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y,z)-fdata(x,y,z))>1e-9) printf("(%i,%i,%i): ref=%.9g != %.9g\n",x,y,z,ref(x,y,z),fdata(x,y,z));
}


void test2D(int N, double sigma, double tol)
{
    MatT data(N,N), ref(N,N), fdata(N,N);
    data.randu();
    VecT kernel=compute_FIR_kernel(sigma, tol);

    Test 2Dx
    gaussFIR_2Dx_small(data, ref, kernel);
    gaussFIR_2Dx(data, fdata, kernel);
    for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y)-fdata(x,y))>1e-9) printf("(%i,%i): ref=%.9g != %.9g\n",x,y,ref(x,y),fdata(x,y));
    gaussFIR_2Dx_arma(data, fdata, kernel);
    for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y)-fdata(x,y))>1e-9) printf("(%i,%i): ref=%.9g != %.9g\n",x,y,ref(x,y),fdata(x,y));

    Test 2Dy
    gaussFIR_2Dy_small(data, ref, kernel);
    gaussFIR_2Dy(data, fdata, kernel);
    for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y)-fdata(x,y))>1e-9) printf("(%i,%i): ref=%.9g != %.9g\n",x,y,ref(x,y),fdata(x,y));
    gaussFIR_2Dy_rowmajor(data, fdata, kernel);
    for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y)-fdata(x,y))>1e-9) printf("(%i,%i): ref=%.9g != %.9g\n",x,y,ref(x,y),fdata(x,y));
    gaussFIR_2Dy_colmajor(data, fdata, kernel);
    for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y)-fdata(x,y))>1e-9) printf("(%i,%i): ref=%.9g != %.9g\n",x,y,ref(x,y),fdata(x,y));
    gaussFIR_2Dy(data.n_rows, data.n_cols, data.memptr(), fdata.memptr(), kernel.n_elem-1,kernel.memptr());
    for(int y=0; y<N; y++) for(int x=0; x<N; x++) if(fabs(ref(x,y)-fdata(x,y))>1e-9) printf("(%i,%i): ref=%.9g != %.9g\n",x,y,ref(x,y),fdata(x,y));
}


void test1D(int N, double sigma, double tol)
{
    VecT data(N), fdata(N);
    data.randu();
        data(0)=1;
        data(1)=.5;
        data(2)=.25;
        data(3)=.125;
        data(4)=1.;
        data(10)=1;
        data(12)=1;
        data(N-1)=1;
    int hw=compute_FIR_width(sigma, tol);
    VecT kernel=compute_FIR_kernel(sigma, tol);

    gaussFIR_1D(data, fdata, kernel);

    std::cout<<"Sigma:"<<sigma<<" Tol:"<<tol<<" HW:"<<hw<<"\n";
    print_vec_row(std::cout, kernel, "Kernel: ", 20, TERM_DIM_MAGENTA);
    print_vec_row(std::cout, data,  "Data: ", 20, TERM_BLUE);
    print_vec_row(std::cout, fdata, "FData: ", 20, TERM_MAGENTA);
    gaussFIR_1D_inplace(data.n_elem, data.memptr(), kernel.n_elem-1, kernel.memptr());
    std::cout<<"Sum(Data):"<<arma::accu(data)<<" Sum(fData):"<<arma::accu(fdata)<<"\n";
    print_vec_row(std::cout, data, "Data: ", 20, TERM_RED);
    std::cout<<"Sum(Data):"<<arma::accu(data)<<"\n";
    std::cout<<((data-fdata)<1e-12)<<"\n";
    assert(arma::all(arma::abs(data-fdata)<1e-12));
}

void testNMS()
{
    int N=100;
    MatT data(N, N, arma::fill::randu);
    data(1,1)=10; data(1,2)=11;
    data(2,2)=22; data(4,4)=3; data(4,3)=2;
    //     std::cout<<"Data: \n"<<data<<"\n";
    arma::umat maxima=detectLocalMaxima_3x3_2D(data);
    //     std::cout<<"Maxima: \n"<<maxima.t()<<"\n";
    checkLocalMaxima_3x3_2D(data, maxima);
    maxima=detectLocalMaxima_3x3_2D_simple(data);
    //     std::cout<<"\n\nMaxima: \n"<<maxima.t()<<"\n";
    checkLocalMaxima_3x3_2D(data, maxima);
    maxima=detectLocalMaxima_3x3_2D_simple2(data);
    //     std::cout<<"\n\nMaxima: \n"<<maxima.t()<<"\n";
    checkLocalMaxima_3x3_2D(data, maxima);
}*/

#include <limits>
// template <class FloatT>
// void test1D(int N, FloatT sigma, FloatT tol)
// {
//     arma::Col<FloatT> data(N), fdata(N), fdata_slow(N);
//     data.randu();
//     GaussFIRFilter<FloatT> filt;
//     int hw=filt.compute_FIR_hw(sigma, tol);
//     arma::Col<FloatT> kernel=filt.compute_FIR_kernel(sigma, hw);
// 
//     gaussFIR_1D<FloatT>(data, fdata, kernel);
//     gaussFIR_1D_small<FloatT>(data.n_elem, data.memptr(), fdata_slow.memptr(), kernel.n_elem-1, kernel.memptr());
//     printf("Testing 1D filter:\n");
//     cout<<"Kernel:"<<kernel<<endl;
//     cout<<"IN:"<<data<<endl;
//     for(int k=1; k<N; k++){
//         if(fabs(fdata(k)-fdata_slow(k))> std::numeric_limits<FloatT>::epsilon())
//             printf("Fast (%i):%.17f  != Slow (%i):%.17f\n",k,fdata(k),k,fdata_slow(k));
//     }
//     cout<<"OUT:"<<fdata<<endl;
// }


void testGaussFilter2D()
{
    typedef float TestFloat;
    GaussFilter2D<TestFloat>::IVecT size={100,100};
    GaussFilter2D<TestFloat>::FVecT sigma={0.8,1.14};
    GaussFilter2D<TestFloat> gauss_filt(size, sigma);
    Maxima2D<TestFloat> maxima2D(size);
    cout<<gauss_filt<<"\n";
    auto image=gauss_filt.make_image();
    image.randu();
//     cout<<"IN:\n"<<image<<"\n";
    
    int Nmaxima=maxima2D.find_maxima(image);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    maxima2D.test_maxima(image);

    auto out=gauss_filt.make_image();
    gauss_filt.test_filter(image);
    gauss_filt.filter(image, out);
//     std::cout<<"OUT:\n"<<out<<"\n";
    int Nmaxima_out=maxima2D.find_maxima(out);
    cout<<"NMaxima_out: "<<Nmaxima_out<<endl;
    maxima2D.test_maxima(out);
}

void testLoGFilter2D()
{
    typedef float TestFloat;
    LoGFilter2D<TestFloat>::IVecT size={100,100};
    LoGFilter2D<TestFloat>::FVecT sigma={0.8,1.14};
    LoGFilter2D<TestFloat> log_filt(size, sigma);
    Maxima2D<TestFloat> maxima2D(size);
    std::cout<<log_filt<<"\n";
    auto image=log_filt.make_image();
    image.randu();
    int Nmaxima=maxima2D.find_maxima(image);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    maxima2D.test_maxima(image);

    auto out=log_filt.make_image();
    log_filt.test_filter(image);
    log_filt.filter(image, out);

    Nmaxima=maxima2D.find_maxima(out);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    maxima2D.test_maxima(out);
}

void testGaussFilter3D()
{
    typedef float TestFloat;
    GaussFilter3D<TestFloat>::IVecT size={100,100,100};
    GaussFilter3D<TestFloat>::FVecT sigma={0.8,1.14,1.0};
    GaussFilter3D<TestFloat> gauss_filt(size, sigma);
    Maxima3D<TestFloat> maxima3D(size);
    std::cout<<gauss_filt<<"\n";
    auto image=gauss_filt.make_image();
    image.randu();

    int Nmaxima=maxima3D.find_maxima(image);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    maxima3D.test_maxima(image);


    auto out=gauss_filt.make_image();
    gauss_filt.test_filter(image);
    gauss_filt.filter(image, out);

    Nmaxima=maxima3D.find_maxima(out);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    maxima3D.test_maxima(out);

}


void testLoGFilter3D()
{
    typedef float TestFloat;
    LoGFilter3D<TestFloat>::IVecT size={100,100,100};
    LoGFilter3D<TestFloat>::FVecT sigma={0.8,1.14,1.0};
    LoGFilter3D<TestFloat> log_filt(size, sigma);
    std::cout<<log_filt<<"\n";
    auto image=log_filt.make_image();
    image.randu();
    auto out=log_filt.make_image();
    log_filt.test_filter(image);
    log_filt.filter(image, out);
}


void testBoxxer2D()
{
    int nT=7;
    int sz=100;
    typedef float TestFloat;
    Boxxer2D<TestFloat>::IVecT size={sz,sz};
    Boxxer2D<TestFloat>::VecT sigma={1.0, 1.0};
    Boxxer2D<TestFloat>::ImageStackT ims(sz,sz,nT);
    ims.randu();
    ims(sz/2,sz/2,0)=1.;
    

    Boxxer2D<TestFloat> boxxer(size, sigma);
//     cout<<"IN: \n"<<ims<<endl;
    auto out=boxxer.make_image_stack(nT);
    Boxxer2D<TestFloat>::filterDoG(ims,out,sigma, 1.6);
//     cout<<"Out: \n"<<out<<endl;
    Boxxer2D<TestFloat>::IMatT maxima;
    Boxxer2D<TestFloat>::VecT max_vals;
    Boxxer2D<TestFloat>::enumerateImageMaxima(out,maxima, max_vals, 5);
    
    int Nmaxima=static_cast<int>(maxima.n_cols);
    cout<<"Boxxer2D: Size:["<<size(0)<<","<<size(1)<<","<<nT<<"]\n";
    cout<<"Nmaxima: "<<Nmaxima<<endl;
    for(int n=0; n<Nmaxima; n++){
        printf("Maxima[%i]: (%i,%i, %i):%.9g\n",n,maxima(0,n),maxima(1,n),maxima(2,n),max_vals(n));
    }
}

void testScaleSpace2D()
{
    int nT=10;
    int sz=16;
    typedef float TestFloat;
    Boxxer2D<TestFloat>::IVecT size={sz,sz};
    Boxxer2D<TestFloat>::MatT sigma;
    sigma << 1.0 << 1.6 << 2.0<<endr
          << 1.0 << 1.6 << 2.0<<endr;
    Boxxer2D<TestFloat> boxxer(size, sigma);
    auto ims=boxxer.make_image_stack(nT);
    ims.randu();
    ims(sz/2,sz/2,0)=100.;

    Boxxer2D<TestFloat>::IMatT maxima;
    Boxxer2D<TestFloat>::VecT max_vals;
    boxxer.scaleSpaceLoGMaxima(ims, maxima, max_vals, 13, 5);
    
    int Nmaxima=static_cast<int>(maxima.n_cols);
    cout<<"Boxxer2D: Size:["<<size(0)<<","<<size(1)<<","<<sigma.n_cols<<","<<nT<<"]\n";
    cout<<"Nmaxima: "<<Nmaxima<<endl;
    for(int n=0; n<Nmaxima; n++){
        printf("Maxima[%i]: (%i,%i, %i, %i):%.9g\n",n,maxima(0,n),maxima(1,n),maxima(2,n),maxima(3,n),max_vals(n));
    }
}


void testBoxxer3D()
{
    int nT=10;
    int sz=70;
    typedef float TestFloat;
    Boxxer3D<TestFloat>::IVecT size={sz,2*sz, 3*sz};
    Boxxer3D<TestFloat>::VecT sigma={1.0, 1.0, 1.0};
    Boxxer3D<TestFloat> boxxer(size, sigma);
    auto ims=boxxer.make_image_stack(nT);
    for(int n=0; n<nT; n++){
        ims.slice(n).randu();
    }
    ims(sz/2,sz/2,sz/2,0)=1.;


    cout<<"IN: \n"<<ims.slice(0)<<endl;
    //     Boxxer2D<TestFloat>::ImageStackT DoG_out(256,256,10000);
    auto LoG_out=boxxer.make_image_stack(nT);
    Boxxer3D<TestFloat>::filterLoG(ims,LoG_out,sigma);
    //     boxxer.filterDoG(ims,DoG_out,1.6);
    cout<<"Out: \n"<<LoG_out.slice(0)<<endl;
    Boxxer3D<TestFloat>::IMatT maxima;
    Boxxer3D<TestFloat>::VecT max_vals;
    Boxxer3D<TestFloat>::enumerateImageMaxima(LoG_out,maxima, max_vals, 3);

    int Nmaxima=static_cast<int>(maxima.n_cols);
    cout<<"Boxxer3D: Size:["<<size(0)<<","<<size(1)<<","<<size(2)<<","<<nT<<"]\n";
    cout<<"Nmaxima: "<<Nmaxima<<endl;
    for(int n=0; n<Nmaxima; n++){
        printf("Maxima[%i]: (%i,%i,%i,%i):%.9g\n",n,maxima(0,n),maxima(1,n),maxima(2,n),maxima(3,n),max_vals(n));
    }
    boxxer.checkMaxima(LoG_out, maxima, max_vals);
}


void testMaxima2D()
{
    typedef float TestFloat;
    Maxima2D<TestFloat>::IVecT size={10,10};
    Maxima2D<TestFloat> maxima2D(size);
    Maxima2D<TestFloat>::ImageT image(size(0), size(1));
    image.randu();
    image*=2;
    image(0,0)=1;
    image(9,0)=1;
    image(0,9)=1;
    image(9,9)=1;
    image(0,4)=1;
    image(4,0)=1;
    image(9,4)=1;
    image(4,9)=1;
    image(4,4)=1;
    Maxima2D<TestFloat>::IMatT maxima;
    Maxima2D<TestFloat>::VecT max_vals;

    int Nmaxima=maxima2D.find_maxima(image,maxima,max_vals);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    for(int n=0; n<Nmaxima; n++)
        cout<<"("<<maxima(0,n)<<","<<maxima(1,n)<<"):"<<max_vals(n)<<endl;
}


void testMaxima3D()
{
    typedef float TestFloat;
    Maxima3D<TestFloat>::IVecT size={8,8,8};
    Maxima3D<TestFloat> maxima3D(size);
    Maxima3D<TestFloat>::ImageT image(size(0), size(1), size(2));
    image.randu();
//     image*=2;
    //Corners
    image(0,0,0)=1;
    image(0,0,7)=1;
    image(0,7,0)=1;
    image(0,7,7)=1;
    image(7,0,0)=1;
    image(7,0,7)=1;
    image(7,7,0)=1;
    image(7,7,7)=1;

    //Forward Edges
    image(0,0,4)=1;
    image(0,4,0)=1;
    image(0,4,7)=1;
    image(0,7,4)=1;

    //Receding Edges
    image(4,0,0)=1;
    image(4,0,7)=1;
    image(4,7,0)=1;
    image(4,7,7)=1;

    //Backward Edges
    image(7,0,4)=1;
    image(7,4,0)=1;
    image(7,4,7)=1;
    image(7,7,4)=1;

    //Faces
    image(0,4,4)=1;
    image(4,0,4)=1;
    image(4,4,0)=1;
    image(7,4,4)=1;
    image(4,7,4)=1;
    image(4,4,7)=1;

    //Center
    image(4,4,4)=1;

    Maxima2D<TestFloat>::IMatT maxima;
    Maxima2D<TestFloat>::VecT max_vals;

    int Nmaxima=maxima3D.find_maxima(image,maxima,max_vals);
    cout<<"NMaxima: "<<Nmaxima<<endl;
    for(int n=0; n<Nmaxima; n++)
        cout<<"("<<maxima(0,n)<<","<<maxima(1,n)<<","<<maxima(2,n)<<"):"<<max_vals(n)<<endl;
}


int main(){
    arma::arma_rng::set_seed_random();
//     test1D<float>(10, 1., 1e-3);
//     test1D<double>(10, 1., 1e-3);
    testGaussFilter2D();
    testGaussFilter3D();
    testLoGFilter2D();
    testLoGFilter3D();
    testMaxima3D();
    testMaxima2D();
    testBoxxer2D();
    testBoxxer3D();
    testScaleSpace2D();
    return 0;
}
