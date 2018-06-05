/** @file Boxxer2D_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-21-2014
 * @brief The class declaration and inline and templated functions for Boxxer2D_Iface.
 */
#include "Boxxer2D_Iface.h"
#include "Handle.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    Boxxer2D_Iface iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}

Boxxer2D_Iface::Boxxer2D_Iface() 
    : Mex_Iface("Boxxer2D_Iface") 
{
    methodmap["setDoGSigmaRatio"]=boost::bind(&Boxxer2D_Iface::objSetDoGSigmaRatio, this);
    methodmap["filterScaledLoG"]=boost::bind(&Boxxer2D_Iface::objFilterScaledLoG, this);
    methodmap["filterScaledDoG"]=boost::bind(&Boxxer2D_Iface::objFilterScaledDoG, this);
    methodmap["scaleSpaceLoGMaxima"]=boost::bind(&Boxxer2D_Iface::objScaleSpaceLoGMaxima, this);
    methodmap["scaleSpaceDoGMaxima"]=boost::bind(&Boxxer2D_Iface::objScaleSpaceDoGMaxima, this);

    staticmethodmap["filterLoG"]=boost::bind(&Boxxer2D_Iface::objFilterLoG, this);
    staticmethodmap["filterDoG"]=boost::bind(&Boxxer2D_Iface::objFilterDoG, this);
    staticmethodmap["filterGauss"]=boost::bind(&Boxxer2D_Iface::objFilterGauss, this);
    staticmethodmap["enumerateImageMaxima"]=boost::bind(&Boxxer2D_Iface::objEnumerateImageMaxima, this);
}


void Boxxer2D_Iface::objConstruct()
{
    checkNumArgs(1,2);
    auto imsize=getIVec();
    auto sigma=getFMat();
    auto boxxer=new Boxxer2D<FloatT>(imsize, sigma);
    outputMXArray( Handle<Boxxer2D<FloatT>>::makeHandle(boxxer));
}

void Boxxer2D_Iface::objDestroy()
{
    checkNumArgs(0,1);
    Handle<Boxxer2D<FloatT>>::destroyObject(rhs[0]);
}

void Boxxer2D_Iface::getObjectFromHandle(const mxArray *mxhandle)
{
    obj=Handle<Boxxer2D<FloatT>>::getObject(mxhandle);
}

void Boxxer2D_Iface::objSetDoGSigmaRatio()
{
    // [in] sigma_ratio: a new sigma_ratio>1
    checkNumArgs(0,1);
    auto sigma_ratio=getDouble();
    obj->setDoGSigmaRatio(sigma_ratio);
}

void Boxxer2D_Iface::objFilterScaledLoG()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [out] fimage: a single stack of imsize x nScales filtered frames. Size 4D: [x y S t]
    checkNumArgs(1,1);
    auto ims=getFStack();
    auto fims=makeFHyperStack(ims.n_rows, ims.n_cols, obj->nScales, ims.n_slices);
    obj->filterScaledLoG(ims,fims);
}

void Boxxer2D_Iface::objFilterScaledDoG()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [out] fimage: a single stack of imsize x nScales filtered frames. Size 4D: [x y S t]
    checkNumArgs(1,1);
    auto ims=getFStack();
    auto fims=makeFHyperStack(ims.n_rows, ims.n_cols, obj->nScales, ims.n_slices);
    obj->filterScaledDoG(ims,fims);
}

void Boxxer2D_Iface::objScaleSpaceLoGMaxima()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [in] neighborhoodSize: an odd integer. Will be converted to uint32.  Acceptable values are in ValidMaximaNeighborhoodSizes.  (default=3)
    // [in] scaleNeighborhoodSize: an odd integer. Will be converted to uint32.  Acceptable values are in ValidMaximaNeighborhoodSizes.  (default=3)
    // [out] maxima: a (dim+1)xN matrix of maxima where rows are X, Y, ..., T and columns are different maxima detected.
    // [out] max_vals; a Nx1 Vector giving the value at each maxima.
    checkNumArgs(2,3);
    auto ims=getFStack();
    auto neighborhoodSize = getInt();
    auto scaleNeighborhoodSize = getInt();
    IMatT maxima;
    VecT max_vals;
    obj->scaleSpaceLoGMaxima(ims, maxima, max_vals, neighborhoodSize, scaleNeighborhoodSize);
    outputIMat(maxima);
    outputFVec(max_vals);
}

void Boxxer2D_Iface::objScaleSpaceDoGMaxima()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [in] neighborhoodSize: an odd integer. Will be converted to uint32.  Acceptable values are in ValidMaximaNeighborhoodSizes.  (default=3)
    // [in] scaleNeighborhoodSize: an odd integer. Will be converted to uint32.  Acceptable values are in ValidMaximaNeighborhoodSizes.  (default=3)
    // [out] maxima: a (dim+1)xN matrix of maxima where rows are X, Y, ..., T and columns are different maxima detected.
    // [out] max_vals; a Nx1 Vector giving the value at each maxima.
    checkNumArgs(2,3);
    auto ims=getFStack();
    auto neighborhoodSize = getInt();
    auto scaleNeighborhoodSize = getInt();
    IMatT maxima;
    VecT max_vals;
    obj->scaleSpaceDoGMaxima(ims, maxima, max_vals, neighborhoodSize, scaleNeighborhoodSize);
    outputIMat(maxima);
    outputFVec(max_vals);
}


void Boxxer2D_Iface::objFilterLoG()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [in] sigma: a col vector [nx1] of sigma size to filter for [sigma_rows, sigma_cols]
    // [out] fimage: a single stack of imsize shaped filtered frames
    checkNumArgs(1,2);
    auto ims=getFStack();
    auto sigma=getFVec();
    auto fims=makeFStack(ims.n_rows, ims.n_cols, ims.n_slices);

    Boxxer2D<FloatT>::filterLoG(ims,fims,sigma);
}

void Boxxer2D_Iface::objFilterDoG()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [in] sigma: a col vector [nx1] of sigma size to filter for [sigma_rows, sigma_cols]
    // [in] sigmaRatio: a scalar giving the ratio of sigmas in the DoG method:
    // [out] fimage: a single stack of imsize shaped filtered frames
    checkNumArgs(1,3);
    auto ims=getFStack();
    auto sigma=getFVec();
    auto sigma_ratio=getDouble();
    auto fims=makeFStack(ims.n_rows, ims.n_cols, ims.n_slices);
    Boxxer2D<FloatT>::filterDoG(ims,fims,sigma,sigma_ratio);
}

void Boxxer2D_Iface::objFilterGauss()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [in] sigma: a col vector [nx1] of sigma size to filter for [sigma_rows, sigma_cols]
    // [out] fimage: a single stack of imsize shaped filtered frames
    checkNumArgs(1,2);
    auto ims=getFStack();
    auto sigma=getFVec();
    auto fims=makeFStack(ims.n_rows, ims.n_cols, ims.n_slices);
    Boxxer2D<FloatT>::filterGauss(ims,fims,sigma);
}


void Boxxer2D_Iface::objEnumerateImageMaxima()
{
    // [in] image: a single stack of imsize shaped frames, last dimension is time
    // [in] neighborhoodSize: an odd integer. Will be converted to uint32.  Acceptable values are in ValidMaximaNeighborhoodSizes.  (default=3)
    // [out] maxima: a (dim+1)xN matrix of maxima where rows are X, Y, ..., T and columns are different maxima detected.
    // [out] max_vals; a Nx1 Vector giving the value at each maxima.
    checkNumArgs(2,2);
    auto ims=getFStack();
    auto neighborhoodSize=getInt();
    IMatT maxima;
    VecT max_vals;
    obj->enumerateImageMaxima(ims,maxima, max_vals, neighborhoodSize);
    outputIMat(maxima);
    outputFVec(max_vals);
}


