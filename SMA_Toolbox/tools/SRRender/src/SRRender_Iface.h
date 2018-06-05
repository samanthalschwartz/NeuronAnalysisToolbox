/** @file SRRender_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-21-2014
 * @brief The class declaration and inline and templated functions for SRRender2D_Iface and SRRenderHS_Iface.
 */

#ifndef _SRRENDER_IFACE
#define _SRRENDER_IFACE

#include "Mex_Iface.h"
#include "srrender.h"

template<class FloatT>
class SRRender2D_Iface : public Mex_Iface
{
public:
    SRRender2D_Iface();
private:
    SRRender2D<FloatT> *obj;
    //Abstract member functions inherited from Mex_Iface
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);

    //Exposed static methods
    void objRenderHist();
    void objRenderGauss();
    void objRenderHistMovie();
    void objRenderGaussMovie();
};


template<class FloatT>
SRRender2D_Iface<FloatT>::SRRender2D_Iface()
: Mex_Iface("SRRender2D_Iface")
{
    staticmethodmap["renderHist"] = boost::bind(&SRRender2D_Iface<FloatT>::objRenderHist, this);
    staticmethodmap["renderGauss"] = boost::bind(&SRRender2D_Iface<FloatT>::objRenderGauss, this);
    staticmethodmap["renderHistMovie"] = boost::bind(&SRRender2D_Iface<FloatT>::objRenderHistMovie, this);
    staticmethodmap["renderGaussMovie"] = boost::bind(&SRRender2D_Iface<FloatT>::objRenderGaussMovie, this);
}

template<class FloatT>
void SRRender2D_Iface<FloatT>::objConstruct()
{
    checkNumArgs(1,0);
    auto *obj = new SRRender2D<FloatT>();
    outputMXArray(Handle<SRRender2D<FloatT>>::makeHandle(obj));
}

template<class FloatT>
void SRRender2D_Iface<FloatT>::objDestroy()
{
    checkNumArgs(0,1);
    Handle<SRRender2D<FloatT>>::destroyObject(rhs[0]);
}

template<class FloatT>
void SRRender2D_Iface<FloatT>::getObjectFromHandle(const mxArray *mxhandle)
{
    obj = Handle<SRRender2D<FloatT>>::getObject(mxhandle);
}

template<class FloatT>
void SRRender2D_Iface<FloatT>::objRenderHist()
{
    // [in] points: mat with n rows and 3 (or more) columns [I, x, y]
    // [in] roi=[xmin, xmax, ymin, ymax] - 4-element vector giving effective region of interest that
    //                                     for the images generated.  This should be exact to the boundaries of im.
    // [in/out] im: a image with arbitrary size but should match aspect ratio of the internal size.  This
    //          image will be modified in-place.

    checkNumArgs(0,3);
    auto points=getMat<FloatT>();
    auto roi=getVec<FloatT>();
    auto im=getMat<FloatT>();
    obj->renderHist(points,roi,im);
}

template<class FloatT>
void SRRender2D_Iface<FloatT>::objRenderGauss()
{
    // [in] points: mat with n rows and 5 (or more) columns [I, x, y, sigma_x, sigma_y]
    // [in] roi=[xmin, xmax, ymin, ymax] - 4-element vector giving effective region of interest that
    //                                     for the images generated.  This should be exact to the boundaries of im.
    // [in] sigmaAccuracy: floating point >0.  Gives accuracy at which gaussians will be rendered
    // [in/out] im: a image with arbitrary size but should match aspect ratio of the internal size.  This
    //          image will be modified in-place.

    checkNumArgs(0,4);
    auto points=getMat<FloatT>();
    auto roi=getVec<FloatT>();
    auto sigmaAccuracy = getDouble();
    auto im=getMat<FloatT>();
    obj->renderGauss(points,roi,im, sigmaAccuracy);
}

template<class FloatT>
void SRRender2D_Iface<FloatT>::objRenderHistMovie()
{
    // [in] points: mat with n rows and 6 (or more) columns [I, x, y, sigma_x, sigma_y, frameIdx]
    //              frame indexs are 0-based.  sigma_x and sigma_y are ignored but must be included,
    // [in] roi=[xmin, xmax, ymin, ymax] - 4-element vector giving effective region of interest that
    //                                     for the images generated.  This should be exact to the boundaries of im.
    // [in/out] im: a image sequence (movie) with arbitrary size but should match aspect ratio of the internal size.  This
    //          image will be modified in-place.  The number of frames should match the frame indexs from points

    checkNumArgs(0,3);
    auto points=getMat<FloatT>();
    auto roi=getVec<FloatT>();
    auto im=getStack<FloatT>();
    obj->renderHistMovie(points,roi,im);
}


template<class FloatT>
void SRRender2D_Iface<FloatT>::objRenderGaussMovie()
{
    // [in] points: mat with n rows and 6 (or more) columns [I, x, y, sigma_x, sigma_y, frameIdx]
    //              frame indexs are 0-based.
    // [in] roi=[xmin, xmax, ymin, ymax] - 4-element vector giving effective region of interest that
    //                                     for the images generated.  This should be exact to the boundaries of im.
    // [in/out] im: a image sequence (movie) with arbitrary size but should match aspect ratio of the internal size.  This
    //          image will be modified in-place.  The number of frames should match the frame indexs from points

    checkNumArgs(0,3);
    auto points=getMat<FloatT>();
    auto roi=getVec<FloatT>();
    auto im=getStack<FloatT>();
    obj->renderGaussMovie(points,roi,im);
}



#endif /* _SRRENDER_IFACE */
