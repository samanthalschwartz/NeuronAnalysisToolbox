/** @file Boxxer2D_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-24-2014
 * @brief The class declaration and inline and templated functions for Boxxer2D_Iface.
 */

#ifndef _BOXXER2D_IFACE
#define _BOXXER2D_IFACE

#include "Mex_Iface.h"
#include "Boxxer2D.h"

class Boxxer2D_Iface : public Mex_Iface
{
public:
    typedef float FloatT;
    typedef Boxxer2D<FloatT>::IMatT IMatT;
    typedef Boxxer2D<FloatT>::VecT VecT;
    Boxxer2D_Iface();
private:
    Boxxer2D<FloatT> *obj;
    /* Inherited from Mex_Iface */
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);

    /* Public method wrappers */
    void objSetDoGSigmaRatio();
    void objFilterScaledLoG();
    void objFilterScaledDoG();
    void objScaleSpaceLoGMaxima();
    void objScaleSpaceDoGMaxima();
    
    /* Static method wrappers */
    void objFilterLoG();
    void objFilterDoG();
    void objFilterGauss();
    void objEnumerateImageMaxima();
    
};

#endif /* _BOXXER2D_IFACE */
