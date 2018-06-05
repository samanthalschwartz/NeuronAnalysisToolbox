/** @file Boxxer3D_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-01-2014
 * @brief The class declaration and inline and templated functions for Boxxer3D_Iface.
 */

#ifndef _BOXXER3D_IFACE
#define _BOXXER3D_IFACE

#include "Mex_Iface.h"
#include "Boxxer3D.h"

class Boxxer3D_Iface : public Mex_Iface
{
public:
    typedef float FloatT;
    typedef Boxxer3D<FloatT>::IMatT IMatT;
    typedef Boxxer3D<FloatT>::VecT VecT;
    Boxxer3D_Iface();
private:
    Boxxer3D<FloatT> *obj;
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



#endif /* _BOXXER3D_IFACE */
