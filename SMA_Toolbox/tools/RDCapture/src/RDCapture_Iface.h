/** @file RDCapture_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief The class declaration and inline and templated functions for RDCapture_Iface.
 */

#ifndef _RCAPTURE_IFACE
#define _RCAPTURE_IFACE

#include "Mex_Iface.h"
#include "rdcapture.h"

class RDCapture_Iface : public Mex_Iface
{
public:
    typedef double FloatT;
    RDCapture_Iface();
private:
    RDCapture *obj;

    //Abstract member functions inherited from Mex_Iface
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);
    //Exposed method calls
    void objSurvivalProb();
    void objMu();
    void objNu();
    
    void objSimulate();
    /* Bessel Function access for debugging */
    void objI0();
    void objI1();
    void objK0();
    void objK1();
    void objlogI0();
    void objlogI1();
    void objlogK0();
    void objlogK1();
    void objI0K0();
    void objI0K1();
    void objI1K0();
    void objI1K1();
};

#endif /* _RCAPTURE_IFACE */
