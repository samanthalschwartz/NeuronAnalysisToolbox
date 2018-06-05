/** @file DEstimator_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @author Peter Relich (physx.grad\@gmail.com)
 * @date 08-2014
 * @brief The class declaration and inline and templated functions for DEstimator_Iface.
 */

#ifndef _DESTIMATOR_IFACE
#define _DESTIMATOR_IFACE

#include "Mex_Iface.h"
#include "destimator.h"

class DEstimator_Iface : public Mex_Iface
{
public:
    typedef double FloatT;
    DEstimator_Iface();
private:
    DEstimator<FloatT> *obj;

    //Abstract member functions inherited from Mex_Iface
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);
    //Exposed method calls
    void objLLH();
    void objLLHdim();
    void objStaticLLH(const std::string &method);
};

#endif /* _DESTIMATOR_IFACE */
