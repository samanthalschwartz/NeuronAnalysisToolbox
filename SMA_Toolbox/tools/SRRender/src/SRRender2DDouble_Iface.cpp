/** @file SRRender2DDouble_Iface.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 03-2015
 *  @brief The entry point for SRRender2DDouble_Iface mex module.
 * 
 */
#include "SRRender_Iface.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    SRRender2D_Iface<double> iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}
