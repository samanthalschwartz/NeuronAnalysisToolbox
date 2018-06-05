/** @file KDTreeDouble_Iface.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 09-2014
 *  @brief The entry point for KDTreeDouble_Iface mex module.
 * 
 */
#include "KDTree_Iface.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    KDTree_Iface<double> iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}
