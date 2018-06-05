/** @file InteractionChangePoint_Iface.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 07-2015
 *  @brief The entry point for InteractionChangePointe mex module.
 *
 */
#include "InteractionChangePoint_Iface.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    InteractionChangePoint_Iface iface;
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}
