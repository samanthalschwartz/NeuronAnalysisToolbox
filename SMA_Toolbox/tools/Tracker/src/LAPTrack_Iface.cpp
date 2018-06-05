/** @file LAPTrack_Iface.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 05-2015
 *  @brief The entry point for LAPTrack_Iface mex module.
 * 
 */
#include "Tracker_Iface.h"
#include "LAPTrack.h"

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    Tracker_Iface<LAPTrack> iface("LAPTrack");
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}

