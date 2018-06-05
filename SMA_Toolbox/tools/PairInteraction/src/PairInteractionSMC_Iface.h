/** @file PairInteractionSMC_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief The class declaration and inline and templated functions for PairInteractionSMC.
 */

#ifndef _PAIRINTERACTIONSMC_IFACE
#define _PAIRINTERACTIONSMC_IFACE

#include "Mex_Iface.h"
#include "pairinteractionsmc.h"

namespace pair_int
{

class PairInteractionSMC_Iface : public Mex_Iface
{
public:
    typedef double FloatT;
    PairInteractionSMC_Iface();
private:
    PairInteractionSMC *obj;

    //Abstract member functions inherited from Mex_Iface
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);
    //Exposed method calls
    void addData();
    void clearData();
    void getData();
    
    void setParams();
    void getParams();
    
    void setPriorState();
    void setPriorPositionGaussian();
    void setPriorPositionRectangular();
    void setPriorPositionCircular();
    
    void runParticleFilter();
    
    void obsLLH();
    void sampleParticle();
    void sampleParticleLLH();
    void getAllParticles();
    
    void computeLLH();
    void computeLLH_debug();

    void simulate();
    void simulateBlur();
    void simulateEC();
    
    void llhG0();
    void llhG();
    void llhH();
    
    void setProposalDist();
    void sampleQ0();
    void sampleQ();

    
    void sampleQ0Z_Y();
    void sampleQZ_X();
    
};

} /* namespace */

#endif /* _PAIRINTERACTIONSMC_IFACE */
