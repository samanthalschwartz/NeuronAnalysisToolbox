/** @file rdcapture_cache.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief  Cache RDCaptrue computations for certain times.
 *
 */
#include "rdcapture_cache.h"


RDCaptureCache::RDCaptureCache(const RDCapture &_rd) 
    : rd(std::make_shared<RDCapture>(_rd)) 
{
    
}
