/** @file rdcapture_cache.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief Cache RDCaptrue computations for certain times.
 *
 */

#ifndef _RDCAPTURE_CACHE_H
#define _RDCAPTURE_CACHE_H


#include "rdcapture.h"

class RDCaptureCache{
    using IdxT = arma::uword;
    std::shared_ptr<RDCapture> rd;
    double cache_t;
    double max_r;
    double delta_r; //spacing between bins
    IdxT Nbins;
    VecT bins;
    
    IdxT getBin(double r);
public:
    RDCaptureCache(const RDCapture &_rd);
    
    double survivalProb(double r0, double t)  const;
    
    double survivalLogProb(double r0, double t)  const;
    double captureLogProb(double r0, double t)  const;
}


inline
double RDCaptureCache::survivalLogProb(double r0, double t) const
{
    return log(survivalProb(r0,t));
}

inline
double RDCaptureCache::captureLogProb(double r0, double t) const
{
    double p = survivalProb(r0,t);
    double log1mp;
    if(p<1e-3) {
        log1mp=log(1-p);
    } else {
        log1mp=log1p(-p);
    }
    return log1mp;
}

#endif /* _RDCAPTURE_CACHE_H */
