/** @file HSData_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-21-2014
 * @brief The class declaration and inline and templated functions for HSData_Iface.
 */

#ifndef _HSDATA_IFACE
#define _HSDATA_IFACE
#include "Mex_Iface.h"

#include "HSData.h"

class HSData_Iface : public Mex_Iface
{
public:
    typedef HSData::PixelT PixelT;
    typedef HSData::RawPixelT RawPixelT;
    HSData_Iface();
private:
    HSData *obj;
    
    void callMethod(std::string name);
    void objConstruct();
    void objDestroy();
    void objLoadFrames();
    void objLoadRawFrames();
    void getObjectFromHandle(const mxArray *mxhandle);
    
    
};


#endif /* _HSDATA_IFACE */
