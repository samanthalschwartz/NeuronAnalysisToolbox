/** @file KDTree_Iface.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 09-2014
 * @brief The class declaration and inline and templated functions for KDTree_Iface.
 */

#ifndef _KDTree_IFACE
#define _KDTree_IFACE

#include "Mex_Iface.h"
#include "kdtree.h"

template<class CoordT>
class KDTree_Iface : public Mex_Iface
{
public:
    KDTree_Iface();
private:
    KDTree<CoordT> *obj;
    void objConstruct();
    void objDestroy();
    void getObjectFromHandle(const mxArray *mxhandle);

    void objQuery();
//     void objStats();
};


template<class CoordT>
KDTree_Iface<CoordT>::KDTree_Iface() 
    : Mex_Iface("KDTree_Iface") 
{
    methodmap["query"]=boost::bind(&KDTree_Iface::objQuery, this);
//     methodmap["stats"]=boost::bind(&KDTree_Iface::objStats, this);
}

template<class CoordT>
void KDTree_Iface<CoordT>::objConstruct() 
{
    // [in] points: NxM set of M different N-dimensional points
    checkNumArgs(1,1);
    auto points=getMat<CoordT>();
    auto *obj=new KDTree<CoordT>(points);
    outputMXArray(Handle<KDTree<CoordT>>::makeHandle(obj));
}

template<class CoordT>
void KDTree_Iface<CoordT>::objDestroy()
{
    checkNumArgs(0,1);
    Handle<KDTree<CoordT>>::destroyObject(rhs[0]);
}

template<class CoordT>
void KDTree_Iface<CoordT>::getObjectFromHandle(const mxArray *mxhandle)
{
    obj=Handle<KDTree<CoordT>>::getObject(mxhandle);
}

template<class CoordT>
void KDTree_Iface<CoordT>::objQuery()
{
    // [in] min_corner: vec of minimum corner of box to query
    // [in] max_corner: vec of maximum corner of box to query
    // [out] points: a N x Q matrix of points
    checkNumArgs(1,2);
    auto min=getVec<CoordT>();
    auto max=getVec<CoordT>();
    auto points=obj->query(min,max);
    outputMat<CoordT>(points);
}

// template<class CoordT>
// void KDTree_Iface<CoordT>::objStats()
// {
//     checkNumArgs(1,0);
//     outputStatsToStruct(obj->stats());
// }

#endif /* _KDTree_IFACE */
