/**
 * @file kdtree.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-05-2014
 * @brief The class declaration for KDTree.
 * 
 * 
 */
#ifndef _KDTREE_H
#define _KDTREE_H

#include <armadillo>
#include <cstdint>
#include <map>
#include <list>
#include <string>

template<class CoordT>
class OrthagonalRange
{
public:
    typedef arma::Col<CoordT> PointT;
    OrthagonalRange(const PointT &min_corner,const PointT &max_corner);
    bool contains(const PointT &p) const;
    bool contains(const OrthagonalRange &o) const;
    bool intersects(const OrthagonalRange &o) const;
    void trimLeft(int d, CoordT axis);
    void trimRight(int d, CoordT axis);
    OrthagonalRange<CoordT> trimedLeft(int d, CoordT axis) const;
    OrthagonalRange<CoordT> trimedRight(int d, CoordT axis) const;
    bool isLeftOf(int d, CoordT axis) const;
    bool isRightOf(int d, CoordT axis) const;
    template<class C> friend std::ostream& operator<<(std::ostream &out, const OrthagonalRange<C> &o);
    bool empty() const;
    const int dim;
private:
    bool isEmpty;
    PointT min_corner, max_corner;
};



template<class CoordT>
class KDTree
{
public:
    typedef uint32_t IndexT;
    typedef arma::Col<IndexT> IndexVecT;
    typedef arma::Col<CoordT> CoordVecT;
    typedef arma::Col<CoordT> PointT;
    typedef arma::Mat<CoordT> PointVecT;
    typedef std::list<IndexT> IndexListT;
    typedef OrthagonalRange<CoordT> RangeT;


    typedef std::map<std::string,double> StatsT;

    KDTree(const PointVecT &points);
    PointVecT query(const PointT &corner0, const PointT &corner1) const;
    StatsT stats() const;
private:
    static const IndexT EMPTY;
    PointVecT points; // dim X Npoints
    int dim; //Dimensionality of points
    IndexT Npoints; //Number of points
    IndexVecT left; //Next node left from current node [intmax=Empty]
    IndexVecT right;//Next node right from current node [intmax=Empty]
    IndexT root; //Root node (median of dim=0)
    RangeT point_range; //Bounding box for all points in tree
    IndexListT queryPoints(const RangeT &Q, const RangeT &C, IndexT node, int depth) const;
    IndexListT allPoints(IndexT node) const;

    IndexT build(const IndexVecT &indexs, int depth);
};

typedef KDTree<int32_t> KDTreeI;
typedef KDTree<uint32_t> KDTreeU;
typedef KDTree<float> KDTreeF;
typedef KDTree<double> KDTreeD;

template<class CoordT>
inline
bool OrthagonalRange<CoordT>::empty() const
{
    return isEmpty;
}

template<class CoordT>
OrthagonalRange<CoordT>
OrthagonalRange<CoordT>::trimedLeft(int d, CoordT axis) const
{
    OrthagonalRange<CoordT> o(*this);
    o.trimLeft(d,axis);
    return o;
}

template<class CoordT>
OrthagonalRange<CoordT>
OrthagonalRange<CoordT>::trimedRight(int d, CoordT axis) const
{
    OrthagonalRange<CoordT> o(*this);
    o.trimRight(d,axis);
    return o;
}

template<class CoordT>
std::ostream& operator<<(std::ostream &out, const OrthagonalRange<CoordT> &o)
{
    out<<"[Empty:"<<o.isEmpty<<" (";
    for(int d=0;d<o.dim-1;d++) out<<o.min_corner(d)<<",";
    out<<o.min_corner(o.dim-1);
    out<<") --- (";
    for(int d=0;d<o.dim-1;d++) out<<o.max_corner(d)<<",";
    out<<o.max_corner(o.dim-1);
    out<<")]";
    return out;
}


#endif /* _KDTREE_H */
