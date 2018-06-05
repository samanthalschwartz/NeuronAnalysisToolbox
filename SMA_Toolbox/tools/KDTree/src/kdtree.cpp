/**
 * @file kdtree.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-05-2014
 * @brief The class method definitions for KDTree.
 */
#include <cassert>
#include <stdexcept>
#include "kdtree.h"

template<class CoordT>
OrthagonalRange<CoordT>::OrthagonalRange(const PointT &_min_corner,const PointT &_max_corner)
    : dim(_min_corner.n_elem), isEmpty(false), min_corner(_min_corner), max_corner(_max_corner)
{
    if (dim<1) throw std::invalid_argument("Dimension must be >=1");
    if (static_cast<int>(max_corner.n_elem) != dim) throw std::invalid_argument("Dimensions of corners do not match");
    for(int d=0; d<dim; d++){
        if (min_corner(d)>max_corner(d)){
            isEmpty=true;
            break;
        }
    }
}
/**
 *
 *
 *
 */
template<class CoordT>
void OrthagonalRange<CoordT>::trimLeft(int d, CoordT axis)
{
    if( isEmpty ) return;
    if( axis<min_corner(d) ) isEmpty=true;
    if( axis>max_corner(d) ) return;
    max_corner(d)=axis;
}

template<class CoordT>
void OrthagonalRange<CoordT>::trimRight(int d, CoordT axis)
{
    if( isEmpty ) return;
    if( axis>max_corner(d) ) isEmpty=true;
    if( axis<min_corner(d) ) return;
    min_corner(d)=axis;
}

template<class CoordT>
bool OrthagonalRange<CoordT>::isLeftOf(int d, CoordT axis) const
{
    if( isEmpty ) return false;
    return axis > min_corner(d);
}

template<class CoordT>
bool OrthagonalRange<CoordT>::isRightOf(int d, CoordT axis) const
{
    if( isEmpty ) return false;
    return axis <= max_corner(d); //Tree should be right biased wher >= axis are on right and < on left
}

template<class CoordT>
bool OrthagonalRange<CoordT>::contains(const PointT &p) const
{
    for(int d=0; d<dim; d++)
        if(p(d)<min_corner(d) || p(d)>max_corner(d)) return false;
    return true;
}

template<class CoordT>
bool OrthagonalRange<CoordT>::contains(const OrthagonalRange &o) const
{
    for(int d=0; d<dim; d++)
        if(min_corner(d)>=o.min_corner(d) || max_corner(d)<o.max_corner(d)) return false;
    return true;
}

template<class CoordT>
bool OrthagonalRange<CoordT>::intersects(const OrthagonalRange &o) const
{
    for(int d=0; d<dim; d++){
        if(min_corner(d)<o.min_corner(d)){
            if(max_corner(d)<o.min_corner(d)) return false;
        } else {
            if(o.max_corner(d)<min_corner(d)) return false;
        }
    }
    return true;
}

template<class CoordT>
const typename KDTree<CoordT>::IndexT
KDTree<CoordT>::EMPTY=std::numeric_limits<typename KDTree<CoordT>::IndexT>::max();

template<class CoordT>
KDTree<CoordT>::KDTree(const PointVecT &points_)
    : points(points_),
      dim(points.n_rows), 
      Npoints(points_.n_cols),
      left(Npoints), right(Npoints), 
      root(-1),
      point_range(arma::min(points_,1),arma::max(points_,1))
{
    assert(point_range.dim==dim);
    left.fill(-1);
    right.fill(-1);
    IndexVecT idx(Npoints);
    for(IndexT i=0;i<Npoints;i++) idx[i]=i;
    root=build(idx, 0);
}

template<class CoordT>
typename KDTree<CoordT>::IndexT
KDTree<CoordT>::build(const IndexVecT &idxs, int depth)
{
    IndexT split_dim=depth%dim;
    IndexT N=idxs.n_elem;
    if(N==1) return idxs(0);


    IndexVecT sorted_indexs=arma::sort_index(points.row(split_dim).eval()(idxs));
    IndexT median_sort_idx=N/2;
    while( median_sort_idx>=1 &&
            points(split_dim,idxs(sorted_indexs(median_sort_idx-1))) == points(split_dim,idxs(sorted_indexs(median_sort_idx)))){
        median_sort_idx--;
    }
    IndexT median_idx=idxs(sorted_indexs(median_sort_idx));

    if(median_sort_idx>0) left(median_idx)=build(idxs(sorted_indexs.subvec(0,median_sort_idx-1)),depth+1);
    else left(median_idx)=EMPTY;
    if(median_sort_idx<N-1) right(median_idx)=build(idxs(sorted_indexs.subvec(median_sort_idx+1,N-1)),depth+1);
    else right(median_idx)=EMPTY;
    return median_idx;
}

template<class CoordT>
typename KDTree<CoordT>::PointVecT
KDTree<CoordT>::query(const PointT &corner0, const PointT &corner1) const
{
    if(Npoints==0) return PointVecT();
    RangeT Q(corner0,corner1);
    IndexListT idxs=queryPoints(Q,point_range,root,0);
    IndexVecT I(idxs.size());
    int count=0;
    for(auto i :idxs) I(count++)=i;
    return points.cols(I);
}

template<class CoordT>
typename KDTree<CoordT>::IndexListT
KDTree<CoordT>::queryPoints(const RangeT &Q, const RangeT &C, IndexT node, int depth) const
{
    assert(not C.empty());
    assert(node<Npoints);
    if(not Q.intersects(C)) return IndexListT();
    if(Q.empty()) return IndexListT();
    if(Q.contains(C)) return allPoints(node);
    IndexListT idxs;
    if(Q.contains(points.col(node))) idxs.push_back(node);
    IndexT split_dim=depth%dim;
    CoordT x=points(split_dim,node);
    IndexT leftN=left(node);
    IndexT rightN=right(node);
    if(leftN!=EMPTY && Q.isLeftOf(split_dim,x))
        idxs.splice(idxs.end(),queryPoints(Q,C.trimedLeft(split_dim,x), leftN, depth+1));
    if(rightN!=EMPTY && Q.isRightOf(split_dim,x))
        idxs.splice(idxs.end(),queryPoints(Q,C.trimedRight(split_dim,x), rightN, depth+1));
    return idxs;
}

template<class CoordT>
typename KDTree<CoordT>::IndexListT
KDTree<CoordT>::allPoints(IndexT node) const
{
    IndexListT idxs;
    idxs.push_back(node);
    IndexT leftN=left(node);
    IndexT rightN=right(node);
    if(leftN!=EMPTY) idxs.splice(idxs.end(),allPoints(left(node)));
    if(rightN!=EMPTY) idxs.splice(idxs.end(),allPoints(right(node)));
    return idxs;
}

/* Explicit Template Instantiation */
template class KDTree<int32_t>;
template class KDTree<uint32_t>;
template class KDTree<float>;
template class KDTree<double>;

template class OrthagonalRange<int32_t>;
template class OrthagonalRange<uint32_t>;
template class OrthagonalRange<float>;
template class OrthagonalRange<double>;


