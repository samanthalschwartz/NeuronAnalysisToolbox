/**
 * @file test_kdtree.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-05-2014
 */

#include "kdtree.h"

using namespace arma;
using namespace std;

int num_contained(mat &points, OrthagonalRange<double> &range)
{
    int count=0;
    for(int n=0; n<(int)points.n_cols; n++) count+=range.contains(points.col(n));
    return count;
}


int main(){
    mat points(2,100,arma::fill::randu);
    points*=10;
    KDTreeD kdt(points);

    vec min_corner1 = {0, 0, 0};
    vec max_corner1 = {2, 2, 2};

    vec min_corner2 = {0, 0, 0};
    vec max_corner2 = {1, 2, 2};

    OrthagonalRange<double> b1(min_corner1, max_corner1);
    OrthagonalRange<double> b2(min_corner2, max_corner2);
    cout<<"B1:"<<b1<<" B2:"<<b2<<endl;
    cout<<"b1 contains b2: "<<b1.contains(b2)<<endl;
    cout<<"b2 contains b1: "<<b2.contains(b1)<<endl;
    cout<<"b1 contains b1: "<<b1.contains(b1)<<endl;
    cout<<"b2 contains b2: "<<b2.contains(b2)<<endl;
    cout<<"b1 intersects b2: "<<b1.intersects(b2)<<endl;
    cout<<"b2 intersects b1: "<<b2.intersects(b1)<<endl;
    cout<<"b1 intersects b1: "<<b1.intersects(b1)<<endl;
    cout<<"b2 intersects b2: "<<b2.intersects(b2)<<endl;
    for(int d=0;d<3;d++){
        b1.trimLeft(d,1);
        b2.trimRight(d,1);
        cout<<"B1:"<<b1<<" B2:"<<b2<<endl;
    }
    cout<<"b1 intersects b2: "<<b1.intersects(b2)<<endl;
    cout<<"b2 intersects b1: "<<b2.intersects(b1)<<endl;
    
    vec minc={3,7};
    vec maxc={3.3,7.3};
    OrthagonalRange<double> Q(minc,maxc);
    mat query=kdt.query(minc,maxc);
    cout<<"Query: "<<Q<<" -->"
        <<"\n  count: "<<query.n_cols
        <<"\n  true_count: "<<num_contained(points,Q)<<"\n";
    if(query.n_cols>0){
        cout<<"\n  min: "<<arma::min(query,1).t();
        cout<<"\n  max: "<<arma::max(query,1).t();
    }
    return 1;
}
