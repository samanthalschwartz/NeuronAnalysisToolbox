#include "srrender.h"
#include <random>

using namespace arma;
using namespace std;

void test2D()
{
    vec roi={0., 256., 0., 256.};
    int nPoints=1e5;
    default_random_engine generator;
    uniform_real_distribution<double> xDist(0.,256.);
    uniform_real_distribution<double> yDist(0.,256.);
    mat points(nPoints,5);
    for(int n=0; n<nPoints; n++){
        points(n,0)=1.;
        points(n,1)=xDist(generator);
        points(n,2)=yDist(generator);
        points(n,3)=0.5;
        points(n,4)=0.5;
    }
    typename SRRender2D<double>::ImageT im(2048,2048);
    SRRender2D<double>::renderGauss(points, roi, im, 5.);
}


int main(){
    test2D();
    return 0;
}
