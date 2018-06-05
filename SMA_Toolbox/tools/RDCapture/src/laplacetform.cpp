/** @file laplacetform.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-2016
 * @brief Implementation of numerical laplace transform computations
 */

#include "laplacetform.h"
#include <cmath>
#include <cassert>

LaplaceInverseGS::LaplaceInverseGS(int order) : order(order)
{
    computeWeights(order,weights);
}


void LaplaceInverseGS::computeWeights(int order, VecT &weights)
{
    assert(!(order<MinOrder || order>MaxOrder));
    weights.set_size(2*order);
    weights.zeros();
    for(int k=1; k<=2*order; k++) {
        for(int j= floor((k+1)/2); j<=std::min(k,order); j++){
            weights(k-1)+= pow(j,order+1) / (pow(factorial(j),2)*factorial(order-j)*factorial(k-j)*factorial(2*j-k)) * factorial(2*j);
        }
        double sign = pow(-1,order+k);
        weights(k-1) *=sign;
    }
}

