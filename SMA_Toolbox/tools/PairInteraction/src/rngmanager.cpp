/** @file rngmanager.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief Fast auto rng for parallel openmp code
 */
#include "rngmanager.h"

namespace pair_int{
    RNGManager<std::mt19937_64> RNG;//The default random number generator
}

