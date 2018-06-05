/** @file rngmanager.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 12-2016
 * @brief Fast auto rng for parallel openmp code
 */

#ifndef _RNGMANAGER_H
#define _RNGMANAGER_H

#include <random>
#include <armadillo>

namespace pair_int {

    
template<typename RngT>
class RNGManager
{
public:
    using IdxT = arma::uword;
    RNGManager();
    RNGManager(IdxT _max_threads);
    RngT& generator();
    double randu();
    double randn();
    arma::vec randu(IdxT N);
    arma::vec randn(IdxT N);
    arma::mat randu(IdxT rows, IdxT cols);
    arma::mat randn(IdxT rows, IdxT cols);
    inline RngT& operator()(){return generator();}
    IdxT resample_dist(const arma::vec &weights);
    arma::uvec resample_dist(const arma::vec &weights, IdxT N);
protected:
    int max_threads;
    std::vector<RngT> thread_local_rng;
    void initializeThreadRNGs(IdxT num_threads);
    std::normal_distribution<double> norm_dist;
    std::uniform_real_distribution<double> uni_dist;
};

extern RNGManager<std::mt19937_64> RNG;



template<typename RngT>
RNGManager<RngT>::RNGManager()
{
    max_threads = omp_get_max_threads();
    initializeThreadRNGs(max_threads);
}

template<typename RngT>
RNGManager<RngT>::RNGManager(IdxT _max_threads) : max_threads(_max_threads)
{
    initializeThreadRNGs(max_threads);
}


template<typename RngT>
void RNGManager<RngT>::initializeThreadRNGs(IdxT num_threads)
{
    thread_local_rng.clear();
    std::random_device true_rnd;
    for(IdxT n=0; n<num_threads; n++) thread_local_rng.emplace_back(true_rnd());
}

template<typename RngT>
inline
RngT& RNGManager<RngT>::generator()
{
    return thread_local_rng[omp_get_thread_num()];
}

template<typename RngT>
inline
double RNGManager<RngT>::randu()
{
    return uni_dist(generator());
}

template<typename RngT>
inline
double RNGManager<RngT>::randn()
{
    return norm_dist(generator());
}

template<typename RngT>
arma::vec RNGManager<RngT>::randu(IdxT N)
{
    arma::vec samp(N);
    RngT &gen = generator();
    for(IdxT n=0;n<N;n++) samp(n) = uni_dist(gen);
    return samp;
}

template<typename RngT>
arma::vec RNGManager<RngT>::randn(IdxT N)
{
    arma::vec samp(N);
    RngT &gen = generator();
    for(IdxT n=0;n<N;n++) samp(n) = norm_dist(gen);
    return samp;
}

template<typename RngT>
arma::mat RNGManager<RngT>::randu(IdxT rows, IdxT cols)
{
    arma::mat samp(rows, cols);
    RngT &gen = generator();
    for(IdxT r=0;r<rows;r++) for(IdxT c=0;c<cols;c++) samp(r,c) = uni_dist(gen);
    return samp;
}

template<typename RngT>
arma::mat RNGManager<RngT>::randn(IdxT rows, IdxT cols)
{
    arma::mat samp(rows, cols);
    RngT &gen = generator();
    for(IdxT r=0;r<rows;r++) for(IdxT c=0;c<cols;c++) samp(r,c) = norm_dist(gen);
    return samp;
}

template<typename RngT>
arma::uword
RNGManager<RngT>::resample_dist(const arma::vec &weights)
{
    std::discrete_distribution<arma::uword> dist(weights.begin(),weights.end());
    return dist(generator());
}

template<typename RngT>
arma::uvec
RNGManager<RngT>::resample_dist(const arma::vec &weights, arma::uword N)
{
    std::discrete_distribution<arma::uword> dist(weights.begin(),weights.end());
    arma::uvec samp(N);
    RngT g = generator();
    for(IdxT n=0; n<N; n++) samp(n) = dist(g);
    return samp;
}




} /* namespace */
#endif /* _RNGMANAGER_H */
