/** @file estimator.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 04-01-2014
 * @brief The class declaration and inline and templated functions for the 
 * Estimator class hierarchy.
 */
#ifndef _ESTIMATOR_H
#define _ESTIMATOR_H

#include <fstream>
#include <string>
#include <limits>
#include <memory>
#include <map>
#include "rng.h"

#include <boost/thread/mutex.hpp>


// TODO fix this once is there is some progress in x-platform timing libraries
// Std chrono works in linux, but boost chrono has a bug somewhere as of 1.52
// Std chrono does not work well in linux and liekly will not cause of GCC
#ifdef WIN32
    #include <boost/chrono.hpp>
    typedef boost::chrono::high_resolution_clock ClockT;
#else
    #include <chrono>
    typedef std::chrono::high_resolution_clock ClockT;
#endif

#include "util.h"

#define DEFAULT_ITERATIONS 2000
#define DEFAULT_CGAUSS_ITERATIONS 50

typedef std::map<std::string,double> StatsT;


template<class Model>
class Estimator{
public:
    /* These improve readabilit, but are (unfortunately) not inherited. */
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ParamVecT ParamVecT;
    typedef typename Model::ImageT ImageT;
    typedef typename Model::ImageStackT ImageStackT;
    Model &model;

    Estimator(Model &model); //Why can't we make this const?
    virtual ~Estimator();

    virtual std::string name() const =0;

    /* These are the main entrypoints for estimation
     * We provide three types of interfaces:
     */
    /* Option 1: Single Image - Estimator returns a stencil at the optimium point.
     * Variants with and without a theta_init.
     */
    Stencil estimate(const ImageT &im);
    Stencil estimate(const ImageT &im, const ParamT &theta_init); //with theta_init

    /* Option 2: Single Image - Estimator returns theta, crlb and llh  at the optimium point.
     * Variants with and without a theta_init.
     */
    void estimate(const ImageT &im, ParamT &theta, ParamT &crlb, double &llh);
    void estimate(const ImageT &im, const ParamT &theta_init, ParamT &theta, ParamT &crlb, double &llh);//with theta_init

    /* Option 3: Single Image Debug mode - Estimator returns theta, crlb and llh  at the optimium point.
     * Variants with and without a theta_init.
     */
    void estimate_debug(const ImageT &im, ParamT &theta, ParamT &crlb, double &llh,
                                MatT &sequence, VecT &sequence_llh);
    void estimate_debug(const ImageT &im, const ParamT &theta_init, ParamT &theta, ParamT &crlb, double &llh,
                        MatT &sequence, VecT &sequence_llh);

    /* Option 3: Parallel Image - Estimator returns theta_stack, crlb_stack and llh_stack at the optimium point
     * for each image in the stack.
     * Variants with and without a theta_init.
     */
    virtual void estimate_stack(const ImageStackT &im, const ParamVecT &theta_init, ParamVecT &theta, ParamVecT &crlb, VecT &llh)=0;
    virtual void estimate_stack(const ImageStackT &im, ParamVecT &theta, ParamVecT &crlb, VecT &llh)=0;

    /* Statistics */
    double mean_walltime() const;
    virtual StatsT get_stats();
    virtual void clear_stats();

    /* I/O */
    template<class T>
    friend std::ostream& operator<<(std::ostream &out, Estimator<T> &estimator);

protected:
    virtual Stencil compute_estimate(const ImageT &im, const ParamT &theta_init)=0;
    virtual Stencil compute_estimate_debug(const ImageT &im, const ParamT &theta_init, ParamVecT &sequence);
    virtual void compute_estimate(const ImageT &im, const ParamT &theta_init, ParamT &theta, ParamT &crlb, double &llh);

    /* statistics */
    int num_estimations=0;
    double total_walltime=0.;

    void record_walltime(ClockT::time_point start_walltime, int nimages);
};

/**
 * We avoid combining Estimator and ThreadedEstimator classes so that a future GPU implementation can
 * inherit directly from Estimator as it will present a differnt method for estimate_stack pure virtual
 * member function.  For now all other (CPU) estimators inherit from ThreadedEstimator.
 *
 */
template<class Model>
class ThreadedEstimator : public Estimator<Model> {
public:
    /* These improve readabilit, but are (unfortunately) not inherited. */
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ParamVecT ParamVecT;
    typedef typename Model::ImageT ImageT;
    typedef typename Model::ImageStackT ImageStackT;

    ThreadedEstimator(Model &model);

    void estimate_stack(const ImageStackT &im, const ParamVecT &theta_init,ParamVecT &theta, ParamVecT &crlb, VecT &llh);
    void estimate_stack(const ImageStackT &im, ParamVecT &theta, ParamVecT &crlb, VecT &llh);

    double mean_thread_walltime();
    virtual StatsT get_stats();
    virtual void clear_stats();

protected:
    int max_threads;
    int num_threads;
    std::vector<double> thread_walltime;
    boost::mutex mtx;

    int thread_start_idx(int nimages, int threadid) const;
    int thread_stop_idx(int nimages, int threadid) const;

    virtual void thread_maximize_stack(int thread_id, const ImageStackT &im, const ParamVecT &theta_init,
                                       ParamVecT &theta, ParamVecT &crlb, VecT &llh);

private:
    void thread_entry(int thread_id, const ImageStackT &im, const ParamVecT &theta_init,
                      ParamVecT &theta, ParamVecT &crlb, VecT &llh);
};

// template<class Mode>
// class GPUEstimator<Model> : public Estimator<Model> {
//
//
// }

template<class Model>
class HeuristicMLE : public ThreadedEstimator<Model> {
public:
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ImageT ImageT;
    HeuristicMLE(Model &model) : ThreadedEstimator<Model>(model) {}

    std::string name() const {return "HeuristicMLE";}
private:
    Stencil compute_estimate(const ImageT &im, const ParamT &theta_init);
};


template<class Model>
class CGaussHeuristicMLE : public ThreadedEstimator<Model> {
public:
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ImageT ImageT;
    CGaussHeuristicMLE(Model &model) : ThreadedEstimator<Model>(model) {}
    
    std::string name() const {return "CGaussHeuristicMLE";}
private:
    Stencil compute_estimate(const ImageT &im, const ParamT &theta_init);
};

template<class Model>
class CGaussMLE : public ThreadedEstimator<Model> {
public:
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ImageT ImageT;
    int max_iterations;
    CGaussMLE(Model &model, int max_iterations=DEFAULT_CGAUSS_ITERATIONS)
        : ThreadedEstimator<Model>(model), max_iterations(max_iterations) {}

    StatsT get_stats();
    inline std::string name() const {return "CGaussMLE";}
protected:
    /* These bring in non-depended names from base classes (only necessary because we are templated) */
    using Estimator<Model>::model;

    Stencil compute_estimate(const ImageT &im, const ParamT &theta_init);
    void compute_estimate(const ImageT &im, const ParamT &theta_init, ParamT &theta, ParamT &crlb, double &log_likelihood);
};

template<class Model>
class SimulatedAnnealingMLE : public ThreadedEstimator<Model> {
public:
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ParamVecT ParamVecT;
    typedef typename Model::ImageT ImageT;
    using Estimator<Model>::model;

    double T_init=100.;
    double cooling_rate=1.02;
    double max_iterations=500;

    inline std::string name() const {return "SimulatedAnnealingMLE";}
    SimulatedAnnealingMLE(Model &model) : ThreadedEstimator<Model>(model) {}
protected:
    Stencil compute_estimate(const ImageT &im, const ParamT &theta_init);
    Stencil compute_estimate_debug(const ImageT &im, const ParamT &theta_init, ParamVecT &sequence);
    Stencil anneal(RNG &rng, const ImageT &im, Stencil &theta_init,
                   ParamVecT &sequence);
};


template<class Model>
class IterativeMaximizer : public ThreadedEstimator<Model> {
public:
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ImageT ImageT;
    typedef typename Model::ParamVecT ParamVecT;
    
    IterativeMaximizer(Model &model, int max_iterations=DEFAULT_ITERATIONS);

    /* Statistics */
    double mean_iterations();
    double mean_backtracks();
    StatsT get_stats();
    void clear_stats();

    void local_maximize(const ImageT &im, const Stencil &theta_init, Stencil &stencil, double &rllh); //This is used by SimulatedAnnealing to clean up max

protected:
    using Estimator<Model>::model;
    using Estimator<Model>::num_estimations;
    using ThreadedEstimator<Model>::mtx;
    int max_iterations;

    /* These parameters control the adaptive convergence testing */
    double epsilon=sqrt(std::numeric_limits<double>::epsilon()); //tolerance for fval
    double delta=sqrt(std::numeric_limits<double>::epsilon()); // tolerance for relative step size
    /* These parameters control backtracking */
    double lambda_min=0.05; //What is the minimum proportion of the step to take
    double alpha=1e-4; //How much drop in f-val do we expect for the step to be OK?
    int max_backtracks=12; //Max # of evaluations to do when backtracking

    /* Statistics: need to be private so they can be mutex protected */
    int total_iterations=0;
    int total_backtracks=0;



    struct MaximizerData {
        typedef typename Model::Stencil Stencil;
        typedef typename Model::ParamT ParamT;
        typedef typename Model::ImageT ImageT;
        typedef typename Model::ParamVecT ParamVecT;
        const ImageT &im;
        ParamT grad;
        ParamT step;
        double rllh;

        MaximizerData(const Model &model, const ImageT &im, const Stencil &s,
                      bool save_seq=false, int max_seq_len=0);

        void record_sequence() {record_sequence(theta());}
        void record_sequence(const ParamT &can_theta);

        ParamVecT theta_sequence() const {return theta_seq.cols(0,seq_len-1);}
        
        Stencil& stencil() {return current_stencil ? s0 : s1;}
        void set_stencil(const Stencil &s) {if(current_stencil) s0=s;  else s1=s;}
        void save_stencil() {current_stencil=not current_stencil;}
        void restore_stencil() {current_stencil=not current_stencil;}
        Stencil& saved_stencil() {return current_stencil ? s1 : s0;}
        ParamT& theta() {return current_stencil ? s0.theta : s1.theta;}
        ParamT& saved_theta() {return current_stencil ? s1.theta : s0.theta;}

    private:
        Stencil s0,s1; //These two stencils will be alternated as the current and old stencil points
        bool current_stencil; //This alternates to indicated weather s0 or s1 is the current stencil

        bool save_seq;
        ParamVecT theta_seq;
        int seq_len=0;
        const int max_seq_len;
    };

    void record_iterations(int iters);
    void record_backtracks(int nbacktracks);

    Stencil compute_estimate_debug(const ImageT &im, const ParamT &theta_init, ParamVecT &sequence);
    Stencil compute_estimate(const ImageT &im, const ParamT &theta_init);

    virtual void maximize(MaximizerData &data)=0;
    inline bool backtrack(MaximizerData &data);
    inline bool convergence_test(MaximizerData &data);
};



template<class Model>
class NewtonRaphsonMLE : public IterativeMaximizer<Model> {
public:
     /* These improve readability, but are (unfortunately) not inherited. */
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ImageT ImageT;
    typedef typename Model::ParamVecT ParamVecT;
    typedef typename IterativeMaximizer<Model>::MaximizerData MaximizerData;


    NewtonRaphsonMLE(Model &model, int max_iterations=DEFAULT_ITERATIONS)
        : IterativeMaximizer<Model>(model,max_iterations) {}

    inline std::string name() const {return "NewtonRaphsonMLE";}
protected:
    /* These bring in non-depended names from base classes (only necessary because we are templated) */
    using Estimator<Model>::model;
    using IterativeMaximizer<Model>::record_iterations;
    using IterativeMaximizer<Model>::max_iterations;
    using IterativeMaximizer<Model>::backtrack;
    using IterativeMaximizer<Model>::convergence_test;

    void maximize(MaximizerData &data);

//     friend class SimulatedAnnealingMLE<Model>;
};

template<class Model>
class NewtonMLE : public IterativeMaximizer<Model> {
public:
     /* These improve readability, but are (unfortunately) not inherited. */
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ParamMatT ParamMatT;
    typedef typename Model::ImageT ImageT;
    typedef typename Model::ParamVecT ParamVecT;
    typedef typename IterativeMaximizer<Model>::MaximizerData MaximizerData;


    NewtonMLE(Model &model, int max_iterations=DEFAULT_ITERATIONS)
        : IterativeMaximizer<Model>(model,max_iterations) {}

    inline std::string name() const {return "NewtonMLE";}
protected:
    /* These bring in non-depended names from base classes (only necessary because we are templated) */
    using Estimator<Model>::model;
    using IterativeMaximizer<Model>::record_iterations;
    using IterativeMaximizer<Model>::max_iterations;
    using IterativeMaximizer<Model>::backtrack;
    using IterativeMaximizer<Model>::convergence_test;


    void maximize(MaximizerData &data);
};

template<class Model>
class QuasiNewtonMLE : public IterativeMaximizer<Model> {
public:
     /* These improve readability, but are (unfortunately) not inherited. */
    typedef typename Model::Stencil Stencil;
    typedef typename Model::ParamT ParamT;
    typedef typename Model::ParamMatT ParamMatT;
    typedef typename Model::ImageT ImageT;
    typedef typename Model::ParamVecT ParamVecT;
    typedef typename IterativeMaximizer<Model>::MaximizerData MaximizerData;



    QuasiNewtonMLE(Model &model, int max_iterations=DEFAULT_ITERATIONS)
        : IterativeMaximizer<Model>(model,max_iterations) {}

    inline std::string name() const {return "QuasiNewtonMLE";}
protected:
    /* These bring in non-depended names from base classes (only necessary because we are templated) */
    using Estimator<Model>::model;
    using IterativeMaximizer<Model>::record_iterations;
    using IterativeMaximizer<Model>::max_iterations;
    using IterativeMaximizer<Model>::backtrack;
    using IterativeMaximizer<Model>::convergence_test;


    void maximize(MaximizerData &data);
};

// template<class Model>
// class TrustRegionMLE : public IterativeMaximizer<Model> {
// public:
//     /* These improve readability, but are (unfortunately) not inherited. */
//     typedef typename Model::Stencil Stencil;
//     typedef typename Model::ParamT ParamT;
//     typedef typename Model::ParamMatT ParamMatT;
//     typedef typename Model::ImageT ImageT;
//     typedef typename Model::ParamVecT ParamVecT;
//     typedef typename IterativeMaximizer<Model>::MaximizerData MaximizerData;
//     
//     
//     TrustRegionMLE(Model &model, int max_iterations=DEFAULT_ITERATIONS)
//         : IterativeMaximizer<Model>(model,max_iterations) {}
//     
//     inline std::string name() const {return "TrustRegionMLE";}
// protected:
//     /* These bring in non-depended names from base classes (only necessary because we are templated) */
//     using Estimator<Model>::model;
//     using IterativeMaximizer<Model>::record_iterations;
//     using IterativeMaximizer<Model>::max_iterations;
//     using IterativeMaximizer<Model>::convergence_test;
//     
//     
//     void maximize(MaximizerData &data);
// };
// 


/**
 * Makes a new estimator by name and returns a pointer and ownership.
 * 
 */
// template<class Model>
// std::shared_ptr<Estimator<Model>> make_estimator(Model &model, std::string ename);

#include "estimator.cpp"

#endif /* _ESTIMATOR_H */
