/** @file estimator.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 01-15-2014
 * @brief 
 */
#ifndef _ESTIMATOR_CPP
#define _ESTIMATOR_CPP

#include <cmath>
#include <functional>
// #include <thread>
#include <boost/thread/thread.hpp>

#include <armadillo>

#include "estimator.h"
#include "display.h"
#include "rng.h"
#include "numerical.h"

#ifdef WIN32
    using namespace boost::chrono;
#else
    using namespace std::chrono;
#endif
/**
 *
 * All models will call for maximization through this virtual function.
 * All non-GPU based maximizers will use this version which spawns threads
 * using a non-virual entry point member function Maximizer::thread_entry.
 * GPU-based maximizers will want to do something custom, so they will declare
 * their own virtual maximize_stack.
 *
 * It is also because of the GPU-based mamixmizers that we are putting initilization,
 * and CRLB/LLH calculations in here even though the Model knows how to do them.
 *
 * We expect that those methods will need to also be paralellized and the GPU will need
 * custom code, and the threaded CPU versions will want to also compute those in parallel,
 * so in order to have a consitent call interface to the Maximizer classes,
 * we put the CRLB/LLH and initialization work within the the maximize_stack method.
 *
 *
 */

template<class Model>
Estimator<Model>::Estimator(Model &model)
    : model(model)
{}

template<class Model>
Estimator<Model>::~Estimator()
{}

/* Just a convenience wrapper that allows a call without a theta_init */
template<class Model>
inline
typename Model::Stencil
Estimator<Model>::estimate(const ImageT &im)
{
    ParamT dummy_theta_init;
    dummy_theta_init.zeros();
   return estimate(im, dummy_theta_init);
}


template<class Model>
inline
typename Model::Stencil
Estimator<Model>::estimate(const ImageT &im, const ParamT &theta_init)
{
    auto start_walltime=ClockT::now();
    Stencil est=compute_estimate(im, theta_init);
    record_walltime(start_walltime, 1);
    return est;
}


/* Just a convenience wrapper that allows a call without a theta_init */
template<class Model>
inline
void Estimator<Model>::estimate(const ImageT &im,
                                ParamT &theta, ParamT &crlb, double &log_likelihood)
{
    ParamT dummy_theta_init;
    dummy_theta_init.zeros();
    estimate(im, dummy_theta_init, theta, crlb, log_likelihood);
}

template<class Model>
inline
void Estimator<Model>::estimate(const ImageT &im, const ParamT &theta_init,
                                ParamT &theta, ParamT &crlb, double &log_likelihood)
{
    auto start_walltime=ClockT::now();
    compute_estimate(im, theta_init, theta, crlb, log_likelihood);
    record_walltime(start_walltime, 1);
}

/* Just a convenience wrapper that allows a call without a theta_init */
template<class Model>
inline
void Estimator<Model>::estimate_debug(const ImageT &im, ParamT &theta, ParamT &crlb, double &llh,
                                      MatT &sequence, VecT &sequence_llh)
{
    ParamT dummy_theta_init;
    dummy_theta_init.zeros();
    estimate_debug(im, dummy_theta_init, theta, crlb, llh, sequence, sequence_llh);
}

template<class Model>
inline
void Estimator<Model>::estimate_debug(const ImageT &im, const ParamT &theta_init, ParamT &theta,
                                      ParamT &crlb, double &llh,
                                      MatT &sequence, VecT &sequence_llh)
{
    auto start_walltime=ClockT::now();
    auto est=compute_estimate_debug(im, theta_init, sequence);
    theta=est.theta;
    crlb=cr_lower_bound(model,est);
    llh=log_likelihood(model,im, est);
    sequence_llh.set_size(sequence.n_cols);
    log_likelihood_stack(model,im,sequence,sequence_llh);
    record_walltime(start_walltime, 1);
}


template<class Model>
inline
double Estimator<Model>::mean_walltime() const
{
    return (num_estimations==0) ? 0. : total_walltime/(double)num_estimations;
}


template<class Model>
StatsT Estimator<Model>::get_stats()
{
    StatsT stats;
    stats["num_estimations"]=num_estimations;
    stats["total_walltime"]=total_walltime;
    stats["mean_walltime"]=mean_walltime();
    return stats;
}

template<class Model>
void Estimator<Model>::clear_stats()
{
    num_estimations=0;
    total_walltime=0.;
}

template<class Model>
std::ostream& operator<<(std::ostream &out, Estimator<Model> &estimator)
{
    auto stats=estimator.get_stats();
    out<<"["<<estimator.name()<<"<"<<estimator.model.name()<<">:";
    for(auto stat: stats) out<<" "<<stat.first<<"="<<stat.second;
    out<<"]";
    return out;
}

template<class Model>
inline
void Estimator<Model>::compute_estimate(const ImageT &im, const ParamT &theta_init,
                                   ParamT &theta, ParamT &crlb, double &llh)
{
    auto est = compute_estimate(im,theta_init);
    crlb = cr_lower_bound(model,est);
    llh = log_likelihood(model,im, est);
    theta = est.theta;
}

/**
 *
 * Estimators that produce a sequence of results (e.g. IterativeEstimators) can override this
 * dummy debug implementation.
 */
template<class Model>
inline
typename Model::Stencil
Estimator<Model>::compute_estimate_debug(const ImageT &im, const ParamT &theta_init, ParamVecT &sequence)
{
    sequence=model.make_param_vec(1);
    auto est=compute_estimate(im,theta_init);
    sequence.col(0)=est.theta;
    return est;
}

template<class Model>
void Estimator<Model>::record_walltime(ClockT::time_point start_walltime, int nimages)
{
    double walltime=duration_cast<duration<double>>(ClockT::now() - start_walltime).count();
    total_walltime+=walltime;
    num_estimations+=nimages;
}



/*
 * 
 * Threaded Estimator
 *
 */ 

template<class Model>
ThreadedEstimator<Model>::ThreadedEstimator(Model &model)
    : Estimator<Model>(model),
      max_threads(boost::thread::hardware_concurrency()),
      num_threads(0),
      thread_walltime(std::vector<double>(max_threads,0.))
{
    char *omp_num_threads_var=getenv("OMP_NUM_THREADS");
    if (omp_num_threads_var) {
        int omp_num_threads=atoi(omp_num_threads_var);
        if (0<omp_num_threads && omp_num_threads<=max_threads) {
            max_threads=omp_num_threads;
        }
    }
}

template<class Model>
inline
void 
ThreadedEstimator<Model>::estimate_stack(const ImageStackT &im,
                                         ParamVecT &theta, ParamVecT &crlb, VecT &log_likelihood)
{
    ParamVecT theta_init;
    estimate_stack(im, theta_init, theta, crlb, log_likelihood);
}

template<class Model>
void 
ThreadedEstimator<Model>::estimate_stack(const ImageStackT &im, const ParamVecT &theta_init,
                                         ParamVecT &theta, ParamVecT &crlb, VecT &log_likelihood)
{
    using std::ref;
    using std::bind;
    using boost::thread;
    auto start_walltime=ClockT::now();
    int min_per_thread=4;
    int nimages = static_cast<int>(im.n_slices);
    //The number of threads we will actually run
    num_threads = std::max(std::min(max_threads, static_cast<int>(floor(nimages/min_per_thread))),1);
    std::vector<thread> threads(num_threads-1);
    for(int i=0; i<num_threads-1; i++) {
        threads[i]=thread(bind(&ThreadedEstimator::thread_entry, this, i,
                               ref(im), ref(theta_init), ref(theta), ref(crlb), ref(log_likelihood)));
    }
    //The main thread assumes the role of num_threads-1
    thread_entry(num_threads-1, im, theta_init, theta, crlb, log_likelihood);
    for(int i=0; i<num_threads-1; i++) threads[i].join();

    this->record_walltime(start_walltime, nimages);
}


template<class Model>
double ThreadedEstimator<Model>::mean_thread_walltime()
{
    double total_thread_time=0.;
    double mean;
    if(num_threads==0)  return this->mean_walltime();

    mtx.lock();
    if(this->num_estimations==0) {
        mean=0;
    } else {
        for(int i=0;i<num_threads; i++) total_thread_time+=thread_walltime[i];
        mean = total_thread_time/static_cast<double>(num_threads);
    }
    mtx.unlock();
    return mean;
}

template<class Model>
StatsT ThreadedEstimator<Model>::get_stats()
{
    auto stats=Estimator<Model>::get_stats();
    double mtwalltime=mean_thread_walltime();
    stats["num_threads"]=num_threads;
    stats["mean_thread_walltime"]=mtwalltime;
    stats["total_thread_walltime"]=mtwalltime*num_threads;
    return stats;
}

template<class Model>
void ThreadedEstimator<Model>::clear_stats()
{
    Estimator<Model>::clear_stats();
    thread_walltime=std::vector<double>(max_threads, 0.0);
}

template<class Model>
int ThreadedEstimator<Model>::thread_start_idx(int nimages, int threadid) const
{
    return ceil(nimages/(double)num_threads*threadid);
}

template<class Model>
int ThreadedEstimator<Model>::thread_stop_idx(int nimages, int threadid) const
{
    return std::min(nimages, (int)std::ceil(nimages/(double)num_threads*(threadid+1)));
}

template<class Model>
void ThreadedEstimator<Model>::thread_maximize_stack(int threadid, const ImageStackT &im, const ParamVecT &theta_init,
                                        ParamVecT &theta, ParamVecT &crlb, VecT &log_likelihood)
{
    int nimages = static_cast<int>(im.n_slices);
    int start = thread_start_idx(nimages, threadid);
    int stop = thread_stop_idx(nimages, threadid);
    auto theta_est = this->model.make_param();
    auto crlb_est = this->model.make_param();
    ParamT init;
    init.zeros();
    for(int n=start; n<stop; n++){
        if(!theta_init.is_empty()) init = theta_init.col(n);
        this->compute_estimate(im.slice(n), init, theta_est, crlb_est, log_likelihood(n));
        theta.col(n) = theta_est;
        crlb.col(n) = crlb_est;
    }
}

/**
 * This is a non-virtual entry point which then calls the virtual function which does the actual
 * maximization.
 *
 */
template<class Model>
void ThreadedEstimator<Model>::thread_entry(int threadid, const ImageStackT &im, const ParamVecT &theta_init,
                               ParamVecT &theta, ParamVecT &crlb, VecT &log_likelihood)
{
    auto start_walltime=ClockT::now();
    thread_maximize_stack(threadid, im, theta_init, theta, crlb, log_likelihood);
    double walltime=duration_cast<duration<double>>(ClockT::now() - start_walltime).count();
    mtx.lock();
    thread_walltime[threadid]+=walltime;
    mtx.unlock();
}


template<class Model>
typename Model::Stencil
HeuristicMLE<Model>::compute_estimate(const ImageT &im, const ParamT &theta_init)
{
    return this->model.initial_theta_estimate(im, theta_init);
}

template<class Model>
StatsT CGaussMLE<Model>::get_stats()
{
    auto stats = ThreadedEstimator<Model>::get_stats();
    stats["max_iterations"] = max_iterations;
    return stats;
}

template<class Model>
typename Model::Stencil
CGaussMLE<Model>::compute_estimate(const ImageT &im, const ParamT &theta_init)
{
    auto crlb = model.make_param();
    auto theta = model.make_param();
    double llh;
    compute_estimate(im, theta_init, theta, crlb, llh);
    return model.make_stencil(theta);
}



/* Iterative Maximizer */

template<class Model>
IterativeMaximizer<Model>::IterativeMaximizer(Model &model, int max_iterations)
    : ThreadedEstimator<Model>(model), 
      max_iterations(max_iterations)
{}

template<class Model>
IterativeMaximizer<Model>::MaximizerData::MaximizerData(const Model &model, const ImageT &im,
                                                const Stencil &s, bool save_seq, int max_seq_len)
    : im(im), 
      grad(model.make_param()), 
      rllh(relative_log_likelihood(model,im,s)), 
      s0(s), 
      current_stencil(true),
      save_seq(save_seq), 
      max_seq_len(max_seq_len)
{
    if (save_seq){
        theta_seq = model.make_param_vec(max_seq_len);
        record_sequence(); //record this initial point
    }
}

template<class Model>
inline
void IterativeMaximizer<Model>::MaximizerData::record_sequence(const ParamT &can_theta)
{
    if(save_seq) {
        if(seq_len>=max_seq_len) {
            std::cout<<"Exceeded MaximizerData sequence limit\n";
            return;
        }
        theta_seq.col(seq_len++) = can_theta;
    }
}


template<class Model>
inline
double IterativeMaximizer<Model>::mean_iterations()
{
    mtx.lock();
    double mean =  num_estimations ? total_iterations/static_cast<double>(num_estimations) : 0;
    mtx.unlock();
    return mean;
}

template<class Model>
inline
double IterativeMaximizer<Model>::mean_backtracks()
{
    mtx.lock();
    double mean =  num_estimations ? total_backtracks/static_cast<double>(num_estimations) : 0;
    mtx.unlock();
    return mean;
}

template<class Model>
StatsT IterativeMaximizer<Model>::get_stats()
{
    auto stats = ThreadedEstimator<Model>::get_stats();
    stats["total_iterations"] = total_iterations;
    stats["total_backtracks"] = total_backtracks;
    stats["mean_iterations"] = mean_iterations();
    stats["mean_backtracks"] = mean_backtracks();
    return stats;
}

template<class Model>
inline
void IterativeMaximizer<Model>::clear_stats()
{
    ThreadedEstimator<Model>::clear_stats();
    total_iterations = 0;
    total_backtracks = 0;
}

template<class Model>
inline
void IterativeMaximizer<Model>::record_iterations(int niters)
{
    mtx.lock();
    total_iterations+=niters;
    mtx.unlock();
}

template<class Model>
inline
void IterativeMaximizer<Model>::record_backtracks(int nbacktracks)
{
    mtx.lock();
    total_backtracks+=nbacktracks;
    mtx.unlock();
}

template<class Model>
bool IterativeMaximizer<Model>::backtrack(MaximizerData &data)
{
    double lambda = 1.0;
    data.save_stencil();
//     if(arma::dot(data.grad, data.step)<=0){
//         //We are maximizing so we should be moving in direction of gradiant not away
//         std::cout<<"ERRRRORRRR: gradient negative.\n";
//         std::cout<<"grad: "<<data.grad.t()<<"\n";
//         std::cout<<"step: "<<data.step.t()<<"\n";
//         std::cout<<"<grad, step>: "<<arma::dot(data.grad, data.step)<<"\n";
//     }
    for(int n=0; n<max_backtracks; n++){
        data.set_stencil(model.make_stencil((data.saved_theta() + lambda*data.step).eval(), false));
        double can_rllh = relative_log_likelihood(model, data.im, data.stencil()); //candidate points log-lh
        bool in_bounds = model.theta_in_bounds((data.saved_theta() + lambda*data.step).eval());
        
//         std::cout<<" [Backtrack:"<<n<<"]\n";
//         std::cout<<" Current Theta: "<<data.saved_theta().t();
//         std::cout<<" Step: "<<data.step.t();
//         std::cout<<" Lambda:"<<lambda<<"\n Scaled Step:"<<(lambda*data.step).t();
//         std::cout<<" Proposed Theta: "<<(data.saved_theta() + lambda*data.step).eval().t();
//         std::cout<<" Proposed Theta Inbounds?: "<<in_bounds<<"\n";
//         std::cout<<" Bounded Proposed Theta: "<<data.theta().t();
//         std::cout<<" CurrentRLLH:"<<data.rllh<<" ProposedRLLH:"<<can_rllh<<"\n";
        double linear_step=lambda*arma::dot(data.grad, data.step); //The amount we would go down if linear in grad
        if (in_bounds) {
            if (can_rllh >= data.rllh + alpha*linear_step) { //Must be a sufficient increase
                //Succes - Found a new point wich is slightly better than where we were before, so update rllh
                data.record_sequence(data.theta());  //We have not changed data.rllh that is still the old rllh
                data.rllh = can_rllh;
                record_backtracks(n);
                data.stencil().compute_derivatives();
                return false; //Tell caller to continue optimizing
            } else if (convergence_test(data)) {
                break;
            }
        }
        double old_lambda = lambda;
        double new_lambda = -.5*lambda*linear_step/(can_rllh - data.rllh - linear_step); //The minimum of a quad approx using cllh and linear_step
        double rho = new_lambda/old_lambda;
        rho = std::min(std::max(rho,0.02),0.25);
        lambda = rho*old_lambda;
//         std::cout<<"***Optimization Step Failure\n\tBacktracking:"<<n+1<<" Old Lambda: "<<old_lambda<<" Proposed lambda:"<<new_lambda<<" Corrected Lambda:"<<lambda<<" linear_step:"<<linear_step<<std::endl;
    }
    record_backtracks(max_backtracks);
    data.restore_stencil();
    return true; //Tell caller to stop and returm the original theta. Backtracking failed to improve it.
}

template<class Model>
bool IterativeMaximizer<Model>::convergence_test(MaximizerData &data)
{
    using arma::norm;
    auto ntheta=data.theta();  //new theta
    auto otheta=data.saved_theta(); //old theta
    double step_size_ratio = norm(otheta-ntheta,2)/std::max(norm(otheta,2),norm(ntheta,2));
    double function_change_ratio = norm(data.grad,2)/fabs(data.rllh);
//     if(step_size_ratio<=delta) std::cout<<"$$$StepSizeRatioTestConvergence! - ratio:"<<step_size_ratio<<" delta:"<<delta<<"\n";
//     if(function_change_ratio<=epsilon) std::cout<<"$$$FunctionChangeRatioTestConvergence! - ratio:"<<function_change_ratio<<" epslion:"<<epsilon<<"\n";
    return step_size_ratio<=delta or function_change_ratio<=epsilon;
}


template<class Model>
typename Model::Stencil
IterativeMaximizer<Model>::compute_estimate(const ImageT &im, const ParamT &theta_init)
{
    auto theta_init_stencil = this->model.initial_theta_estimate(im, theta_init);
    assert(theta_init_stencil.derivatives_computed);
    MaximizerData data(model, im, theta_init_stencil);
    maximize(data);
    return data.stencil();
}

template<class Model>
typename Model::Stencil
IterativeMaximizer<Model>::compute_estimate_debug(const ImageT &im, const ParamT &theta_init, ParamVecT &sequence)
{
    auto theta_init_stencil = this->model.initial_theta_estimate(im, theta_init);
    assert(theta_init_stencil.derivatives_computed);
    MaximizerData data(model, im, theta_init_stencil, true, max_iterations*max_backtracks+1);
    maximize(data);
    sequence = data.theta_sequence();
    return data.stencil();
}

/* This is called to clean up simulated annealing */
template<class Model>
void IterativeMaximizer<Model>::local_maximize(const ImageT &im, const Stencil &theta_init, Stencil &stencil, double &rllh)
{
    MaximizerData data(model, im, theta_init);
    maximize(data);
    stencil = data.stencil();
    rllh = data.rllh;
}

template<class Model>
void NewtonRaphsonMLE<Model>::maximize(MaximizerData &data)
{
    auto grad2 = model.make_param();
    int n;
    for(n=0; n<max_iterations; n++) { //Main optimization loop
        model_grad2(model, data.im, data.stencil(), data.grad, grad2);
        data.step = -data.grad/grad2;
        assert(data.step.is_finite());
        if(backtrack(data) || convergence_test(data)) {
            //Backing up did not help.  Just quit.
            record_iterations(n+1);
            return;
        }
    }
    record_iterations(n);
}


template<class Model>
void NewtonMLE<Model>::maximize(MaximizerData &data)
{
    auto hess = model.make_param_mat();
    auto C = model.make_param_mat();
    auto Cstep = model.make_param();
    int n;
    for(n=0; n<max_iterations; n++) { //Main optimization loop
        model_hessian(model, data.im, data.stencil(), data.grad, hess);
        data.step = arma::solve(arma::symmatu(hess), -data.grad);
        C=-hess;
        copy_Usym_mat(C);
//         std::cout<<"{Newton ITER:"<<n<<"}\n";
//         std::cout<<"Theta:\n"<<data.theta().t()<<"\n";
//         std::cout<<"RLLH: "<<relative_log_likelihood(model, data.im, data.stencil())<<"\n";
//         std::cout<<"LLH: "<<log_likelihood(model, data.im, data.stencil())<<"\n";
//         std::cout<<"Hess:\n"<<arma::symmatu(hess)<<"\n";
//         std::cout<<"Positive-definite:"<<is_positive_definite(C)<<"\n";
        modified_cholesky(C);
        Cstep = cholesky_solve(C,data.grad);
        auto Cfull = C;
        cholesky_convert_full_matrix(Cfull);
//         std::cout<<"C:\n"<<Cfull<<"\n";
//         std::cout<<"Grad: "<<data.grad.t()<<"\n";
//         std::cout<<"HessStep: "<<data.step.t()<<"\n";
//         std::cout<<"HessStepDir: "<<arma::dot(data.step,data.grad)<<"\n";
//         std::cout<<"CholStep: "<<Cstep.t()<<"\n";
//         std::cout<<"CStepDir: "<<arma::dot(Cstep,data.grad)<<"\n";
        data.step = Cstep;
        assert(data.step.is_finite());
        assert(arma::dot(data.grad, data.step)>0);
        if(backtrack(data) || convergence_test(data)) {
            //Backing up did not help.  Just quit.
            record_iterations(n+1);
            return;
        }
    }
    record_iterations(n);
}


template<class Model>
void QuasiNewtonMLE<Model>::maximize(MaximizerData &data)
{
    auto grad_old=model.make_param();
    auto H=model.make_param_mat(); //This is our approximate hessian
    int n;
    for(n=0; n<max_iterations; n++) { //Main optimization loop
        if(n==0) {
            auto hess=model.make_param_mat();
            model_hessian(model, data.im, data.stencil(), data.grad, hess);
            H=arma::inv(arma::symmatu(hess));
        } else {
            //Approx H
            model_grad(model, data.im, data.stencil(), data.grad);
            ParamT delta_grad=data.grad-grad_old;
            double rho=1./arma::dot(delta_grad, data.step);
            ParamMatT K=rho*data.step*delta_grad.t();//This is our approximate inverse hessian
            K.diag()-=1;
            H=K*H*K.t()+rho*data.step*data.step.t();
        }
        data.step=-H*data.grad;
        assert(data.step.is_finite());
        grad_old=data.grad;
        if(backtrack(data) || convergence_test(data)) {
            //Backing up did not help.  Just quit.
            record_iterations(n+1);
            return;
        }
        data.step=data.theta()-data.saved_theta();//If we backtracked then the step may have changed
    }
    record_iterations(n);
}

template<class Model>
typename Model::Stencil
SimulatedAnnealingMLE<Model>::compute_estimate(const ImageT &im, const ParamT &theta_init)
{
    auto rng = make_parallel_rng_stream(make_seed());
    ParamVecT sequence;
    auto theta_init_stencil = model.initial_theta_estimate(im,theta_init);
    return anneal(rng, im, theta_init_stencil, sequence);
}

template<class Model>
typename Model::Stencil
SimulatedAnnealingMLE<Model>::compute_estimate_debug(const ImageT &im, const ParamT &theta_init, ParamVecT &sequence)
{
    auto rng = make_parallel_rng_stream(make_seed());
    auto theta_init_stencil = model.initial_theta_estimate(im,theta_init);
    return anneal(rng, im, theta_init_stencil, sequence);
}


template<class Model>
typename Model::Stencil
SimulatedAnnealingMLE<Model>::anneal(RNG &rng, const ImageT &im, Stencil &theta_init, ParamVecT &sequence)
{
    NewtonRaphsonMLE<Model> nr(model);
    UnitRNG u;
    int niters=max_iterations*model.num_candidate_sampling_phases;
    sequence=model.make_param_vec(niters+1);
    VecT sequence_rllh(niters+1);
    sequence.col(0)=theta_init.theta;
    sequence_rllh(0)=relative_log_likelihood(model, im, theta_init);
    double max_rllh=sequence_rllh(0);
    int max_idx=0;
    Stencil max_s;
    double T=T_init;
    int naccepted=1;
    for(int n=1; n<niters; n++){
        ParamT can_theta=sequence.col(naccepted-1);
        model.sample_candidate_theta(n, rng, can_theta);
        if(!model.theta_in_bounds(can_theta)) { //OOB
            n--;
            continue;
        }
        double can_rllh=relative_log_likelihood(model, im, can_theta);
        assert(std::isfinite(can_rllh));
        double old_rllh=sequence_rllh(naccepted-1);
        if(can_rllh < old_rllh && u(rng)>exp((can_rllh-old_rllh)/T) ){
            continue;//Reject
        }
        //Accept
        T/=cooling_rate;
        sequence.col(naccepted)=can_theta;
        sequence_rllh(naccepted)=can_rllh;
        if(can_rllh>max_rllh){
            max_rllh=can_rllh;
            max_idx=naccepted;
        }
        naccepted++;
        assert(naccepted<niters+1);
    }

    //Run a NR maximization
    nr.local_maximize(im, model.make_stencil(sequence.col(max_idx)), max_s, max_rllh);
    //Fixup sequence to return
    sequence.resize(sequence.n_rows, naccepted+1);
    sequence.col(naccepted)=max_s.theta;
    return max_s;
}


#endif /* _ESTIMATOR_CPP */
