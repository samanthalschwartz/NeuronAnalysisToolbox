/**
 * @file destimator.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @author Peter Relich (physx.grad\@gmail.com)
 * @date 08-05-2014
 * @brief The class declaration for DEstimator.
 *
 * This class is designed to be used either directly from  C++ or indirectly in Matlab via the DEstimator_Iface class interface.
 */
#ifndef _DESTIMATOR_H
#define _DESTIMATOR_H

#include<armadillo>
#include<list>


/**
 * @brief A class for computing the log-likelihood of diffusion constants from observed noisy trajectories
 *
 * This class includes the C++ implementation of the algorithms from "Estimation of the Diffusion Constant from Intermittent Trajectories with
 * Variable Position Uncertainties", Relich et. al. 2015
 *
 * This class is templated on the floating point type so that it may be implemented with either float or double as the primary type.
 *
 * The class is instantiated with the observations, standard error and times for a single observed trajectory.  The trajectory can be
 * of arbitrary dimension Ndim.  The default method used by the class is the recursive algorithm for computing the log-likelihood.  This computation
 * is done in parallel for each separate D value when multiple D values are queried.
 *
 * The class also provides static methods that operate on 1D observed trajectories and implement  the
 * recursive, Laplace, and Markov methods as described in the main text and supplement.  These static
 * member functions can be used directly without instantiating the class.
 *
 * A note on units.  This implementation maintains a generic unit approach.  The observed positions and observation times can use any
 * units and the diffusion constant will then have units: (dist_unit)^2/(time_unit).  For this reason we do not explicitly mention
 * units anywhere else in the documentation and allow the users to use whatever units are the most natural.  Note that with this in mind
 * the exposureT variable should use the same time units as T, and V should have units of (dist_unit)^2.
 *
 */
template<class FloatT>
class DEstimator
{
public:
    typedef arma::Col<FloatT> VecT; //A vector type
    typedef arma::Mat<FloatT> MatT; //A matrix type
    
    int Ndim; //Number of dimensions in trajectory
    int N;    //Number of points in trajectory
    MatT Obs; // size:[N x Ndim] the observed positions of the trajectory
    VecT T;  // size:[N x 1] Time of the observations
    MatT SE;  // size:[N x Ndim] The estimated standard error of each observed position in each dimension
    FloatT exposureT; //Frame exposure time.
    
    /** Create a DEstimator object for a single trajectory in Ndim dimensions.
     * @param[in] Obs - size:[N x Ndim] the observed positions of the trajectory
     * @param[in] T - size:[N x 1] Time of the observations
     * @param[in] SE - size:[N x Ndim] The estimated standard error of each observed position in each dimension
     * @param[in] exposureT - The exposure time over which each observation was collected
     */
    DEstimator(const MatT &Obs, const VecT &T, const MatT &SE, FloatT exposureT);
    

    /** Compute the log-likelihood of the trajectory resulting from a pure diffusion process in Ndim dimensions
     * with diffusion constant D.
     * @param[in] D - The diffusion constant of interest >0
     * @returns The log-likelihood: log(P(O|D))
     */
    FloatT LLH(FloatT D);

    /** Compute the log-likelihood of the trajectory resulting from a pure diffusion process with diffusion constant D,
     * considering only a specific dimension dim out of the Ndim dimensions of the trajectory.
     *
     * Since each dimension is separable, the likelihood can be computed separately.  This is used by the LLH() member function
     * which computes the LLH in each dimension separately, then sums the results.
     * .
     * @param[in] D - The diffusion constant of interest >0
     * @param[in] dim - The index of the dimension to compute for (0<=dim< Ndim)
     * @returns The log-likelihood: log(P(O|D))
     */
    FloatT LLHdim(FloatT D, int dim); //Parallelized over Dim

    /** Parallel computation of the log-likelihood for multiple D values.
     *
     * This is exactly like the non-parallel version except we compute the LLH for multiple D values at once.
     *
     * @param[in] D - A vector of D values to compute the LLH for.  Each D>0.
     * @param[out] LLH - A vector of log-likelihood values for each of the given D values
     */
    void LLH(const VecT &D, VecT &LLH); //Parallelized over D

    /** Parallel computation of the log-likelihood for multiple D values over a single specific dimension.
     *
     * This is exactly like the non-parallel version except we compute the LLH for multiple D values at once.
     * This is used by the LLH() member function which computes the LLH in each dimension separately, then sums the results.
     * @param[in] D - A vector of D values to compute the LLH for.  Each D>0.
     * @param[in] dim - The index of the dimension to compute for (0<=dim< Ndim)
     * @param[out] LLH - A vector of log-likelihood values for each of the given D values
     */
    void LLHdim(const VecT &D, int dim, VecT &LLH); //Parallelized over D

    /** Use the recursive method to compute the log-likelihood of a given 1D trajectory, for a vector of candidate D values.
     * @param[in] D - A vector of D values to compute the LLH for.  Each D>0.
     * @param[in] Obs - A vector of 1D position observations
     * @param[in] T - A vector of observation times
     * @param[in] SE - A vector of estimated standard errors for the 1D observations
     * @param[in] exposureT - The exposure time over which each observation was collected.
     * @param[out] LLH - The computed log-likelihood for each of the given D values
     */
    static void LLH_recursive1D( const VecT &D, const VecT &Obs, const VecT &T,  const VecT &SE, FloatT exposureT, VecT &LLH);

    /** Use the Laplace method to compute the log-likelihood of a given 1D trajectory, for a vector of candidate D values.
     * @param[in] D - A vector of D values to compute the LLH for.  Each D>0.
     * @param[in] Obs - A vector of 1D position observations
     * @param[in] T - A vector of observation times
     * @param[in] SE - A vector of estimated standard errors for the 1D observations
     * @param[in] exposureT - The exposure time over which each observation was collected.
     * @param[out] LLH - The computed log-likelihood for each of the given D values
     */
    static void LLH_laplace1D( const VecT &D, const VecT &Obs, const VecT &T, const VecT &SE, FloatT exposureT, VecT &LLH);

    /** Use the Markov method to compute the log-likelihood of a given 1D trajectory, for a vector of candidate D values.
     * @param[in] D - A vector of D values to compute the LLH for.  Each D>0.
     * @param[in] Obs - A vector of 1D position observations
     * @param[in] T - A vector of observation times
     * @param[in] SE - A vector of estimated standard errors for the 1D observations
     * @param[in] exposureT - The exposure time over which each observation was collected.
     * @param[out] LLH - The computed log-likelihood for each of the given D values
     */
    static void LLH_markov1D( const VecT &D, const VecT &Obs, const VecT &T, const VecT &SE, FloatT exposureT, VecT &LLH);

private:
    static const FloatT min_variance; //The smallest a variance is allowed to get to zero in absolute value in the recursive and Markov methods
    static const FloatT min_laplace_variance; //The smallest a variance is allowed to get to zero in absolute value in the Laplace method
    static const FloatT log2pi; // log(2*pi)
    VecT dT; //length N-1 computed once on initialization
    VecT vD; //Variance due to diffusion (omega in paper) temporary space, computed for each D
    VecT vM; //Variance due to Measurement after motion blur correction (eps in paper), temporary space, computed for each D
    
    /** This is the core of the recursive algorithm that is used internally by the LLH_recursive1D static member function
     * as well as the LLH() member functions.  This code is written to use 1D arrays, and assumes that the variances due to measurement
     * and diffusion have been computed and appropriately normalized away from 0.
     *
     * In this method we have already used a given D value to compute vD and vM with the computeVariance() private static method.  Thus
     * there is no need to take in D explicitly
     *
     * @param[in] N - The number of points in the observed trajectory
     * @param[in] Obs - Length=N: The observed positions in 1D
     * @param[in] T - Length=N: The observation times
     * @param[in] vD - Length=N: The variance due to diffusion (omega in paper) vD>min_variance.
     *                 This should be computed by computeVariance() for a given D.
     * @param[in] vM - Length=N: The variance due to measurement and including the motion-blur correction (epsilon in paper). |vM|>min_variance.
     *                 This should be computed by computeVariance() for a given D.
     * @returns The log-likelihood of D as computed with the recursive method.
     */
    static FloatT coreLLH(int N, const FloatT Obs[], const FloatT T[], const FloatT vD[], const FloatT vM[]);

    /** Computes the variance due to diffusion vD (omega in paper) and the variance due to measurement with motion-blur correction, vM (eps in paper).
     * This is used by all of the computational methods to ensure the vM and vD are computed and bounded the same way for any of th
     * computational methods.
     * @param[in] D - The D value to compute the variances for
     * @param[in] dT - The pointwise difference of T.  These are the dT's between each observation.  Length=N-1.
     * @param[in] SE - The estimated observation standard errors
     * @param[in] exposureT - The exposure time over which each observation was collected.
     * @param[in] min_variance - The smallest any variance value should be in absolute value.  This prevents getting variances arbitrarily close
     *                           to 0 which cause numerical difficulty.
     * @param[out] vD - The computed and bounded variance due to diffusion at for each observation (omega in paper)
     * @param[out] vM - The computed and bounded variance due to measurement with motion-blur correction at for each observation (eps in paper)
     */
    static void computeVariance(FloatT D, const VecT &dT, const VecT &SE, FloatT exposureT, FloatT min_variance, VecT &vD, VecT &vM);

    /** Check that the arguments given are a valid 1D trajectory.  This is run by each of the static methods to avoid problems with
     * bad inputs.  The arguments are the same as any of the LLH_*1D() static methods.
     */
    static void checkLLHargs(const VecT &D, const VecT &Obs, const VecT &dT, const VecT &SE, VecT &LLH);
};


#endif /* _DESTIMATOR_H */
