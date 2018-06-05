% File: DEstimator.m
%
%
% Copyright (c) 2015, Mark J. Olah, Peter. K. Relich, Keith A. Lidke,
% Regents of University of New Mexico
% All rights reserved.

% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:

%1. Redistributions of source code must retain the above copyright notice, 
%   this list of conditions and the following disclaimer.

%2. Redistributions in binary form must reproduce the above copyright notice, 
%   this list of conditions and the following disclaimer in the documentation 
%   and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.%

% The DEstimator class is the Matlab interface to the C++ DEstimator class.  It uses the IfaceMixin superclass
% to create a connection to the C++ DEstimator_IFace Mex module which exposes an interface to the C++ class.
% The DEstimator_IFace mex function should not be used directly but rather only through this DEstimator class,
% which ensures that arguments are given to the C++ side in the correct format.


classdef DEstimator < IfaceMixin
    % DEstimator - A class to estimate the Diffusion constant log-likelihood given noisy observations of a particle
    % moving in N dimensions.  Each instance of the DEstimator class works with a single track that is specified by
    % the observation locations, observation times, observation measurement standard error, and observation motion blur
    % coefficient.  An instance of the class can then have log-likelihood (LLH) values evaluated at multiple D values
    % in parallel in C++ using openMP.
    %
    % This code is described in detail in the paper "Estimation of the Diffusion Constant from Intermittent
    % Trajectories with Variable Position Uncertainties", Relich et. al., 2015
    %
    % In addition this class supplies static methods to implement the 1D LLH calculations for the several different
    % algorithms as described in the paper.  The interface for these calculations is through the computeLLH static
    % method.
    %
    % Finally there are static methods to simulate a 1D observational trace with a given diffusion constant, and
    % to graphically compare the LLH calculations with the various algorithms both for numerical accuracy and for speed.
    %
    % Usage:
    %   dest_obj = DEstimator(Obs, T, V, exposureT)
    %
    % Where Inputs are:
    %  Obs - size=[N,Ndim] - Observed locations in each of the Ndim dimensions, where N is the number of localizations.
    %  T - size=N  - time vector (must be increasing)
    %  SE - size=[N,Ndim] - standard error of observed localizations in each dimension
    %  exposureT - size=1 or N - Exposure Time scalar or vector [default min(dT)]
    %
    % *Note: Floating Point Types - All inputs should be of type double.  The underlying C++ class is templated to produce versions
    % that work with single or double precision floating point numbers, but we have found especially for small 
    % diffusion constants that the computations are much more accurate with double precision.
    %
    % *Note: Units - This implementation maintains a generic unit approach.  The observed positions and observation times can use any
    % units and the diffusion constant will then have units: (dist_unit)^2/(time_unit).  For this reason we do not explicitly mention
    % units anywhere else in the documentation and allow the users to use whatever units are the most natural.  Note that with this in mind
    % the exposureT variable should use the same time units as T, and V should have units of (dist_unit)^2.
    %
    properties (Constant=true)

        % min_variance -- The minimum variance to use for either the observational or diffusion variances.
        %   This prevents numerical issues.  We default to machine epsilon.
        min_variance=eps;

        % min_laplace_variance -- The minimum variance to use for either the observational or
        %   diffusion variances for the Laplace methods.  This prevents numerical inaccuracies in the
        %   Laplace tri-diagonal computations.
        min_laplace_variance=1e-8;

        % Algorithms - The valid algorithm names that the static methods can use to compute the LLH.
        Algorithms={'recursiveC', ... % The most efficient recursive solution in C++, but not parallel.
                    'recursive', ...  % A Matlab implementation of the most efficient recursive solution.
                    'LaplaceC', ... % Laplace's method solution in C++ using tri-diagonal symmetric matrix algorithms.
                    'Laplace', ...  % Laplace's method solution in Matlab using tri-diagonal symmetric matrix algorithms.
                    'LaplaceSparse', ... % Laplace's method solution in Matlab using sparse matrices.
                    'LaplaceDirect', ... % A Matlab implementation of the Laplace solution using a direct approach to the motion blur solution.
                    'MarkovC', ... % Markov's method using tri-diagonal matrix algorithms.
                    'Markov',... % Markov's method  in Matlab using tri-diagonal matrix algorithms.
                    'MarkovSparse',... % Markov's method using sparse matrices.
                    'DSTApproach'... % Original DSTApproach method rewritten to be as faithful to original DSTApproach as possible, but correcting
                    };                    % for constant LLH offset, and our variance terms which have a minimum_variance to prevent numeric difficulties.
    end
    properties
        Ntraj; %Number of trajectories
        Ndim; %Number of dimensions
        exposureT;    %Time over which each observation's data is captured.  Normally the camera exposure time.  This affects motion-blur.
        algorithm;    % Algorithm used for Likelihood calculations (default = recursiveC)
        N;    %cell array [Ntrajx1]: Number of observations
        Obs;  %cell array [Ntrajx1]: Observed positions size=[N(i), Ndim]
        T;    %cell array [Ntrajx1]: Observation times size=N(i)
        SE;   %cell array [Ntrajx1]: Observed position measurement standard errors size=[N, Ndim]
    end
    
    properties (Access = protected, Transient=true)
        initialized=false; %True once the object is correctly initialized by the initializeTrack method
    end
    
    methods
        function obj=DEstimator(varargin)
            % Construct a new DEstimator object to work with a particular track.
            %
            % Usage:
            %   obj=DEstimator();
            % [or]
            %   obj=DEstimator(Obs,T,SE,exposureT);
            %
            % Inputs:
            % [either]
            %   empty - Constructs an empty DEstimator object which will need to have initializeTrack called
            %      later on this object before using it.
            % [or] 
            %  Calls initialize track which accepts:
            %   Obs - size=[N,Ndim] - Observed locations in each of the Ndim dimensions, where N is the number of localizations.
            %   T - size=N  - time vector (must be increasing)
            %   SE - size=[N,Ndim] - standard error of observed localizations in each dimension
            %   exposureT [optional] - scalar double - Time over which each observation's data is captured.  Normally the camera exposure time. [default: min(diif(T))]
            obj=obj@IfaceMixin(@DEstimator_Iface);
            if ~isempty(varargin)
                obj.initializeTrack(varargin{:});
            end
        end
        
        function success=initializeTrack(obj, Obs, T, SE, exposureT, algorithm)
            % Initialize this object with a track and create a new c++ Iface object.
            %
            % This can be called on a currently initialized DEstimator object and it will cause the old
            % object's data to be overwritten, the old C++ object to be deleted and a new C++ to be created.
            %
            % Usage:
            %  success=obj.initializeTrack(Obs, T, SE, exposureT);
            %
            % Inputs: 
            %   Obs - size=[N,Ndim] - Observed locations in each of the Ndim dimensions, where N is the number of localizations.
            %   T - size=N  - time vector (must be increasing)
            %   SE - size=[N,Ndim] - standard error of observed localizations in each dimension
            %   exposureT  [optional]- scalar double - Time over which each observation's data is captured.  Normally the camera exposure time. [default: min(diif(T))]
            %   algorithm  [optional] - name of algorithm to use [Default: recursiveC]
            %
            % Output:
            %   success - True if there were no initialization errors.
            if obj.initialized
                obj.closeIface()
            end
            obj.N = size(Obs,1);
            obj.Ndim = size(Obs,2);
            obj.Obs = double(Obs);
            obj.T = double(T(:));
            assert(length(obj.T)==obj.N);
            assert(all((T(2:end)-T(1:end-1))>0));
            assert(length(SE)==obj.N);
            if size(SE,2)==1 && obj.Ndim>1 % if SE is 1D assume all dimensions have same SE
                obj.SE = repmat(double(SE),1,obj.Ndim);
            else
                obj.SE = double(SE);
            end
            if nargin==4
                %Default exposureT is normally min(dT) in the dataset
                %This is a sane default value.
                exposureT = min(diff(T));
            end
            if nargin < 6
                algorithm = 'recursiveC';
            end
            obj.algorithm = algorithm;
            obj.exposureT = double(exposureT);
            obj.initialized = obj.openIface(obj.Obs, obj.T, obj.SE, obj.exposureT);
            success = obj.initialized;
        end
        
        function llh=LLH(obj,D)
            % Evaluate the log-likelihood (llh) using the most efficient and accurate method (recursiveParallel)
            %
            % Usage:
            %  llh=obj.LLH(D);
            %
            % This works over all Ndim dimensions and computes D LLH values in parallel in C with openMP.  This is
            % the fastest way to evaluate multiple D values LLH's simultaneously.
            %
            % Inputs:
            %   D - scalar or vector of diffusion values to evaluate LLH at.
            % Output:
            %   llh - scalar or vector of log-likelihood calculated at each given D.
            if ~obj.initialized
                error('DEstimator:NotInitialized','Object not initialized');
            end
            
            if strcmp(obj.algorithm,'recursiveC')
                llh=obj.call('LLH',double(D(:)));
            else
                llh = zeros(length(D),1);
                for ii = 1:size(obj.Obs,2)
                    tObs = obj.Obs(:,ii);
                    tSE = obj.SE(:,ii);
                    llh=llh+obj.computeLLH(D, tObs, obj.T, tSE, obj.exposureT, obj.algorithm) ;
                end
            end
        end
        
        function inform=fisher(obj,D)
            % Evaluate the fisher information (inform) using a recursive
            % approach
            %
            % Usage:
            %  llh=obj.LLH(D);
            %
            % This works over all Ndim dimensions and computes D LLH values in parallel in C with openMP.  This is
            % the fastest way to evaluate multiple D values LLH's simultaneously.
            %
            % Inputs:
            %   D - scalar or vector of diffusion values to evaluate LLH at.
            % Output:
            %   llh - scalar or vector of log-likelihood calculated at each given D.
            if ~obj.initialized
                error('DEstimator:NotInitialized','Object not initialized');
            end
            % we currently only calculate the fisher information in the
            % recursive method
            inform = DEstimator.recursiveFisher1D(D, obj.T, obj.SE, obj.exposureT);

        end
        
        function llh=LLHdim(obj,D,dim)
            % Evaluate the log-likelihood (llh) using the most efficient and accurate method (recursiveC),
            % considering only a single dimension dim.
            %
            % Usage: 
            %  llh=obj.LLHdim(D, dim);
            %
            % This works over a single dimension dim, and computes D LLH values in parallel in C with openMP.
            %
            % Inputs:
            %   D - scalar or vector of diffusion values to evaluate LLH at.
            %   dim - dimension to evaluate must be >1 and <=Ndim
            % Output:
            %   llh - scalar or vector of log-likelihood calculated at each given D.
            if ~obj.initialized
                error('DEstimator:NotInitialized','Object not initialized');
            end
            if ~isscalar(dim) || dim<=0 || dim > obj.Ndim
                error('DEstimator:ValueError','dim is not valid');
            end
            dim=dim-1;            
            if strcmp(obj.algorithm,'recursiveC')
                llh=obj.call('LLHdim',double(D(:)),int32(dim));
            else
                ii = dim+1;
                tObs = obj.Obs(:,ii);
                tSE = obj.SE(:,ii);
                llh=obj.computeLLH(D, tObs, obj.T, tSE, obj.exposureT, obj.algorithm);
            end
        end

        function [MLE,MLE_LLH]=MLE(obj)
            % Find the Maximum likelihood estimate using the most efficient and accurate LLH computational
            % method (recursiveC), along with Brent's method for 1D maximization without derivative information.
            %
            % Usage:
            %   [max_D, max_llh, ncalls]=obj.MLE();
            %
            % This works over all Ndim dimensions.
            %
            % Inputs:
            %    None
            % Output:
            %   Dmle - scalar: The maximum likelihood estimate of D
            %   Dmle_llh - The log-likelihood at max_D
            minD = 1e-8;
            maxD = 1e8;
            maxEvals = 1000;
            opts = optimset('MaxFunEvals', maxEvals, 'MaxIter', maxEvals);
            [MLE, MLE_LLH] = fminbnd(@(d) -obj.LLH(d),minD, maxD, opts); 
            MLE_LLH = -MLE_LLH;
        end

        function [posMLE,posSE,llh]=MLEPositions(obj,D)
            % Find the Maximum likelihood estimate of the particles true positions given the track data and
            % a particular D value. This uses the Laplace method to calculate the LLH for each D.  A side effect
            % of this computation is the computation of the maximum likelihood positions of the particles true
            % locations given the Obs, T, SE, exposureT for this object.
            %
            % Usage:
            %   [mle_pos,llh]=obj.MLEPositions(D)
            %
            % This works over all Ndim dimensions.
            %
            % Inputs:
            %    D - size=ND: [Vector or scalar] The number of D values to calculate MLE positions for.
            % Output:
            %   posMLE - size=[N,Ndim,ND]: The maximum likelihood position estimate of the particles true location
            %         for each D value queried.  For each D(j) the corresponding MLE positions are
            %         posMLE(:,:,j), where the rows correspond to observations and the columns to spatial
            %         dimensions.
            %   llh - length=ND: The log-likelihood at each D value queried.          
            ND=length(D);
            posMLE=zeros(obj.N, obj.Ndim, ND);
            posSE=zeros(obj.N, obj.Ndim, ND);
            llh=zeros(ND,1);
            for dim=1:obj.Ndim
                [posMLE_dim, posSE_dim, llh_dim] = obj.posMLE_laplace1D(double(D(:)), obj.Obs, obj.T, obj.SE, obj.exposureT);
                llh = llh+llh_dim;
                posMLE(:,dim,:) = posMLE_dim;
                posSE(:,dim,:) = posSE_dim;
            end
        end
        
        function [MAP,MAP_LPP]=MAP(obj)
            % Find the Maximum A Posteriori estimate using the most efficient and accurate LLH computational
            % method (recursiveC), along with Brent's method for 1D maximization without derivative information.
            %
            % Usage:
            %   [max_D, max_llh, ncalls]=obj.MLE();
            %
            % This works over all Ndim dimensions.
            %
            % Inputs:
            %    None
            % Output:
            %   MAP - scalar: The maximum a posteriori estimate of D
            %   MAP_LPP - The log-proportional posterior at the MAP
            minD = 1e-8;
            maxD = 1e8;
            maxEvals = 1000;
            opts = optimset('MaxFunEvals', maxEvals, 'MaxIter', maxEvals);
            [MAP, MAP_LPP] = fminbnd(@(d) -obj.LLH(d)-0.5*log(obj.fisher(d)),minD, maxD, opts);
            MAP_LPP = -MAP_LPP;
        end
    
    end

    methods (Access=public, Static=true)
        function LLH = computeEnsembleLLH(Ds, Ts, exposureT,  algorithm)
            % This is an extra function to generate the log likelihood of a
            % single D value given ensembles of trajectories presented  in a cell array format.
            %
            % Usage: llH = DEstimator.computeEnsembleLLH(D, Obs, T, SE,
            % exposure T, algorithm);
            %
            %   INPUTS:
            %       Ds - Vector: sampled diffusion constants
            %       Ts - Cell array (Mx1) cell array of trajectory
            %            columns [t x y SEx SEy]
            %       exposureT - the exposure time.  REQUIRED.
            %       algorithm - [optional] the algorithm used to generate the
            %       likelihood. DEstimator.Algorithms for the list, the
            %       (default 'recursiveC')
            if nargin < 4
                algorithm = 'recursiveC';
            end
            if nargin < 3
                error('A scalar exposure time is required as the third input');
            end
            if ~iscell(Ts)
                error('Trajectory input must be a cell array with Nx3 elements in each cell')
            end
            LLH = zeros(length(Ds),1);
            % loop over each trajectory to get the log likelihood value
            Ndim = (size(Ts{1},2)-1)/2;
            for ii = 1:length(Ts)
                T = Ts{ii}(:,1);
                for n=1:Ndim
                    Obs = Ts{ii}(:,1+n);
                    SE = Ts{ii}(:,1+n+Ndim);
                    llh_temp = DEstimator.computeLLH(Ds, Obs, T, SE, exposureT, algorithm);
                    LLH = LLH + llh_temp;
                end
            end            
        end
        
        function [MLE, MLE_LLH, confInt] = computeEnsembleMLE(Ts, exposureT, conf_alpha, algorithm)
            % This function computes the maximum likelihood estimate of an
            % ensemble of trajectories stored in a cell array.
            %
            % We use a LLH ratio test to return approximate confidence intervals assuming a chi^2
            % distribution.  100*(1-conf_alpha)% intervals are returned
            %
            % Usage: llH = DEstimator.computeEnsembleLLH(D, Obs, T, SE,
            % exposure T, algorithm);
            %
            %   INPUTS:
            %       Ts - Cell array (Mx1) cell array of trajectory
            %            columns [t x y SEx SEy]
            %       exposureT - the exposure time.  REQUIRED.
            %       conf_alpha - [optional] [default=0.05] the 100*(1-alpha)% confidence itervals to return.
            %                                              this should be 0<conf_alpha<1
            %       algorithm - [optional] the algorithm used to generate the
            %       likelihood. DEstimator.Algorithms for the list, the
            %       (default 'recursiveC')
            if nargin < 4
                algorithm = 'recursiveC';
            end
            if nargin <3
                if nargout >=3
                    conf_alpha=0.05;
                else
                    conf_alpha=[];
                end
            end
                
            if nargin <2 || isempty(exposureT)
                exposureT = min(cellfun(@(T) min(diff(T(:,1))),Ts)); % minimum time gap in any track
            end
            tempTraj = Ts{1}; % Looking at the first trajectory to get a D upper bound
            maxD = 100*sum(abs(diff(tempTraj(:,2)))) / (tempTraj(2,1)-tempTraj(1,1)); %Maximum possible D value.
            [MLE, MLE_LLH] = fminbnd(@(d) -DEstimator.computeEnsembleLLH(d,Ts,exposureT,algorithm),0, maxD,optimset('MaxFunEvals', 1000, 'MaxIter', 1000)); 
            MLE_LLH = -MLE_LLH; %Change from minimization to maximization

            %Find confidence itervals
            delta = .5*chi2inv(1-conf_alpha,1); %looking for theta solving: 0 = -llh(theta) + llh(theta_mle) - delta 
            CI_func = @(d) -DEstimator.computeEnsembleLLH(d,Ts,exposureT,algorithm) + MLE_LLH - delta;
            theta_lb = fzero(CI_func,[1e-6,MLE]);
            theta_ub = fzero(CI_func,[MLE,1e6]);
            confInt = [theta_lb, theta_ub];
        end
        
        function [Obs, T, SE]=simulate1D(D, N, SE, dT, exposureT, var_gen_method)
            % Simulate a 1D diffusing particle observation data considering the motion-blur as
            % parameterized by exposureT.
            %
            % Usage:
            %   [Obs, T, SE]=DEstimator.simulate1D(D, N, SE, dT, eposureT);
            %
            %   INPUTS:
            %       D  - scalar: diffusion constant
            %       N  - Number of observations to simulate
            %       SE  - Vector or scalar of standard errors
            %       dT - scalar: The time step
            %       exposureT [optional] - The exposure time.  [Default = dT].
            %   OUTPUT:
            %       Obs - vector length=N: simulated observed localizations
            %       T   - vector length=N: Times
            %       SE   - vector length=N: standard error of observed localizations            
            if nargin<5
                exposureT = dT;
            end
            if nargin<6
                var_gen_method = 'gamma';
            end         
            X = DEstimator.simulate1DDiffusion(N+1, D, dT); %make N+1 X values
            [Obs, T, SE] = DEstimator.simulate1DObservationVar(X, SE, D, dT, exposureT,var_gen_method);
        end        
        
        function X=simulate1DDiffusion(N, D, dT)
            %
            % Simulate a 1D diffusion with no observational error
            % 
            %   INPUTS:
            %       D  - scalar: diffusion constant
            %       N  - Number of observations to simulate
            %       dT - scalar: time step
            %   OUTPUT:
            %       X - vector length=N: simulated observed localizations
            X = cumsum(randn(N,1).*sqrt(2*D*dT));
        end
        
        function [Obs, T, SE]=simulate1DObservationVar(X, SE, D, dT, exposureT, var_gen_method)
            % Add observational variance to a 1D diffusion X
            %   INPUTS:
            %       X - vector length=N+1: simulated observed localizations
            %       SE  - Vector or scalar of standard errors
            %       D  - scalar: diffusion constant
            %       N  - Number of observations to simulate
            %       dT - scalar: The time step
            %       exposureT [optional] - The exposure time.  [Default = dT].
            %       var_gen_method [optional] - String: the method name for
            %       how to generate the variable standard errors
            %                                   [default: gamma]
            %   OUTPUT:
            %       Obs- vector length=N: simulated observed localizations
            %       T  - vector length=N: Times
            %       SE - vector length=N: standard error of observed localizations   
            if nargin<6
                var_gen_method='gamma';
            end
            N = numel(X)-1;
            if isscalar(SE)
                %If a scalar V is given 
                SE = DEstimator.generateObsStandardError(SE,N,var_gen_method);
            else
                SE = SE(1:N);
                SE = SE(:);
            end
            T = (0:N-1)'.*dT;
            alpha = exposureT/(2*dT);
            Xbar_mean = (1-alpha)*X(1:end-1) + alpha*X(2:end);
            Xbar_var = 2*D*exposureT*(1/3-alpha/2);
            Obs = Xbar_mean+randn(N,1).*sqrt(SE.^2 + Xbar_var);
        end
        
        function SE=generateObsStandardError(meanSE, N, method)
            if nargin==2
                method = 'gamma';
            end
            switch method
                case 'gamma'
                    %model the distribution of sqrt(V) as a gamma distribution
                    % with given sqrt(V) as mean and shape parameter=4.  This gives a reasonable distribution
                    % of V values, although not really based on any theory.
                    gamscale = 4;
                    SE = gamrnd(gamscale,meanSE/gamscale,N,1);
                case 'constant'
                    % model V as constant.  This is the assumption of the DSTApproach algorithm
                    SE = ones(N,1)*meanSE;
                case 'uniform'
                    %model sqrt(V) as uniform on [0.5*meanSE, 1.5*meanSE)]
                    SE = ((rand(N,1)+0.5)*meanSE);
                otherwise
                    error('DEstimator:generateObsStandardError','Unknown standard error generation method: "%s"',method);
            end
        end
        
        function [MLE,MLE_LLH]=computeMLE(Obs, T, SE, exposureT, algorithm)
            %
            % This is the main entry point for static LLH computations  MLE estimation 
            % This method works with any of the algorithm names in DEstimator.Algorithms.
            %
            % Usage:
            %   llh=DEstimator.computeMLE(D,Obs, T, SE, exposureT, algorithm);
            %
            % Works with abs(D) so LLH(D)=LLH(-D).  This is nice for
            % optimization where the peak is at 0.
            % 
            %   INPUTS:
            %       Obs - size=N: observed points
            %       T - size=N: time vector (must be incremental)
            %       SE - size=N: observed localization errors
            %       exposureT - size=N: Exposure Time vector or scalar (default 1/6 dt)
            %       algorithm - algorithms name from DEstimator.Algorithms
            %                   (default recursiveC)
            %   OUTPUT:
            %       LLH - size=ND: log likelihood of observing Obs with standard error SE
            %       given true diffusion constant is D.
            %
            % Internally this ensures that all the inputs have the correct types and shapes and are
            % passed in a uniform way to each of the computational methods which are implemented as separate
            % C or Matlab static methods.
            %
            if nargin<5
                algorithm = DEstimator.Algorithms{1};
            end
            maxD = 100*sum((diff(Obs)).^2) / (T(2)-T(1)); %Maximum possible D value.
            [MLE, MLE_LLH] = fminbnd(@(d) -DEstimator.computeLLH(d, Obs, T, SE, exposureT,algorithm),0, maxD,optimset('MaxFunEvals', 1000, 'MaxIter', 1000)); 
            MLE_LLH = -MLE_LLH; %Change from minimization to maximization
        end
        
        function LLH=computeLLH(D, Obs, T, SE, exposureT,  algorithm)
            %
            % This is the main entry point for static LLH computations and
            % comparison.  This method works with any of the algorithm
            % names in DEstimator.Algorithms.
            %
            % Usage:
            %   llh=DEstimator.computeLLH(D,Obs, T, SE, exposureT, algorithm);
            %
            % Works with abs(D) so LLH(D)=LLH(-D).  This is nice for
            % optimization where the peak is at 0.
            % 
            %   INPUTS:
            %       D -   vector length=ND: sample diffusion values vector or scalar
            %       Obs - vector length=N: observed points
            %       T -   vector length=N: time vector (must be incremental)
            %       SE -   vector length=N: observed localization standard error (sigma^2)
            %       exposureT - [optional] scalar: Exposure Time [default: min(diff(T))]
            %       algorithm - [optional] algorithms name from DEstimator.Algorithms
            %                   [default: 'recursiveC']
            %   OUTPUT:
            %       LLH - size=ND: log likelihood of observing Obs with standard error SE
            %                      given true diffusion constant is D.
            %
            % Internally this ensures that all the inputs have the correct types and shapes and are
            % passed in a uniform way to each of the computational methods which are implemented as separate
            % C or Matlab static methods.
            %
            if nargin==5
                algorithm = DEstimator.Algorithms{1};
            end
            if nargin==4
               exposureT = min(T(2:end)-T(1:end-1));
            end
            %Check shapes and sizes of input variables
            N = length(Obs);
            if length(Obs)~=numel(Obs)
                error('DEstimator:ValueError','Obs is not 1-D');
            end
            if length(SE)~=numel(SE)
                error('DEstimator:ValueError','SE is not 1-D');
            end
            if length(T)~=N
                error('DEstimator:ValueError','T is incorrect size');
            end
            if length(SE)~=N
                error('DEstimator:ValueError','SE is incorrect size');
            end
            if length(Obs)~=N
                error('DEstimator:ValueError','Obs is incorrect size');
            end
            if ~isscalar(exposureT) && ~length(exposureT)==N
                error('DEstimator:ValueError','exposureT is incorrect size');
            end
            %Ensure input variables are correct orientation (column vector) and type (double)
            Obs = double(Obs(:));
            T = double(T(:));
            SE = double(SE(:));
            D = double(D(:));
            exposureT = double(exposureT);

            if strcmpi(algorithm,'recursiveC')
                LLH = DEstimator.callstatic(@DEstimator_Iface,'LLH_recursive1D', D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'recursive')
                LLH = DEstimator.LLH_recursive1D(D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'LaplaceC')
                LLH = DEstimator.callstatic(@DEstimator_Iface,'LLH_laplace1D', D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'Laplace')
                LLH = DEstimator.LLH_laplace1D(D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'LaplaceSparse')
                LLH = DEstimator.LLH_laplaceSparse1D(D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'LaplaceDirect')
                LLH = DEstimator.LLH_laplaceDirect1D(D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'MarkovC')
                LLH = DEstimator.callstatic(@DEstimator_Iface,'LLH_markov1D', D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'DSTApproach')
                LLH = DEstimator.LLH_DSTApproach1D( D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'Markov')
                LLH = DEstimator.LLH_markov1D(D, Obs, T, SE, exposureT);
            elseif strcmpi(algorithm,'MarkovSparse')
                LLH = DEstimator.LLH_markovSparse1D(D, Obs, T, SE, exposureT);
            else
                error('DEstimator:computeLLH',['Unknown algorithm: ' algorithm]);
            end
        end
        
        function results=compareLLHspeed(Nsizes, max_size, trials, NDvals, flag)
            %
            % Test and record the walltime for execution of each of the
            % algorithms ins DEstimator.Algorithms on many different sized 
            % data sets.
            %
            % Usage:
            %   results=DEstimator.compareLLHspeed(Nsizes, max_size, trials, NDvals);
            %
            % Inputs:
            %   Nsizes -  number of different sizes N to test speed on (suggest ~50)
            %   max_size -  maximum size to test speed on suggest <10000
            %   trials - number of repetitions for each measurements (suggest 2 or 3). Default 3.
            %   NDvals - the number of different D's to test on (suggest 1000 - 10000). Default 1000.
            %   flag - A parameter to test the 4 algorithms (0) or just
            %   recursive vs. DST (1).  Default 0
            % Output:
            %   results - A struct with all the info on the measurements.  Send this struct to plotLLHspeed().
            if nargin<4
                NDvals=1000;
            end
            if nargin==2
                trials=3;
            end
            if nargin < 5
                flag = 0;
            end
            min_size=4;
            
            % flag added to make simple plot for computational tools paper
            if flag
                algs={'recursiveC','DSTApproach'};
            else
                algs={'recursiveC','DSTApproach','LaplaceC','MarkovC'};
            end
            simD=1.0;
            simSE=1.0;
            simdT=1.0;
            exposureT=simdT;
            results.sizes=min(max_size,round(logspace(log10(min_size),log10(max_size),Nsizes)));
            results.trials=trials;
            results.D=linspace(0,5,NDvals)';
            [results.Obs, results.T, results.SE]=DEstimator.simulate1D(simD,max_size,simSE,simdT);
            Nalg=length(algs);
            H=waitbar(0);
            for n=1:Nalg
                alg=algs{n};
                waitbar(n/(Nalg+2),H,alg);
                fun=(@(sz) DEstimator.computeLLH(results.D, results.Obs(1:sz),...
                                    results.T(1:sz),results.SE(1:sz),exposureT,alg));
                results.algorithms.(alg)=timeit(trials, fun, results.sizes);
            end
            
            waitbar(1,H);
            delete(H);
        end
        
        function plotLLHspeed(results)
            % Plot the speed of each LLH method on a log-log plot.
            %
            % Usage:
            %   DEstimator.plotLLHspeed(results);
            %
            % Inputs:
            %   results - A struct with all the info on the timing measurements from DEstimator.compareLLHspeed.
            figure('position',[100,100,600,350]);
            hold('on');
            markers={'p','o','s','d','^','<','>','v','h','x','*','.','+'};
            colors=[0,0,0; 1,0,0; 0,1,0; 0,0,1; 1.0,0.5,0.2; 1,0,1; 0,1,1; 0.5, 0.5, 0.5; 0.5, 1, 0; 0.7, 0.5,0.2];
            algs=fieldnames(results.algorithms);
            Nalgs=length(algs);
            for n=1:Nalgs
                alg=algs{n};
                plot(results.sizes,results.algorithms.(alg).mean,'markerfacecolor',colors(n,:),...
                    'markersize',6,'markeredgecolor','k','Color',colors(n,:),'Marker',markers{n}, 'DisplayName',alg);
            end
            % find the minimum valued result to scale the linear line to
            tempMeanResults = table2array(results.algorithms.(algs{1})(:,1));
            for n = 2:Nalgs
                otherMeanResults = table2array(results.algorithms.(algs{n})(:,1));
                tempMeanResults = min([tempMeanResults otherMeanResults],[],2);
            end
            % find minimum possible slope given all minimum data points
            minslopecands = tempMeanResults./results.sizes';
            minslope = min(minslopecands);
            
            ys=minslope/2*results.sizes;
            plot(results.sizes,ys,'k--','DisplayName','(linear)');
            hold('off');
            ylim([1E-4,3E1]);
            xlim([4, results.sizes(end)]);
            set(gca(),'XScale','log','YScale','log');
            set(gca(),'XTick',[1e2,1e3,1e4,1e5]);
            set(gca(),'FontSize',12);
            grid('on');
            lh=legend('show', 'location','bestoutside');
            set(lh,'interpreter','latex','FontSize',12);
            xlabel('$N$','interpreter', 'latex');
            ylabel('Wall time (s)','interpreter', 'latex')
            title('Execution time','interpreter', 'latex');
        end

        function compareLLHcomputations(ND, N, D, SE, dT, exposureT)
            % Compare the numerical accuracy of each LLH method on linear and logarithmic scales.
            %
            % Usage:
            %   DEstimator.compareLLHcomputations(ND, N, D, SE, dT);
            %
            % This simulates a particle and computes D LLH with each of the algorithms in 
            % DEstimator.Algorithms.
            %
            %In:
            %   ND - scalar the number of D's to test on (suggest 100-1000)
            %   N - The size of the data sets, (#number of localizations)
            %   D - scalar Diffusion constant for the particle we simulate
            %   SE - the standard error scalar or vector for the particle we simulate
            %   dT - the time step for the simulated particle. [default 1.0]
            if nargin<5
                dT=1;
            end
            if nargin<6
                exposureT=dT;  
            end
            
            [Obs, T, SE]=DEstimator.simulate1D(D, N, SE, dT); %#ok<*PROP>
            dest=DEstimator(Obs,T,SE,exposureT);
            mle=dest.MLE()
            minD=1e-8;
            maxD=1e8;
            f=figure();
            position=[0.2 0.2 10.6, 8.2];
            set(f,'PaperType','usletter','PaperOrientation','landscape',...
                  'PaperPosition',position, 'Position',position*96);
            ax=subplot(2,2,1);
            ax2=subplot(2,2,2);
            logax=subplot(2,2,3);
            logax2=subplot(2,2,4);
            axs=[ax,ax2,logax,logax2];
            for i=1:4
                hold(axs(i));
            end

            algs={'recursiveC', 'recursive', 'Laplace', 'LaplaceC', 'LaplaceDirect','Markov', 'MarkovC'};
            Nalgs=length(algs);
            colors=lines(Nalgs);
            linD=linspace(minD,maxD,ND);
            logD=logspace(log10(minD),log10(maxD),ND);
            ref_llh=DEstimator.computeLLH(linD, Obs, T, SE, exposureT,'recursive');
            log_ref_llh=DEstimator.computeLLH(logD, Obs, T, SE, exposureT,'recursive');
            leg_lineseries=cell(1,Nalgs);
            for n=1:Nalgs
                alg=algs{n};
                llh=DEstimator.computeLLH(linD, Obs, T, SE, exposureT, alg);
                logllh=DEstimator.computeLLH(logD, Obs, T, SE, exposureT, alg);
                lin_err=sqrt((llh-ref_llh).^2);
                log_err=sqrt((logllh-log_ref_llh).^2);
                
                leg_lineseries{n}=plot(ax, linD, llh,'-','Color', colors(n,:),'DisplayName', alg);             
                plot(ax2, linD, max(eps,lin_err),'-','Color', colors(n,:),'DisplayName', alg);
                plot(logax, logD, logllh,'-','Color', colors(n,:),'DisplayName', alg);
                plot(logax2, logD, max(eps,log_err),'-','Color', colors(n,:),'DisplayName', alg);
            end
            for i=1:4
                hold(axs(i));
                xlabel(axs(i),'$D$','interpreter','latex');
                set(axs(i),'FontSize',12);
                set(axs(i),'YScale','log');
                if i<=2 %linear axes
                    xlim(axs(i),[-0.1,maxD]);
                else %linear axes
                    set(axs(i),'XScale','log');
                end
            end

            %legend
            [~,objh] =legend(ax,[leg_lineseries{:}],'location','Best');
            set(objh,'linewidth',3);
            
            %titles and axis labels
            title(ax,'LLH Linear D scale','interpreter','latex');
            title(ax2,'LLH Linear D scale: absolute errors','interpreter','latex');
            title(logax,'LLH Log D scale' ,'interpreter','latex');
            title(logax2,'LLH Log D scale: absolute errors','interpreter','latex');
            ylabel(ax,'$\mathrm{LLH}(D)$','interpreter','latex');
            ylabel(logax,'$\mathrm{LLH}(D)$','interpreter','latex');
            ylabel(ax2,'$\Delta\mathrm{LLH}(D)$','interpreter','latex');
            ylabel(logax2,'$\Delta\mathrm{LLH}(D)$','interpreter','latex');
            s=sprintf('$N=%i,D=%.3g(\\mu\\mathrm{m}^2/\\mathrm{s}),\\sqrt{V}=%.3f (\\mu\\mathrm{m}),\\delta t=%.3f$',N,D,mean(SE(:)),dT);
            annotation('textbox', [0 0.9 1 0.1],'string',s,'interpreter','latex',...
                       'EdgeColor','none','HorizontalAlignment', 'center',...
                       'FontSize',12);
        end
                
        function results=compareMLEerror(Nsizes, max_size, D, meanSE, dT, nTrials, var_gen_method)
            % Compare the estimation accuracy of the MLE using the the DSTApproach method 
            % to the reference method [recursiveC] for many simulated trajectories.
            % This makes a structure which is plotted with DEstimator.plotMLEerror
            %
            % Usage:
            %   DEstimator.compareMLEerror(Nsizes, max_size, D, meanSE, nTrials);
            %
            %In:
            %   Nsizes - The number of track sizes to measure (suggest 10-100)
            %   max_size - The maximum size to measure estimation accuracy for (suggest 1e3 - 1e5)
            %   meanSE - The mean of the standard error to set.
            %   dT - The time step (s).  This is also used as the exposure time.
            %   Ntrials - the number of trials at each size to average over (suggest 100-1000)
            %   var_gen_method - the distribution to pull the standard
            %   error
            %   values from (default = gamma, available: uniform, constant)
            if nargin < 7
                var_gen_method = 'gamma';
            end
            
            referenceAlg='recursiveC';
            min_size=4;

            results.nTrials=nTrials;
            results.meanSE=meanSE;
            results.D=D;
            results.dT=dT;
            results.exposureT=results.dT;
            results.sizes=min(max_size,round(logspace(log10(min_size),log10(max_size),Nsizes)));
            algorithms={referenceAlg, 'DSTApproach'};
            Nalgs=length(algorithms);
            Loss=zeros(Nsizes, Nalgs);
            
            h=waitbar(0,'Initializing');
            vmle = zeros(Nsizes,nTrials,Nalgs);
            for m=1:Nsizes
                waitbar(m/Nsizes, h, sprintf('Size:%i',results.sizes(m)));
                for t=1:nTrials
                    [Obs, T, obsSE]=DEstimator.simulate1D(results.D, results.sizes(m), results.meanSE, results.dT,results.exposureT,var_gen_method);
                    for n=1:Nalgs
                        alg=algorithms{n};
                        mle=DEstimator.computeMLE(Obs, T, obsSE, results.exposureT, alg);
                        % Squared Log Loss
                        Loss(m, n)=Loss(m, n)+ (log(results.D)-log(mle))^2;
%                         % Entropy Loss 
%                         Loss(m, n)=Loss(m, n)+ (results.D/mle + log(mle) - log(results.D) -1);
%                         % Stein Loss
%                         Loss(m, n)=Loss(m, n)+ (mle/results.D - log(mle) + log(results.D) -1);
                        vmle(m,t,n) = mle;
                    end
                end
            end
            close(h);
            Loss=Loss./nTrials;
            for n=1:Nalgs
                alg=algorithms{n};
                results.Loss.(alg)=Loss(:,n);
                results.mle.(alg) = vmle(:,:,n);
            end
        end
        
        function plotMLEerror(results)
            % This plots the results produced by the compareMLEerror() method.
            % This should only be called with the output of compareMLEerror().
            figure();
            ax=axes();
            algs=fieldnames(results.Loss);
            Nalgs=length(algs);
            colors={'red', 'blue', 'green', 'black', 'cyan', 'magenta'};
            markers={'o','s','>', '<', 'v', '^'};
            hold('on');
            for n=1:Nalgs
                alg=algs{n};
                plot(results.sizes, results.Loss.(alg), 'markerfacecolor',colors{n},'markeredgecolor','none',...
                        'Color',colors{n},'Marker',markers{n}, 'DisplayName',['\textbf{' alg '}']); 
            end
            % find the minimum valued result to scale the linear (fiducial) line to
            tempMeanResults = results.Loss.(algs{1})(:,1);
            for n = 2:Nalgs
                otherMeanResults = results.Loss.(algs{n})(:,1);
                tempMeanResults = min([tempMeanResults otherMeanResults],[],2);
            end
            % find minimum possible slope given all minimum data points
            minslopecands = tempMeanResults.*results.sizes';
            minslope = min(minslopecands)*0.99;
            % generate proportional fidicual line
            fiducial= minslope*results.sizes.^(-1.);
            plot(results.sizes, fiducial,'k--','DisplayName', '$\mathbf{\propto N^{-1}}$');
            hold('off');
            lh=legend('show','location', 'southwest');
            set(lh,'interpreter', 'latex','FontSize',14);
            xlim([results.sizes(1),results.sizes(end)])
            ta = sprintf('$D=%.3f(\\mu\\mathrm{m}^2/\\mathrm{s}),\\quad \\langle \\sqrt{V} \\rangle=%.3f(\\mu\\mathrm{m})$,',results.D,results.meanSE);
            tb = sprintf('$\\quad \\delta t=%.3f(\\mathrm{s}),\\quad N_{\\mathrm{trials}}=%i$', results.dT, results.nTrials);
            title([ta tb], 'interpreter','latex');
            xlabel('\textbf{Trajectory Length} $\mathbf{N}$', 'interpreter','latex');
            ylabel('\textbf{Estimator Risk} $\mathbf{\langle \ell(\widehat{D},D) \rangle}$', 'interpreter','latex');
            ax.XScale='log';
            ax.YScale='log';
            ax.XGrid='on';
            ax.YGrid='on';
            ax.FontSize=12;
            set(ax,'FontWeight','bold')
            % figure out limits for y
            % find minimum and maximum ylimits
            minY = 1;
            maxY = 1;
            for n = 1:Nalgs
                minY = min(minY, min(results.Loss.(algs{n})));
                maxY = max(maxY, max(results.Loss.(algs{n})));
            end
            ylim([minY*0.99, maxY*1.01]);
        end        
    end % Public static methods
    
    % methods for interval analysis, will integrate better later
    methods (Access=public, Static=true)
        
        function Information = DSTFisher1D(D, T, SE, exposureT)
            N = length(T);
            ND = length(D);
            D = abs(D);
            Information = zeros(ND,1);
            dT = mean(diff(T)); %dT is the scalar used for DSTApproach dT parameter
            rV = mean(SE.*SE); %rV is the mean variance used as the DSTApproach sig2 scalar parameter
            R = exposureT/6/dT; % We fix R to 1/6 * Exposure Time / Avg. Frame Time
            for nn=1:ND
                x = rV/D(nn)/dT-2*R;
                for kk = 1:N-1
                    psi = 2*D(nn)*dT*(1+x*(1-cos(pi*kk/N)));
                    dpsi = 2*dT*(1-2*R*(1-cos(pi*kk/N)));
                    Information(nn) = Information(nn) + (dpsi/psi)^2;
                end
                Information(nn) = 0.5*Information(nn);
            end            
        end
        
        function Information = recursiveFisher1D(D, T, SE, exposureT)
            % Recursive solution to the Fisher Information, pure Matlab implementation,
            % Returns a vector of Fisher information values for a given D
            N = length(SE);
            ND = length(D);
            dt = diff(T);
            Information = zeros(ND,1);            
            for nn = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(nn), dt, SE, exposureT);
                % initial expectation values
                expdmu2 = 0;                
                % initial derivative values
                dvD = 2*dt;
                dvM = -exposureT/3*ones(length(SE),1);
                dalpha = dvM(1)+dvD(1)+dvM(2);
                % initial recursive values
                alpha = vM(1)+vD(1)+vM(2);
                for k = 2:N-1
                    Information(nn) = Information(nn) + (dalpha^2)/2/alpha^2 ...
                        + expdmu2/alpha;
                    % propagate expectation values, higher values first
                    expdmu2 = 1/alpha*(dalpha*vM(k)/alpha - dvM(k))^2 ...
                        + (vM(k)/alpha)^2*expdmu2;
                    % propagate recursive derivatives, higher derivatives first.
                    dalpha = dvD(k)+dvM(k+1)+dvM(k)-2*vM(k)*dvM(k)/alpha ...
                        +(vM(k)/alpha)^2*dalpha;
                    alpha = vD(k)+vM(k+1)+vM(k)-vM(k)^2/alpha;
                end
                Information(nn) = Information(nn) + (dalpha^2)/2/alpha^2 ...
                    + expdmu2/alpha;
            end            
        end 
        
        function [LLH, ObsInformation] = LLH_ObsI_recursive1D(D, Obs, T, SE, exposureT)
            % Recursive solution, pure Matlab implementation, returns information
            % vector as well
            N = length(Obs);
            ND = length(D);
            dt = diff(T);
            LLH = zeros(ND,1);
            ObsInformation = zeros(ND,1);           
            for nn = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(nn), dt, SE, exposureT);
                % initial derivative values
                dvD = 2*dt;
                dvM = -exposureT/3*ones(length(SE),1);
                ddmu = 0;
                ddalpha = 0;
                dalpha = dvM(1) + dvD(1) + dvM(2);
                dmu = 0;
                % initial recursive values
                eta = vM(1)+vD(1);
                mu = Obs(1);
                for k = 2:N-1
                    alpha = eta+vM(k);
                    LLH(nn) = LLH(nn) + log(alpha)+ (Obs(k)-mu)^2/alpha;
                    ObsInformation(nn) = ObsInformation(nn) + (ddalpha-dalpha^2/alpha)/alpha/2 ...
                        + (Obs(k)-mu)^2/alpha^2*( (dalpha)^2/alpha - ddalpha/2) + ...
                        (Obs(k)-mu)/alpha*(2*dalpha*dmu/alpha - ddmu) + (dmu)^2/alpha;
                    % propagate recursive derivatives, higher derivatives first.
                    ddmu = ddmu*vM(k)/alpha + 2/alpha*dmu*(dvM(k) - vM(k)*dalpha/alpha)...
                        + (Obs(k)-mu)/alpha^2*(2*dvM(k)*dalpha - 2*vM(k)*dalpha^2/alpha ...
                        + ddalpha*vM(k));
                    ddalpha = (vM(k)/alpha)^2*ddalpha + 4*vM(k)/alpha^2*dvM(k)*dalpha ...
                        - 2*(vM(k)/alpha)^2/alpha*dalpha^2 - 2/alpha*(dvM(k))^2;
                    dmu = dmu*vM(k)/alpha+ (Obs(k)-mu)/alpha*(vM(k)/alpha*dalpha-dvM(k));
                    dalpha = dvD(k) + dvM(k+1) + dvM(k) + vM(k)/alpha*(vM(k)/alpha*dalpha - 2*dvM(k));
                    % update recursive variables
                    mu = (mu*vM(k)+(Obs(k))*eta)/alpha;
                    eta = vD(k)+vM(k)*eta/alpha;
                end
                alpha = eta+vM(N);
                LLH(nn) = LLH(nn) + (N-1)*log(2*pi) + log(alpha) + (Obs(N)-mu)^2/alpha;
                ObsInformation(nn) = ObsInformation(nn) + (ddalpha-dalpha^2/alpha)/alpha/2 ...
                    + (Obs(N)-mu)^2/alpha^2*( (dalpha)^2/alpha - ddalpha/2) + ...
                    (Obs(N)-mu)/alpha*(2*dalpha*dmu/alpha - ddmu) + (dmu)^2/alpha;
            end            
            LLH = -0.5*LLH;  %factored out -0.5 from everything above            
        end
        
    end
        
    methods (Access=public, Static=true)
        function [vM,vD] = computeVariance(D, dt, SE, exposureT, min_variance)
            %
            % Given the parameters of the observed track, compute the effective (motion blur corrected) variances
            % due to measurement (epsilon in paper) and due to diffusion (omega in paper).  Because all LLH
            % computation algorithms implemented here use this common function we can ensure that they all see
            % equivalently normalized and corrected variances.
            %
            % Inputs:
            %   D - length=1: scalar diffusion constant
            %   dt - length=N-1: delta T, time steps (can be non uniform)
            %   SE - length=N-1: delta T, time steps (can be non uniform)
            %   exposureT - length=N-1: delta T, time steps (can be non uniform)
            %   min_variance - length=1: The smallest absolute value tolerated for any value in vD or vM.
            %        if a computed variance is less than min_variance, we set it to min_variance and preserve
            %        the sign.  [Default: DEstimator.min_variance]
            % Outputs:
            %   vM - length=N: Variance due to measurement (motion blur corrected)
            %   vD - length=N-1: Variance due to diffusion
            %
            %
            % Bounding the minimum absolute value on each variance value prevents numerical difficulties.  We
            % recommend using units that keep variances large compared to machine epsilon.
            % All methods use normally call this with the static value DEstimator.min_variance, except the
            % Laplace method which uses DEstimator.min_laplace_variance, which is slightly larger, as smaller values
            % produce numerical errors near D=0.
            %
            if nargin==4
                min_variance=DEstimator.min_variance;
            end
            D = abs(D); %This makes the estimates symmetric around the origin with respect to D, and prevents errors from negative D
            vD = max(min_variance,2*D*dt);
            vM = SE.*SE - 2*D*exposureT/6;
            vM(vM>=0 & vM<min_variance) =  min_variance;
            vM(vM<0  & vM>-min_variance) = -min_variance;
        end
        
        function LLH = LLH_recursive1D(D, Obs, T, SE, exposureT)
            % Recursive solution, pure Matlab implementation
            N = length(Obs);
            ND = length(D);
            dt = diff(T);
            LLH = zeros(ND,1);
            for n = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(n), dt, SE, exposureT);
                eta = vM(1)+vD(1);
                mu = Obs(1);
                for k = 2:N-1
                    alpha = eta+vM(k);
                    LLH(n) = LLH(n) + log(alpha)+ (Obs(k)-mu)^2/alpha;
                    mu = (mu*vM(k)+Obs(k)*eta)/alpha;
                    eta = vD(k)+vM(k)*eta/alpha;
                end
                alpha = eta+vM(N);
                LLH(n) = LLH(n) + (N-1)*log(2*pi) + log(alpha) + (Obs(N)-mu)^2/alpha;
            end
            LLH = -0.5*LLH;  %factored out -0.5 from everything above
        end

        function LLH = LLH_laplaceSparse1D(D, Obs, T, SE, exposureT)
            % Laplace's method using symmetric sparse matrix algorithms
            dt = diff(T);
            N = length(Obs);
            ND = length(D);
            LLH = zeros(ND,1);
            for n = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(n), dt, SE, exposureT, DEstimator.min_laplace_variance/N^2);
                hess = DEstimator.makeSparseLaplaceMatrix(vM,vD);
                theta = hess \ (Obs./vM);
                LLH(n) = -0.5*( (N-1)*log(2*pi) ...
                               +sum((Obs-theta).^2 ./ vM) ...
                               +sum(diff(theta).^2 ./ vD) ...
                               +sum(log(abs(vM))) ...
                               +sum(log(vD)) ...
                               +real(DEstimator.sparse_logdet(hess)));
            end
        end
    
        function LLH = LLH_markov1D(D, Obs, T, SE, exposureT)
            % Markov's solution to a generalized version of DSTApproach's approach.
            % Like the recursive solution but unlike DSTApproach's method this algorithm gives
            % the exact solution for systems with non-constant SE and dT.
            % Uses symmetric tri-diagonal matrix algorithms
            dt = diff(T);
            N = length(Obs);
            ND = length(D);
            LLH = zeros(ND,1);
            for n = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(n), dt, SE, exposureT);
                [hessA,hessB] = DEstimator.makeTriDiagMarkovMatrix(vM,vD);
                dObs = diff(Obs);
                theta = dot(dObs,DEstimator.tridiagonal_solve(hessA, hessB, dObs));
                LLH(n) = -0.5*( (N-1)*log(2*pi) ...
                                +DEstimator.tridiagonal_logdet(hessA,hessB) ...
                                +theta);
            end
        end

        function LLH = LLH_markovSparse1D(D, Obs, T, SE, exposureT)
            % Markov's solution to a generalized version of DSTApproach's approach.
            % Like the recursive solution but unlike DSTApproach's method this algorithm gives
            % the exact solution for systems with non-constant SE and dT.
            % Uses sparse matrix algorithms
            dt = diff(T);
            dObs = diff(Obs);
            N = length(Obs);
            ND = length(D);
            LLH = zeros(ND,1);
            for n = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(n), dt, SE, exposureT);
                hess = DEstimator.makeSparseMarkovMatrix(vM,vD);
                theta = dot(dObs, hess\dObs);
                LLH(n) = -0.5*( (N-1)*log(2*pi) ...
                                +DEstimator.sparse_logdet(hess) ...
                                +theta);
            end
        end

        function LLH = LLH_DSTApproach1D(D, Obs, T, SE, exposureT)
            % DSTApproach original code, for comparison.  
            % This method cannot handle non-constant SE or dt, and uses the mean V value and first dT for
            % all observations.
            % 
            % The core of the code is based on Likelihood_subfunction1D in
            % DSigma2_Est.m supplemental to 
            % Berglund and Michalet [Phys. Rev. E 85, 061916]
            % Note: Minor corrections made for constant offset and to map to real values only.
            N = length(Obs);
            ND = length(D);
            D = abs(D);
            dX = diff(Obs);%-repmat(mean(diff(Obs)),N-1,1);
            dXX = dst(dX).^2;
            LLH = zeros(ND,1);
            dT = mean(diff(T)); %dT is the scalar used for DSTApproach dT parameter
            rV = mean(SE.*SE); %rV is the mean variance used as the DSTApproach sig2 scalar parameter
            R = exposureT/6/dT; % We fix R to 1/6 * Exposure Time / Frame Time
            for n=1:ND
                alpha = 2*D(n)*dT*(1-2*R) + 2*rV;
                beta = 2*R*D(n)*dT - rV;

                eigvec = alpha + 2*beta*cos(pi*(1:N-1)./N);
                eigvec = reshape(eigvec, size(dXX));

                LLH(n) = -.5*sum(log(eigvec) + (2/N).*dXX./eigvec);
            end
            LLH = LLH-(.5*(N-1)*log(2*pi)); %Added this constant correction to original computation
            LLH = real(LLH); %Added to map to real part
        end

        function [LLH,  posMLE]= LLH_laplaceDirect1D(D, Obs, T, SE, exposureT)
            % This is a different version of the Laplace method that is based on a 
            % Eq.(5) in the main text which is the direct description of the motion-blurred
            % probability where each observation is parametrized on both x_i and x_{i+1}.
            % The probability distributions for P(O_i|X_i,X_i+1) is based on the derivations
            % in the supplemental text sections 2.1 and 2.2.
            %
            % This is included mainly to confirm that the computation is equivalent to the
            % transformed Eq. 13 in the main text from which the other Laplace-based methods
            % are derived.  This is probably not the fastest or best algorithm for general use.
            %
            % We use the symmetric tri-diagonal matrix algorithms for the computation
            dt = diff(T);
            dt(end+1) = dt(end); % add an extra time point since we are in the formulation with N+1 X values
            alpha = exposureT./(2*dt);
            N = length(Obs);
            ND = length(D);
            LLH = zeros(ND,1);
            posMLE = zeros(N+1, ND);
            min_var=DEstimator.min_laplace_variance;
            for n = 1:ND
                %Compute the variances vD and vM which are different for this method
                vD = max(min_var,2*abs(D(n))*dt);
                vM = SE.*SE+2*abs(D(n)).*exposureT.*( 0.5*(1-alpha) - 1/6 );
                vM(vM>=0 & vM<min_var) =  min_var;
                vM(vM<0  & vM>-min_var) = -min_var;
                %Formulate the matrix equation
                vMsign = prod(sign(vM));
                [hessA, hessB] = makeTriDiagLaplaceDirectMatrix(vM,vD);
                rhs_solve = zeros(length(Obs)+1,1);
                rhs_solve(1:end-1) = rhs_solve(1:end-1) + (Obs./vM).*(1-alpha);
                rhs_solve(2:end) = rhs_solve(2:end) + (Obs./vM).*(alpha);
                %Solve the matrix equation and compute LLH
                theta = DEstimator.tridiagonal_solve(hessA, hessB, rhs_solve);
                posMLE(:,n) = theta; % save MLE position
                LLH(n) = -0.5*( (N-1)*log(2*pi)...
                    +sum((Obs-theta(1:end-1).*(1-alpha) - theta(2:end).*alpha).^2 ./ vM)...
                    +sum(diff(theta).^2 ./ vD) ...
                    +sum(log(abs(vM)))...
                    +sum(log(vD))...
                    + vMsign*(DEstimator.tridiagonal_logdet(hessA,hessB)));
            end
            function [a,b] = makeTriDiagLaplaceDirectMatrix(vM,vD )
                %inner helper function to make the correct matrix for this Direct Laplace version
                a = -1./vD + alpha.*(1-alpha)./vM;
                b = zeros(length(vD)+1,1);
                r = (1-alpha).^2./vM + 1./vD;
                c = alpha.^2./vM + 1./vD;
                b(1:end-1) = b(1:end-1) + r;
                b(2:end) = b(2:end) + c;
            end
        end      
        
        function [LLH, posMLE] = LLH_laplace1D(D, Obs, T, SE, exposureT)
            % Laplace's method using symmetric tri-diagonal matrix algorithms
            %
            % [out] LLH - Length:ND - The log-likelihood at each given D value
            % [out] posMLE - Size:[N,ND] - The maximum likelihood estimates of true positions X
            %                 at each observation point, given each of the D values.  The j-th
            %                 column corresponds to the MLE positions when D is D(j).
            dt = diff(T);
            N = length(Obs);
            ND = length(D);
            LLH = zeros(ND,1);
            posMLE = zeros(N, ND);
            for n = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(n), dt, SE, exposureT, DEstimator.min_laplace_variance/N^2);
                vMsign = prod(sign(vM));
                [hessA, hessB] = DEstimator.makeTriDiagLaplaceMatrix(vM,vD);
                theta = DEstimator.tridiagonal_solve(hessA, hessB, Obs./vM);
                posMLE(:,n) = theta; %save MLE position
                LLH(n) = -0.5*( (N-1)*log(2*pi)...
                               +sum((Obs-theta).^2 ./ vM) ...
                               +sum(diff(theta).^2 ./ vD) ...
                               +sum(log(abs(vM))) ...
                               +sum(log(vD)) ...
                               + vMsign*(DEstimator.tridiagonal_logdet(hessA,hessB)));
            end
        end
        
        function [posMLE, posSE, LLH] = posMLE_laplace1D(D, Obs, T, SE, exposureT)
            % Laplace's method using symmetric tri-diagonal matrix algorithms
            %
            % [out] posMLE - Size:[N,ND] - The maximum likelihood estimates of true positions X
            %                 at each observation point, given each of the D values.  The j-th
            %                 column corresponds to the MLE positions when D is D(j).
            % [out] posSE - Size:[N,ND] - the standard error of each
            %                 maximum likelihood estimate
            % [out] LLH - Length:ND - The log-likelihood at each given D value
            dt = diff(T);
            N = length(Obs);
            ND = length(D);
            LLH = zeros(ND,1);
            posMLE = zeros(N, ND);
            posSE = zeros(N, ND);
            for n = 1:ND
                [vM,vD] = DEstimator.computeVariance(D(n), dt, SE, exposureT, DEstimator.min_laplace_variance/N^2);
                vMsign = prod(sign(vM));
                [hessA, hessB] = DEstimator.makeTriDiagLaplaceMatrix(vM,vD);
                theta = DEstimator.tridiagonal_solve(hessA, hessB, Obs./vM);
                posMLE(:,n) = theta; %save MLE position
                posCRLB = DEstimator.calculate_invhessdiag(hessA, hessB);
                posSE(:,n) = sqrt(posCRLB); % save position SEs
                LLH(n) = -0.5*( (N-1)*log(2*pi)...
                               +sum((Obs-theta).^2 ./ vM) ...
                               +sum(diff(theta).^2 ./ vD) ...
                               +sum(log(abs(vM))) ...
                               +sum(log(vD)) ...
                               + vMsign*(DEstimator.tridiagonal_logdet(hessA,hessB)));
            end
        end
        
        function posCRLB = calculate_invhessdiag(hessA, hessB)
           % function to calculate the diagonal element of the inverse of
           % the hessian matrix from the Laplace method
           % [out] posCRLB - Inverse of the hessian matrix, or CRLB of the
           % true particle positions given their MLE estimates
           % uses algorithm from Emrah Kilic: Explicit formula for the
           % inverse of a tridiagonal matrix by backward continued
           % fractions, Applied Mathematics and Computation: 2008
           numelem = length(hessB);
           posCRLB = zeros(numelem,1); % our output, the CRLB of the positions
           Bc = zeros(numelem,1); % the backwards fraction
           Pnum = zeros(numelem+2,1); % used to get the backwards fraction
           % put in first two elements of the recursion.
           Pnum(1) = 1;
           Pnum(2) = hessB(1);
           Bc(1) = Pnum(2)/Pnum(1);
           % calculate the backwards fractions.
           for ii = 2:numelem
               Pnum(ii+1) = hessB(ii)*Pnum(ii)-hessA(ii-1)^2*Pnum(ii-1);
               Bc(ii) = Pnum(ii+1)/Pnum(ii);
           end
           % calculate the positional CRLBs from the Hessian
           for ii = 1:numelem
              posCRLB(ii) = 1/Bc(ii);
              for kk = ii+1:numelem
                 tempval = 1/Bc(kk);
                 for tt = ii:kk-1
                     tempval = tempval*hessA(tt)^2/Bc(tt)^2;
                 end
                 posCRLB(ii) = posCRLB(ii) + tempval;
              end
           end
        end
               
        function [a,b]=makeTriDiagLaplaceMatrix(vM,vD)
            % Make the diagonal vectors for the symmetric tri-diagonal matrix for the Laplace method.
            % Inputs:
            %   vM - length=N: Variance due to measurement (motion blur corrected)
            %   vD - length=N-1: Variance due to diffusion
            % Output:
            %   a - length=N-1: +1/-1 diagonal elements (i.e. the diagonals +1 and -1 away from the central diagonal)
            %   b - length=N: Central diagonal elements
            a = -1./vD;
            b = 1./vM;
            b(1:end-1) = b(1:end-1)-a;
            b(2:end) = b(2:end)-a;
        end

        function M=makeSparseLaplaceMatrix(vM,vD)
            % Make the sparse matrix for the Laplace method.
            % Inputs:
            %   vM - length=N: Variance due to measurement (motion blur corrected)
            %   vD - length=N-1: Variance due to diffusion
            % Output:
            %   M - size=[N,N]: The sparse matrix
            N = length(vM);
            [hessA,hessB] = DEstimator.makeTriDiagLaplaceMatrix(vM,vD);
            M = spdiags([[hessA;0], hessB, [0;hessA]],[-1,0,1],N,N);
        end

        function [a,b]=makeTriDiagMarkovMatrix(vM,vD)
            % Make the diagonal vectors for the symmetric tri-diagonal matrix for the Markov method.
            % Inputs:
            %   vM - length=N: Variance due to measurement (motion blur corrected)
            %   vD - length=N-1: Variance due to diffusion
            % Output:
            %   a - length=N-2: +1/-1 diagonal elements (i.e. the diagonals +1 and -1 away from the central diagonal)
            %   b - length=N-1: Central diagonal elements
            a = -vM(2:end-1);
            b = vM(1:end-1)+vM(2:end)+vD;
        end

        function M=makeSparseMarkovMatrix(vM,vD)
            % Make the sparse matrix for the Markov method.
            % Inputs:
            %   vM - length=N: Variance due to measurement (motion blur corrected)
            %   vD - length=N-1: Variance due to diffusion
            % Output:
            %   M - size=[N-1,N-1]: The sparse matrix
            N = length(vM);
            [hessA,hessB] = DEstimator.makeTriDiagMarkovMatrix(vM,vD);
            M = spdiags([[hessA;0], hessB, [0;hessA]],[-1,0,1],N-1,N-1);
        end

        function Ldet=tridiagonal_logdet(a, b)
            % Compute the log determinant of a symmetric tri-diagonal matrix of size=[k, k].
            %
            % This runs in linear time.
            %
            % Inputs:
            %   a - length=k-1: +1/-1 diagonal elements (i.e. the diagonals +1 and -1 away from the central diagonal)
            %   b - length=k: Central diagonal elements
            % Outputs:
            %  Ldet - The log determinant of the matrix
            N = length(b);
            Sdet = sign(b(1));
            Ldet = log(abs(b(1)));
            D = -a(1)/b(1);
            for n = 2:N-1
                B = b(n)+a(n-1)*D;
                D = -a(n)/B;
                Sdet = Sdet*sign(B);
                Ldet = Ldet+log(abs(B));
            end
            B = b(N)+a(N-1)*D;
            Sdet = Sdet*sign(B);
            Ldet = Ldet+log(abs(B));
            Ldet = Ldet*Sdet;
        end

        function x=tridiagonal_solve(a, b, y)
            % Solve the linear system Mx=y, where M is a symmetric tri-diagonal matrix of size=[k, k], given by diagonals a & b.
            %
            % This runs in linear time
            %
            % Inputs:
            %   a - length=k-1: +1/-1 diagonal elements of M (i.e. the diagonals +1 and -1 away from the central diagonal)
            %   b - length=k: Central diagonal elements of M
            %   y - length=k: Right hand side of linear system to solve.
            % Outputs:
            %   x - length=k: The solution to the linear system.
            N = length(b);
            D = zeros(N-1,1);
            E = zeros(N-1,1);
            x = zeros(N,1);
            D(1) = -a(1)/b(1);
            E(1) = y(1)/b(1);
            for n = 2:N-1
                B = b(n)+a(n-1)*D(n-1);
                D(n) = -a(n)/B;
                E(n) = (y(n)-a(n-1)*E(n-1))/B;
            end
            x(N) = (y(N)-a(N-1)*E(N-1))/(b(N)+a(N-1)*D(N-1));
            for n = N-1:-1:1
                x(n) = E(n)+D(n)*x(n+1);
            end
        end

        function Ldet=sparse_logdet(M)
            % In:
            %   M - sparse matrix
            % Out:
            %  Ldet - The log determinant.
            %
            % det is product of diag elements of triangular matrix, logdet follows simply
            % This uses LU factorization.
            Ldet = sum(log(diag(lu(M))));
        end       
    end %private static methods
end %classdef
