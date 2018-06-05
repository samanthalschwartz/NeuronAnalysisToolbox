% Mark J. Olah (mjo@cs.unm.edu)
% 03-2014

classdef MappelBase < IfaceMixin
    %
    %
    % Design decisions
    %  * All computational functions are vectorized.  We take in inputs where the
    %    size of the last array dim represents the number of inputs to process. 
    %    this allows the memory layout to be contiguous for in-place 
    %    parallel thread execution.  Thus paramter vectors are column vectors, because
    %    matlab is column-major oriented.

    properties (Constant=true)
        %MinSize: The minimum imsize of an image in pixels
        MinSize=4;
        EstimationMethods={'Heuristic',...          %[C++(OpenMP)] Hueristic starting point guesstimate
                           'Newton',...             %[C++(OpenMP)] Newtons method with modified cholesky step
                           'NewtonRaphson',...      %[C++(OpenMP)] Newtons method using a diagonal hessian approximation
                           'QuasiNewton',...        %[C++(OpenMP)] Quasi-Newton method
                           'SimulatedAnnealing',... %[C++(OpenMP)] Simulated Annealing global maximization method
                           'CGauss',...             %[C(OpenMP)] C NewtonRaphson implementation from Smith et.al. Nat Methods (2010)
                           'CGaussHeuristic',...    %[C(OpenMP)] C Hueristic starting guestimate implementation from Smith et.al. Nat Methods (2010)
                           'GPUGauss',...           %[C(CUDA)] CUDA implementation from Smith et.al. Nat Methods (2010)
                           'trust-region-reflective',... %[MATLAB (fminunc)] Single threaded trust-region based constrained optimization from Matlab Optimization toolbox
                           'interior-point'...           %[MATLAB (fminunc)] Single threaded interior-point based constrained optimization from Matlab Optimization toolbox
                           };
    end
    
    properties (SetAccess = protected)
        imsize; % 2D:[X Y] or Hyperspectral:[X Y L]
        psf_sigma; % 2D:[X Y] or Hyperspectral:[X Y L]
    end
    
    properties 
        GPUGaussMLE_Iterations=15;
    end

    properties (Abstract=true,Access=protected)
        GPUGaussMLEFitType;
    end
    
    methods
        function obj = MappelBase(iface,imsize,psf_sigma)
            % obj = MappelBase(iface,imsize,psf_sigma) - Make a new MappelBase for
            % point localization in 2D.
            % (in) imsize: int [2x1] - size of image in pixels on each side (min: obj.MinSize)
            % (in) psf_sigma: double [2x1] >0 - size of PSF in pixels
            % (out) obj - A new object
            obj=obj@IfaceMixin(iface);
            if isscalar(imsize)
                imsize=[imsize, imsize];
            end
            if any(imsize<obj.MinSize)
                error('MappelBase:constructor','Invalid imsize')
            end
            if isscalar(psf_sigma)
                psf_sigma=[psf_sigma, psf_sigma];
            end
            if any(psf_sigma<=0)
                error('MappelBase:constructor','Invalid psf_sigma')
            end
            obj.psf_sigma=double(psf_sigma(:)');
            obj.imsize=double(imsize(:)');

            initialized = obj.openIface(int32(obj.imsize), obj.psf_sigma);
            if(~initialized) 
                error('MappelBase:Constructor','Bad initialization');
            end
        end
                
        function setHyperparameters(obj, prior)
            % Set the hyperparamters of the model prior
            % (in) prior - The vecotor of hyperparamters (model specific format)
            prior=double(prior(:));
            if size(prior,2) ~= 1 || size(prior,1) ~= obj.nHyperParams
                error('MappelBase:setHyperparameters', 'Invalid prior size.');
            elseif any(prior<0)
                error('Gauss2DModel:setHyperparameters', 'Priors must be non-negative');
            end
            obj.call('setHyperparameters', prior);
        end
        
        function theta=samplePrior(obj, count)
            % theta=obj.samplePrior(count)
            % (in) count: int (default 1) number of thetas to sample
            % (out) theta: A (nParams X n) double of theta values
            if nargin==1
                count=1;
            end
            count=obj.checkCount(count);
            theta= obj.call('samplePrior', count);
        end
        
        function bounded=boundTheta(obj, theta)
            % (in) theta - A theta to correct to ensure it is in boundsa
            % (out) bounded - A corrected theta that is now inbounds
            if size(theta,1) ~= obj.nParams
                if length(theta) == obj.nParams
                    theta=theta';
                else
                    error('MappelBase:boundTheta', 'Invalid theta shape');
                end
            end
            bounded=obj.call('boundTheta', theta);
        end
        
        function inbounds=thetaInBounds(obj, theta)
            % (in) theta - A theta to check if it is valid
            % (out) iboundsa - True if the theta is valid
            if size(theta,1) ~= obj.nParams
                if length(theta) == obj.nParams
                    theta=theta';
                else
                    error('MappelBase:checkTheta', 'Invalid theta shape');
                end
            end
            inbounds=obj.call('thetaInBounds', theta);
        end

        function image=modelImage(obj, theta)
            % image=obj.modelImage(theta)
            % (in) theta: an (nParams X n) double of theta values
            % (out) image: a double (imsize X imsize X n) image stack
            theta=obj.checkTheta(theta);
            image= obj.call('modelImage', theta);
        end
        
        function image=modelDipImage(obj, theta)
            % image=obj.modelDipImage(theta)
            % This just calls obj.modelImage(theta), then rotates the image so that x
            % goes from left to right, and returns a DIP image (or image
            % sequence.)
            % (in) theta: an (nParams X n) double of theta values
            % (out) image: a (imsize X imsize X n) DIP image stack of model
            %              (mean) value at each pixel.
            image = dip_image(obj.modelImage(theta));
        end

        function image=simulateImage(obj, theta, count)
            % image=obj.simulateImage(theta, count)
            % If theta is size (nParams X 1) then count images with that theta are
            % simulated.  Default count is 1.  If theta is size (nParams X n) with n>1
            % then n images are simulated, each with a seperate theta, and count is ignored.
            % (in) theta: a single (nParams X 1) or (nParams X n) double theta value.
            % (in) count[optional]: the number of independant images to
            %       generate only used if theta is a single vector
            % (out) image: a double size:[imsizeY, imsizeX, n] image stack, all sampled with params theta
            theta=obj.checkTheta(theta);
            if nargin==2
                count=1;
            elseif nargin==3
                if(count>1 && size(theta,2)>1)
                    error('MappelBase:simulateImage', 'Cannot set count>1 when theta_stack has more than one theta');
                end
            end
                    
            count=obj.checkCount(count);
            image= obj.call('simulateImage', theta, count);
        end

        function image=simulateDipImage(obj, varargin)
            % image=obj.simulateDipImage(theta)
            % This just calls obj.simulateImage(theta), then rotates the image so that x
            % goes from left to right, and returns a DIP image (or image
            % sequence.)
            % (in) theta: an (nParams X n) double of theta values
            % (in) count[optional]: the number of independant images to
            %       generate only used if theta is a single vector
            % (out) image: a size:[imsizeY, imsizeX, n] DIP image stack of simulated
            %              data
            image = dip_image(obj.simulateImage(varargin{:}));
        end

        
        function llh=LLH(obj, image, theta)
            % llh=obj.LLH(image, theta)
            % This takes in a N images and M thetas.  If M=N=1, 
            % then we return a single LLH.  If there are N=1 images
            % and M>1 thetas, we return M LLHs of the same image with each of 
            % the thetas.  Otherwise, if there is M=1 thetas and N>1 images,
            % then we return N LLHs for each of the images given theta
            % (in) image: a double (imsize X imsize X N) image stack
            % (in) theta: an (nParams X M) double of theta values
            % (out) llh: a (1 X max(M,N)) double of log_likelihoods
            theta=obj.checkTheta(theta);
            image=obj.checkImage(image);
            if size(image,3) ~= size(theta,2)  && size(theta,2)~=1 && ~ismatrix(image)
                error('MappelBase:LLH','image and theta must match dims');
            end
            llh=obj.call('LLH', image, theta);
        end
        
        function grad=modelGrad(obj, image, theta)
            % grad=obj.modelGrad(obj, image, theta) - Compute the model gradiant.
            % This takes in a N images and M thetas.  If M=N=1, 
            % then we return a single Grad.  If there are N=1 images
            % and M>1 thetas, we return M Grads of the same image with each of 
            % the thetas.  Otherwise, if there is M=1 thetas and N>1 images,
            % then we return N Grads for each of the images given theta
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) theta: an (nParams X M) double of theta values
            % (out) grad: a (nParams X max(M,N)) double of gradiant vectors
            theta=obj.checkTheta(theta);
            image=obj.checkImage(image);
            if size(image,3) ~= size(theta,2)  && size(theta,2)~=1 && ~ismatrix(image)
                error('MappelBase:modelGrad','image and theta must match dims');
            end
            grad=obj.call('modelGrad', image, theta);
        end

        function hess=modelHessian(obj, image, theta)
            % hess=obj.modelHessian(obj, image, theta) - Compute the model hessian
            % This takes in a N images and M thetas.  If M=N=1, 
            % then we return a single Hessian.  If there are N=1 images
            % and M>1 thetas, we return M Hessian of the same image with each of 
            % the thetas.  Otherwise, if there is M=1 thetas and N>1 images,
            % then we return N Hessians for each of the images given theta
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) theta: an (nParams X M) double of theta values
            % (out) hess: a (nParams X nParams X max(M,N)) double of hessian matricies
            theta=obj.checkTheta(theta);
            image=obj.checkImage(image);
            if size(image,3) ~= size(theta,2)  && size(theta,2)~=1 && ~ismatrix(image)
                error('MappelBase:modelHessian','image and theta must match dims');
            end
            hess=obj.call('modelHessian', image, theta);
        end
                
        function hess=modelPositiveHessian(obj, image, theta)
            % hess=obj.modelPositiveHessian(obj, image, theta) - Compute the model hessian
            % This takes in a N images and M thetas.  If M=N=1, 
            % then we return a single Hessian.  If there are N=1 images
            % and M>1 thetas, we return M Hessian of the same image with each of 
            % the thetas.  Otherwise, if there is M=1 thetas and N>1 images,
            % then we return N Hessians for each of the images given theta
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) theta: an (nParams X M) double of theta values
            % (out) hess: a (nParams X nParams X max(M,N)) double of hessian matricies
            theta=obj.checkTheta(theta);
            image=obj.checkImage(image);
            if size(image,3) ~= size(theta,2)  && size(theta,2)~=1 && ~ismatrix(image)
                error('MappelBase:modelPositiveHessian','image and theta must match dims');
            end
            hess=obj.call('modelPositiveHessian', image, theta);
        end

        function crlb=CRLB(obj, theta)
            % crlb=obj.CRLB(obj, theta) - Compute the Cramer-Rao Lower Bound at theta
            % (in) theta: an (nParams X n) double of theta values
            % (out) crlb: a  (nParams X n) double of the cramer rao lower bounds.
            %             these are the lower bounds on the variance at theta 
            %             for any unbiased estimator.
            theta=obj.checkTheta(theta);
            crlb= obj.call('CRLB', theta);
        end

        function e_std=estimationAccuracy(obj, theta)
            % e_std=obj.estimationAccuracy(obj, theta) - Compute the estimation accuracy as the sqrt(CRLB)
            % (in) theta: an (nParams X n) double of theta values
            % (out) e_std: a  (nParams X n) double of the estimation accuracy standard deviation based on the CRLB.
            crlb= obj.call('CRLB', theta);
            e_std=sqrt(crlb);
        end

        function fisherI=fisherInformation(obj, theta)
            % fisherI=obj.FisherInformation(obj, theta) - Compute the Fisher Information matrix
            % at theta
            % (in) theta: an (nParams X n) double of theta values
            % (out) fisherI: a  (nParams X nParms) matrix of the fisher information at theta 
            theta=obj.checkTheta(theta);
            fisherI= obj.call('fisherInformation', theta);
        end
    
        function obsI=observedInformation(obj, im, theta)
            % fisherI=obj.ObservedInformation(obj, theta) - Compute the Fisher Information matrix
            % at theta
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) theta: an (nParams X 1) double of theta values
            % (out) obsI: the observed information matrix (nParams X nParms) for the data sequence of 
            %               images that are all drawn from the distribution with parameters theta 
            im=obj.checkImage(im);
            theta=obj.checkTheta(theta);
            obsI=zeros(obj.nParams);
            for n=1:size(im,ndims(im))
                obsI = obsI - obj.call('modelHessian',im(:,:,n), theta);
            end
        end

        function score = scoreFunction(obj, im, theta)
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) theta: an (nParams X M) double of theta values
            % (out) grad: a (nParams X max(M,N)) double of gradiant vectors
            theta=obj.checkTheta(theta);
            image=obj.checkImage(im);
            grad=obj.call('modelGrad', image, theta);
        end

        function [theta, crlb, llh, stats]=estimate(obj, image, estimator_name, theta_init)
            % [theta, crlb, llh]=obj.estimate(image, name) - estimate theta's
            % crlb's and llh's for each image in stack.  
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) estimator_name: (optional) name for the optimization method. (default = 'Newton')  
            %           Valid names are in obj.EstimationMethods
            % (in) theta_init: (optional) Initial theta guesses size (nParams x n).  Values of 0 indicate
            %            that we have no initial guess and the estimator should form its own guess.
            % (out) theta: a (nParams X n) double of estimated theta values
            % (out) crlb: a (nParams X n) estimate of the CRLB for each parameter estimate.
            %             This gives the approximiate variance in the theta parameters
            % (out) llh: a (1 X n) double of the log likelihood at each theta estimate.
            % (out) stats: A 1x1 struct of fitting statistics.
            image=obj.checkImage(image);
            if nargin==2
                estimator_name='Newton';
            end
            if nargin<4
                theta_init=[];
            end
            %Check to make sure we have a theta_init for each image
            nIms = size(image, ndims(obj.imsize)+1);
            theta_init=obj.checkThetaInit(theta_init, nIms);

            if ~ischar(estimator_name)
                error('MappelBase:estimate', 'Invalid estimation method name');
            elseif strcmpi(estimator_name,'GPUGauss')
                [theta, crlb, llh, stats]=obj.estimate_GPUGaussMLE(image);
            elseif strncmpi(estimator_name,'trust-region-reflective',5) %matches case insensitive 'trust*'
                [theta, crlb, llh, stats]=obj.estimate_fmincon(image, theta_init, 'trust-region-reflective');
            elseif strncmpi(estimator_name,'interior-point',8) %matches case insensitive 'interior*'
                [theta, crlb, llh, stats]=obj.estimate_fmincon(image, theta_init, 'interior-point');
            elseif ismember(lower(estimator_name),cellmap(@lower,obj.EstimationMethods))
                [theta, crlb, llh, stats]=obj.call('estimate',image, estimator_name, theta_init);
            else
                error('MappelBase:estimate', 'Invalid estimation method name');
            end
        end

        function [theta, crlb, llh, stats, sample, sample_llh]=estimateDebug(obj, image, estimator_name, theta_init)
            % [theta, crlb, llh, stats, sample, sample_rllh]=estimatedebug(image, name) - estimate theta's
            % crlb's and llh's for each image in stack.  
            % (in) image: a size:[imsizeY, imsizeX] image
            % (in) estimator_name: (optional) name for the optimization method. (default = 'Newton')  
            %           Valid names are in obj.EstimationMethods
            % (in) theta_init: (optional) Initial theta guess size (nParams x 1).  Values of 0 indicate
            %            that we have no initial guess and the estimator should form its own guess.
            % (out) theta: a (nParams X 1) double of estimated theta values
            % (out) crlb: a (nParams X 1) estimate of the CRLB for each parameter estimate.
            %             This gives the approximiate variance in the theta parameters
            % (out) llh: a (1 X 1) double of the log likelihood at each theta estimate.
            % (out) stats: A 1x1 struct of fitting statistics.
            % (out) sample: A (nParams X n) array of thetas that were searched as part of the maximization process
            % (out) sample_rllh: A (1 X n) array of relative log likelyhoods at each sample theta
            image=obj.checkImage(image);
            if nargin==2
                estimator_name='Newton';
            end
            if nargin<4
                theta_init=[];
            end
            %Check to make sure we have a theta_init for each image
            theta_init=obj.checkThetaInit(theta_init, 1);
            
            if ~ischar(estimator_name)
                error('MappelBase:estimateDebug', 'Invalid estimation method name');
            elseif strcmpi(estimator_name,'GPUGauss')
                [theta, crlb, llh, stats, sample, sample_llh]=obj.estimateDebug_GPUGaussMLE(image);
            elseif strcmpi(estimator_name,'fminsearch')
                [theta, crlb, llh, stats, sample, sample_llh]=obj.estimateDebug_fminsearch(image, theta_init);
            elseif strcmpi(estimator_name,'fminunc')
                [theta, crlb, llh, stats, sample, sample_llh]=obj.estimateDebug_fminunc(image, theta_init);
            elseif strncmpi(estimator_name,'trust-region-reflective',5) %matches case insensitive 'trust*'
                [theta, crlb, llh, stats, sample, sample_llh]=obj.estimateDebug_fmincon(image, theta_init,'trust-region-reflective');
            elseif strncmpi(estimator_name,'interior-point',8) %matches case insensitive 'interior*'
                [theta, crlb, llh, stats, sample, sample_llh]=obj.estimateDebug_fmincon(image, theta_init,'interior-point');
            elseif ismember(estimator_name,obj.EstimationMethods)
                [theta, crlb, llh, stats, sample, sample_llh]=obj.call('estimateDebug',image, estimator_name, theta_init);
            else
                error('MappelBase:estimateDebug', 'Invalid estimation method name');
            end
        end    

        function [mean, cov]=estimatePosterior(obj, image, max_samples, theta_init)
            % [mean, cov, stats]=estimatePosterior(obj, image, max_iterations)
            % (in) image: a double size:[imsizeY, imsizeX, n] image stack
            % (in) max_samples: A maximum number of samples to take
            % (in) theta_init: (optional) Initial theta guesses size (nParams x n).  Values of 0 indicate
            %            that we have no initial guess and the estimator should form its own guess.
            % (out) mean: a (nParams X n) double of estimated posterior mean values
            % (out) cov: a (nParams X nParams X n) estimate of the posterior covarience.
            % (out) stats: A 1x1 struct of fitting statistics.
            if nargin==2 || isempty(max_samples) || ~isfinite(max_samples) || max_samples<=1
                max_samples=int64(3000);
            elseif max_samples<10;
                error('MappelBase:estimatePosterior', 'Too few samples');
            else
                max_samples=int64(max_samples);
            end
            if nargin<4
                theta_init=[];
            end
            %Check to make sure we have a theta_init for each image
            nIms = size(image, ndims(obj.imsize)+1);
            theta_init=obj.checkThetaInit(theta_init, nIms);

            obj.checkCount(max_samples);
            image=obj.checkImage(image);
            [mean, cov]=obj.call('estimatePosterior',image, max_samples, theta_init);
        end

        function [mean, cov, sample, sample_llh, candidates, candidate_llh]=estimatePosteriorDebug(obj, image, max_samples, theta_init)
            % [mean, cov, stats]=estimatePosterior(obj, image, max_iterations)
            % (in) image: a double size:[imsizeY, imsizeX] image
            % (in) max_samples: A maximum number of samples to take
            % (in) theta_init: (optional) Initial theta guess size (nParams x 1).  Values of 0 indicate
            %            that we have no initial guess and the estimator should form its own guess.
            % (out) mean: a (nParams X 1) double of estimated posterior mean values
            % (out) cov: a (nParams X nParams X 1) estimate of the posterior covarience.
            % (out) stats: A 1x1 struct of fitting statistics.
            % (out) sample: A (nParams X nsamples) array of thetas samples
            % (out) sample_llh: A (1 X nsmaples) array of log likelyhoods at each sample theta
            % (out) candidates: A (nParams X nsamples) array of candidate thetas
            % (out) candidate_llh: A (1 X nsmaples) array of log likelyhoods at each candidate theta
            if nargin==2
                max_samples=int64(3000);
            elseif max_samples<10;
                error('MappelBase:estimatePosterior', 'Too few samples');
            else
                max_samples=int64(max_samples);
            end
            if nargin<4
                theta_init=[];
            end
            %Check to make sure we have a theta_init for each image
            theta_init=obj.checkThetaInit(theta_init, 1);

            obj.checkCount(max_samples);
            image=obj.checkImage(image);
            [mean, cov, sample, sample_llh, candidates, candidate_llh]=obj.call('estimatePosteriorDebug',image, max_samples, theta_init);
        end

        function [theta, crlb, llh, stats]=estimate_GPUGaussMLE(obj, image)
            if ~ispc()
                error('MappelBase:estimateGPUGaussMLE','Unable to run GPUGaussMLE on this archetecture');
            end
            if obj.psf_sigma(1)~=obj.psf_sigma(2) || obj.imsize(1)~=obj.imsize(2)
                error('MappelBase:estimateGPUGaussMLE','Unable to run GPUGaussMLE as boxsize or psf_sigma is not uniform');            
            end
            if obj.GPUGaussMLEFitType<1
                error('MappelBase:estimateGPUGaussMLE','Unable to run GPUGaussMLE on this model type');
            end
            data=single(image);
            psf=obj.psf_sigma(1);
            iters=obj.GPUGaussMLE_Iterations;
            type=obj.GPUGaussMLEFitType;
            [P, CRLB, LL]=gpugaussmlev2(data, psf, iters, type);
            theta=double(P');
            crlb = double(CRLB'); 
            
            crlb([1,2],:)=crlb([2,1],:); %swap dims
            theta([1,2],:) = theta([2,1],:); %Swap dims
            theta([1,2],:) = theta([1,2],:)+0.5; %Correct for 1/2 pixel

            llh = LL;
%             stats.time=t;
            stats.iterations=iters;
            stats.fittype=type;
        end
        
        function [theta, crlb, llh, stats, sequence, sequence_llh]=estimateDebug_GPUGaussMLE(obj, image)
            if ~ispc()
                error('MappelBase:estimateGPUGaussMLE','Unable to run GPUGaussMLE on this archetecture');
            end
            if obj.psf_sigma(1)~=obj.psf_sigma(2) || obj.imsize(1)~=obj.imsize(2)
                error('MappelBase:estimateGPUGaussMLE','Unable to run GPUGaussMLE as boxsize or psf_sigma is not uniform');            
            end
            if obj.GPUGaussMLEFitType<1
                error('MappelBase:estimateGPUGaussMLE','Unable to run GPUGaussMLE on this model type');
            end
            data=single(image);
            psf=obj.psf_sigma(1);
            iters=obj.GPUGaussMLE_Iterations;
            type=obj.GPUGaussMLEFitType;
            [P, CRLB, LL]=gpugaussmlev2(data, psf, iters, type);
            theta=double(P');
            crlb = double(CRLB'); 
            
            crlb([1,2],:)=crlb([2,1],:); %swap dims
            theta([1,2],:) = theta([2,1],:); %Swap dims
            theta([1,2],:) = theta([1,2],:)+0.5; %Correct for 1/2 pixel

            llh = LL;
%             stats.time=t;
            stats.iterations=iters;
            stats.fittype=type;
            sequence=theta;
            sequence_llh=LL;
        end
        
        function [theta, crlb, llh, stats, sequence, sequence_llh]=estimateDebug_fminsearch(obj, image, theta_init)
            max_iter=5000;
            sequence=zeros(obj.nParams, max_iter);
            function stop = output(theta, opt, state)
                sequence(:,opt.iteration+1)=theta;
                stop = strcmp(state,'done');
            end
            if isempty(theta_init)
                theta_init = obj.estimate(image,'Heuristic', theta_init);
            end
            objective=@(t) -obj.LLH(image, t);
            opts = optimset('MaxIter',max_iter, 'OutputFcn', @output);
            opts.Diagnostics='on';
            opts.Display='final-detailed';
            [theta, fval, flag, out] = fminsearch(objective, theta_init, opts);
            crlb=obj.CRLB(theta);
            llh=obj.LLH(image,theta);
            stats.method = out.algorithm;
            stats.iterations = out.iterations;
            stats.fval = fval;
            stats.flag = flag;
            sequence = sequence(:,1:out.iterations);
            sequence_llh = obj.LLH(image, sequence);
        end
        
        function [theta, crlb, llh, stats, sequence, sequence_llh]=estimateDebug_fminunc(obj, image, theta_init)
            max_iter=5000;
            sequence=zeros(obj.nParams, max_iter);
            if isempty(theta_init)
                theta_init = obj.estimate(image,'Heuristic', theta_init);
            end
            function stop = output(theta, opt, state)
                sequence(:,opt.iteration+1)=theta;
                stop = strcmp(state,'done');
            end
            objective = @(theta) deal(-obj.LLH(image,theta), -obj.modelGrad(image,theta), -obj.modelHessian(image,theta));
            opts = optimset('MaxIter',max_iter,'MaxFunEvals',max_iter, 'OutputFcn', @output);
            opts.Algorithm='quasi-newton';
            opts.Diagnostics='on';
            opts.Display='final-detailed';
            opts.GradObj='on';
            opts.Hessian='on';
            [theta, fval, flag, out, grad, H] = fminunc(objective, theta_init,opts);
            crlb=obj.CRLB(theta);
            llh=obj.LLH(image,theta);
            stats.method = out.algorithm;
            stats.iterations = out.iterations;
            stats.hessian = H;
            stats.grad = grad;
            stats.fval = fval;
            stats.flag = flag;
            sequence = sequence(:,1:out.iterations);
            sequence_llh = obj.LLH(image, sequence);
        end
        
        function [theta, crlb, llh, stats]=estimate_fmincon(obj, image, theta_init, algorithm)
            %
            % Uses matlab's fmincon to maximize the LLH for a stack of images.  Uses LLH grad and hessian calculations from C++
            % interface.
            % [in] image - a double size:[imsizeY, imsizeX] image
            % [in] theta_init - size:[nParams, 1] vector giving initial theta value.  Use: [] or omit to use Hueristric
            %                   starting point
            % [in] algorithm - string giving algorithm to choose.  Options: 'interior-point' or 'trust-region-reflective'
            %                  [default: 'interior-point']
            if nargin==3
                algorithm = 'interior-point';
            end
            if nargin==2 || isempty(theta_init)
                theta_init = obj.estimate(image,'Heuristic', theta_init);
            end
            max_iter=500;
            lb = obj.boundTheta(-ones(obj.nParams,1)); %lower bound on theta
            ub = obj.boundTheta(inf(obj.nParams,1)); %upper bound on theta
            opts = optimset('MaxIter',max_iter,'MaxFunEvals',max_iter);
            opts.GradObj='on';
            opts.Diagnostics='off';
            opts.Display='off';
            hess_objective = @(im,theta) -obj.modelHessian(im,theta);
            switch algorithm
                case 'interior-point'
                    opts.Algorithm='interior-point';
                    opts.Hessian='user-supplied';
                    objective =@(im,theta) deal(-obj.LLH(im,theta), -obj.modelGrad(im,theta));
                case 'trust-region-reflective'
                    opts.Algorithm='trust-region-reflective';
                    opts.Hessian='user-supplied';
                    objective = @(im,theta) deal(-obj.LLH(im,theta), -obj.modelGrad(im,theta), -obj.modelHessian(im,theta));
                otherwise
                    error('MappelBase:estimateDebug_fmincon','Unknown maximization method: "%s"', algorithm);
            end

            nIms = size(image, numel(obj.imsize)+1); %number of images to process
            theta=zeros(obj.nParams, nIms);
            niters = zeros(nIms,1);
            for n=1:nIms
                if numel(obj.imsize)==2
                    im=image(:,:,n);
                else
                    im=image(:,:,:,n);
                end
%                 ob = objective(im);
                opts.HessFcn = @(theta,~) hess_objective(im,theta);
                [theta(:,n), ~, ~, out, ~, ~] = fmincon(@(x) objective(im,x), theta_init(:,n), [],[],[],[], lb, ub,[], opts);
                niters(n) = out.iterations;
            end
            crlb=obj.CRLB(theta);
            llh=obj.LLH(image,theta);
            stats.method = out.algorithm;
            stats.iterations = mean(niters);
        end
        
        function [theta, crlb, llh, stats, sequence, sequence_llh]=estimateDebug_fmincon(obj, image, theta_init, algorithm)
            %
            % Uses matlab's fmincon to maximize the LLH for a single image giving debugging infdormation.
            % Uses LLH grad and hessian calculations from C++ interface.
            %
            % [in] image - a double size:[imsizeY, imsizeX] image
            % [in] theta_init - size:[nParams, 1] vector giving initial theta value.  Use: [] or omit to use Hueristric
            %                   starting point
            % [in] algorithm - string giving algorithm to choose.  Options: 'interior-point' or 'trust-region-reflective'
            %                  [default: 'interior-point']
            if nargin==3
                algorithm = 'interior-point';
            end
            if nargin==2 || isempty(theta_init)
                theta_init = obj.estimate(image,'Heuristic', theta_init);
            end
            max_iter=500;
            sequence=zeros(obj.nParams, max_iter); %sequence of points evaluated at
            lb = obj.boundTheta(-ones(obj.nParams,1)); %lower bound on theta
            ub = obj.boundTheta(inf(obj.nParams,1)); %upper bound on theta
            function stop = output(theta, opt, state)
                sequence(:,opt.iteration+1)=theta;
                stop = strcmp(state,'done');
            end
            opts = optimset('MaxIter',max_iter,'MaxFunEvals',max_iter, 'OutputFcn', @output);
            opts.Diagnostics='on';
            opts.Display='final-detailed';
            opts.GradObj='on';
            switch algorithm
                case 'interior-point'
                    opts.Algorithm='interior-point';
                    opts.Hessian='user-supplied';
                    opts.HessFcn=@(theta,lambda) -obj.modelHessian(image,theta);
                    objective = @(theta) deal(-obj.LLH(image,theta), -obj.modelGrad(image,theta));
                case 'trust-region-reflective'
                    opts.Algorithm='trust-region-reflective';
                    opts.Hessian='user-supplied';
                    objective = @(theta) deal(-obj.LLH(image,theta), -obj.modelGrad(image,theta), -obj.modelHessian(image,theta));
                otherwise
                    error('MappelBase:estimateDebug_fmincon','Unknown maximization method: "%s"', algorithm);
            end
            [theta, fval, flag, out, grad, H] = fmincon(objective, theta_init, [],[],[],[], lb, ub,[], opts);
            crlb=obj.CRLB(theta);
            llh=obj.LLH(image,theta);
            stats.method = out.algorithm;
            stats.iterations = out.iterations;
            stats.hessian = H;
            stats.grad = grad;
            stats.fval = fval;
            stats.flag = flag;
            sequence = sequence(:,1:out.iterations);
            sequence_llh = obj.LLH(image, sequence);
        end

        function stats=getStats(obj)
            % stats=obj.getStats() - Statistics and parameters of this model
            % (out) stats - 1x1 Struct describing this model.
            stats=obj.call('getStats');
        end
        
        %%Model Testing
        function [llh, theta_bg_mle] = uniformBackgroundModelLLH(obj, ims)
            % Test the model fit of a 1-parameter constant background model to the stack of images.
            % The mle estimate for a 1-parameter baground parameter is just the mean of the image.
            % The log-likelihood is calculated at this MLE estimate.
            % (in) ims: a double size:[imsizeY, imsizeX, n] image stack
            % (out) llh: a length N vector of the LLH for each image for the consant-background model
            % (out) theta_bg_mle: a length N vector of the estimated MLE constant background.
            npixels = prod(obj.imsize);
            ims = reshape(ims,npixels,[]);
            theta_bg_mle = mean(ims)';
            llh =  log(theta_bg_mle).*sum(ims)' - npixels*theta_bg_mle - sum(gammaln(ims+1))';
        end

        function [pass, LLRstat] = modelComparisonUniform(obj, alpha, ims, theta_mle)
            % Do a LLH ratio test to compare the emitter model to a single parameter constant background model
            % The images are provided along with the estimated theta mle values for the emitter model.
            % The LLH ratio test for nested models can be used and we compute the test statistice
            % LLRstat = -2*llh_const_bg_model + 2*llh_emitter_model.  This should be chisq distributed
            % with number of degrees of freedom given by obj.nParams-1 since the const bg model has 1
            % param.
            % (in) alpha: 0<=alpha<1 - the certainty with which we should be sure to accept an emitter fit
            %                          vs. the uniform background model.  Values close to 1 reject more
            %                          fits.  Those close to 0 accept most fits.  At 0 only models where
            %                          the constant bg is more likely (even though it has only 1 free parametet)
            %                          will be rejected.  These arguably should always be rejected.  It
            %                          indicates an almost certainly bad fit.
            % (in) ims: a double size:[imsizeY, imsizeX, n] image stack
            % (in) theta_mle: a double size:[nParams, n] sequence of theta MLE values for each image in
            %                 ims
            % (out) pass: a boolean length N vector which is true if the emitter model passes for this test.
            %             in otherwords it is true for images with good fits that should be kept.  Images
            %             that fail could just as easily have been just random noise.
            % (out) LLRstat: -2*log(null_model_LH / emitter_model_LH);
            assert(0<=alpha && alpha<1);
            if nargin<4
                theta_mle = obj.estimate(ims);
            end
            llh_emitter_model = obj.LLH(ims,theta_mle);
            llh_const_bg_model = obj.uniformBackgroundModelLLH(ims);
            LLRstat = 2*(llh_emitter_model-llh_const_bg_model);
            threshold = chi2inv(alpha,obj.nParams-1);
            pass = LLRstat > threshold;
        end

        function llh = noiseBackgroundModelLLH(obj, ims)
            % Test the model fit of an npixels-parameter all noise background model to the stack of images.
            % In this model each pixel has its own parameter and that pixels mle will of course be the
            % value of the pixel itself.  Unlike the constant bg model there is no point to return the
            % mle values themselves since they are just the images.
            % (in) ims: a double size:[imsizeY, imsizeX, n] image stack
            % (out) llh: a length N vector of the LLH for each image for the consant-background model
            npixels = prod(obj.imsize);
            ims = reshape(ims,npixels,[]);
            llh = sum(ims.*log(ims)-ims-gammaln(ims+1))';
        end

        function [pass, bg_prob] = modelComparisonNoise(obj, alpha, ims, theta_mle)
            if nargin<4
                theta_mle = obj.estimate(ims);
            end
            emitter_aic = 2*obj.nParams - 2*obj.LLH(ims,theta_mle);
            npixels = prod(obj.imsize);
            apparent_bg_aic = 2*npixels - 2*obj.noiseBackgroundModelLLH(ims);
            pass = emitter_aic < apparent_bg_aic;
            bg_prob = min(1,exp((emitter_aic - apparent_bg_aic)/2));
            pass = pass & (bg_prob<(1-alpha));
        end
        %% Evaluation and testing methods
    
        function [rmse, stddev, error] = evaluateEstimatorAt(obj, estimator, theta, nTrials, theta_init)
            %Evaluate this 2D estimator at a particular theta by running a fixed number of trials
            %
            % (in) estimator - Estimator name.  Can be any of the MAP estimator names or 'Posterior N' where N is a count
            % (in) theta - Parameter location to test at
            % (in) nTrials - Number of simulated images to fit
            % (out) rmse - [1,obj.nParams]: The root-mean-squared error
            if nargin<5
                theta_init=[];
            end
            if nargin<4
                nTrials=1000;
            end
            theta = theta(:);
            ims = obj.simulateImage(theta,nTrials);
            
            if strncmpi(estimator,'posterior',9)
                count = str2double(estimator(10:end));
                theta_est=obj.estimatePosterior(ims,count,theta_init);
            else
                theta_est=obj.estimate(ims,estimator, theta_init);
            end
            error = theta_est-repmat(theta,1,nTrials);
            rmse = sqrt(mean(error.*error,2));
            stddev = std(error')';
        end

        function [theta_est, est_var]=evaluateEstimatorOn(obj, estimator, images)
            %Evaluate this 2D estimator at a particular theta using the given samples which may have
            % been generated using different models or parameters.
            %
            % (in) estimator - String. Estimator name.  Can be any of the MAP estimator names or 'Posterior N' where N is a count
            % (in) images - size[imsizeY,imsize,N] -  An array of sample images to test on
            % (out) theta_est - size:[nParams,N]: the estimated thetas
            % (out) est_var - size:[nParams,N]: the estimated variance of the estimate at each theta
            if strncmpi(estimator,'posterior',9)
                count = str2double(estimator(10:end));
                [theta_est,est_cov]=obj.estimatePosterior(images,count);
                est_var=zeros(obj.nParams,size(est_cov,3));
                for i=1:size(est_cov,3)
                    est_var(:,i)=diag(est_cov(:,:,i));
                end
            else
                [theta_est,est_var]=obj.estimate(images,estimator);
            end
        end
      
        function [theta_est_grid,est_var_grid]=mapEstimatorAccuracy(obj,estimator, sample_grid)
            % (in) estimator - String. Estimator name.  Can be any of the MAP estimator names or 'Posterior N' where N is a count
            % (in) sample_grid - size:[sizeY,sizeX,nTrials,gridsizeX,gridsizeY]
            %                     This should be produced by makeThetaGridSamples()
            % (out) theta_est_grid - size:[nParams, nTrials, gridsizeX,gridsizeY] - estimated theta for each grid
            %                           image
            % (out) est_var_grid - size:[nParams, nTrials, gridsizeX,gridsizeY] - estimated variance at each theta
            nTrials=size(sample_grid,3);
            gridsize=[size(sample_grid,4), size(sample_grid,5)];
            theta_est_grid=zeros([obj.nParams, nTrials, gridsize]);
            est_var_grid=zeros([obj.nParams, nTrials, gridsize]);
            h=waitbar(0,sprintf('Maping Accuracy Model:%s Estimator%s gridsize:%s',obj.Name, estimator, mat2str(gridsize)));
            for x=1:gridsize(1)
                for y=1:gridsize(2)
                    [theta_est, est_var]=obj.evaluateEstimatorOn(estimator, sample_grid(:,:,:,x,y));
                    theta_est_grid(:,:,x,y)=theta_est;
                    est_var_grid(:,:,x,y)=est_var;
                end
                waitbar(x/gridsize(1),h);
            end
            close(h);
        end

        function [theta_grid,sample_grid]=makeThetaGridSamples(obj, theta, gridsize, nTrials)
            % Make a grid of many theta values and corresponding sampled
            % images to test on for a particular theta.  This allows the
            % testing of the sensitivity of the estimator to the absolute
            % position in the image.  This generates the input for
            % mapEstimatorAccuacy().
            %
            % (in) theta - A theta to test over a spatial grid.  The x an y
            % are irrelevent as they will be modified for each grid point.
            % (in) gridsize - [X Y] The size of the grid of test locations.
            %                   suggest: [30 30]
            % (out) theta_grid - size:[obj.nParams, nTrials,gridsize(1),gridsize(2)]
            % (out) sample_grid - size:[sizeY,sizeX,nTrials,gridsize(1),gridsize(2)]
            theta=theta(:);
            if isscalar(gridsize)
                gridsize=[gridsize gridsize];
            end
            theta_grid=zeros([obj.nParams, nTrials, gridsize]);
            sample_grid=zeros([obj.imsize, nTrials, gridsize]);
            grid_edges.x=linspace(0,obj.imsize(1),gridsize(1)+1);
            grid_edges.y=linspace(0,obj.imsize(2),gridsize(2)+1);
            for x=1:gridsize(1)
                for y=1:gridsize(2)
                    pixel_thetas=repmat(theta,1,nTrials);
                    e0=[grid_edges.x(x) grid_edges.y(y)]';
                    e1=[grid_edges.x(x+1) grid_edges.y(y+1)]';
                    %Make pixel thetas uniformly distributied over the pixel
                    pixel_thetas(1:2,:)=rand(2,nTrials).*repmat(e1-e0,1,nTrials)+repmat(e0,1,nTrials);
                    theta_grid(:,:,x,y)=pixel_thetas;
                    sample_grid(:,:,:,x,y)=obj.simulateImage(pixel_thetas);
                end
            end
        end
        
        %% Visualization methods
        
        function srim=superResolutionModel(obj, theta, theta_err, res_factor)
            % Generate a super-res image showing the estimated location and
            % error as a gaussian
            % (in) theta - The parameter value to visualize
            % (in) theta_err - [optional] The RMSE error or sqrt(CRLB) that represents
            %                   the localization error. [default = run a simulation to estimate error] 
            % (in) res_factor - [optinal] integer>1 - the factor to blow up the image
            %                   [default = 100]
            % (out) srim - double size:[sizeY*res_factor,sizeX*res_factor]
            %             A super-res rendering of the emitter fit postion.  
            if nargin<4
                res_factor=100;
            end
            if nargin<3
                theta_err = obj.evaluateEstimatorAt('Newton',theta);
            end
            srimsize=obj.imsize*res_factor;
            theta(1)=theta(1)*res_factor;
            theta(2)=theta(2)*res_factor;
            theta_err(1)=theta_err(1)*res_factor;
            theta_err(2)=theta_err(2)*res_factor;

            xs=0:srimsize(1)-1;
            X=0.5*(erf(((xs+1)-theta(1))/(sqrt(2)*theta_err(1)))-erf((xs-theta(1))/(sqrt(2)*theta_err(1))));
            ys=0:srimsize(2)-1;
            Y=0.5*(erf(((ys+1)-theta(2))/(sqrt(2)*theta_err(2)))-erf((ys-theta(2))/(sqrt(2)*theta_err(2))));
            srim=theta(3)*Y'*X;           
        end
        
        function plotAccuracyMap(obj, grid)
            x=[0.5,obj.imsize(1)-0.5];
            y=[0.5,obj.imsize(2)-0.5];
            imagesc(x,y,grid');
            xlabel('x (pixels)');
            ylabel('y (pixels)');
        end

    end% public methods
    
    methods (Static = true)
        function [C,D]=cholesky(A)
            [valid, C, d] = MappelBase.callstatic(@Gauss2DMLE_Iface,'cholesky',A);
            if(~valid)
                error('MappelBase:cholesky','Matrix not positive definite symmetric');
            end
            D = diag(d);
        end

        function [C,D, is_positive]=modifiedCholesky(A)
            [is_positive, C, d] = MappelBase.callstatic(@Gauss2DMLE_Iface,'modifiedCholesky',A);
            D = diag(d);
        end

        function [x, modified]=choleskySolve(A,b)
            [modified, x] = MappelBase.callstatic(@Gauss2DMLE_Iface,'choleskySolve',A,b);
        end

        function fig = viewDipImage(image, fig)
            if nargin==1
                fig=figure();
            end
            if ~isa(image,'dip_image')
                image = dip_image(image);
            end
            dipshow(fig,image);
            diptruesize(100*1000./max(size(image)))
            dipmapping(fig,'colormap',hot);
        end
    end %public static methods
    
    methods (Access=protected)        
        function image=checkImage(obj, image)
            ndim=length(obj.imsize);
            if (ndims(image) < ndim) || (ndims(image) > ndim+1)
                error('MappelBase:checkImage', 'Invalid image dimension');
            end
            if ndim==2
                if size(image,1)~=obj.imsize(2) || size(image,2)~=obj.imsize(1) 
                    error('MappelBase:checkImage', 'Invalid image shape');
                end
            else
                imsz=size(image);
                if any(imsz(ndim:-1:1)~=obj.imsize)
                    error('MappelBase:checkImage', 'Invalid image shape');
                end
            end
            if ~all(image>=0)
                error('MappelBase:checkImage', 'Invalid image');
            end
            image=double(image);
        end
        
         function count=checkCount(~, count)
            if ~isscalar(count) || count<=0
                error('MappelBase:checkCount', 'Invalid count');
            end
            count=int32(count);
        end

        function theta=checkTheta(obj, theta)
            if size(theta,1) ~= obj.nParams
                if length(theta) == obj.nParams
                    theta=theta';
                else
                    error('MappelBase:checkTheta', 'Invalid theta shape. Expected %i params',obj.nParams);
                end
            end
            inbounds=obj.thetaInBounds(theta);
            if ~inbounds
                theta=obj.boundTheta(theta);
            else
                theta=double(theta);
            end
        end
        
        function theta_init=checkThetaInit(obj, theta_init, nIms)
            if isempty(theta_init)
                return;
            elseif isvector(theta_init)
                if numel(theta_init)~=obj.nParams
                    error('Mappel:thetaInitValue','Invalid theta init shape');
                end
                theta_init=repmat(theta_init(:),1,nIms);
            elseif any(size(theta_init) ~= [obj.nParams, nIms])
                error('Mappel:thetaInitValue','Invalid theta init shape');
            end    
        end
        
        function mask=paramMask(obj, names)
            mask=zeros(obj.nParams,1);
            for n=1:obj.nParams;
                param=obj.ParamNames{n};
                mask(n)=any(cellfun(@logical, cellfun(@(s) strcmp(s,param), names, 'Uniform', 0)));
            end
        end
    end % protected methods

end %classdef
