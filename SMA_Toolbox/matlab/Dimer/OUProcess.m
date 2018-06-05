classdef OUProcess < handle
    %This does simulation and estimation of the OU process
    
    properties
        mu;    % Long-run mean        
        sigma; % Long-run variance
        gamma; % Relaxation rate
        x0;    % Starting position
    end
    properties (Dependent=true)
        D; % The effective diffusion constant
    end

    methods
        function obj=OUProcess(mu, sigma, gamma, x0)
            obj.mu = mu;
            obj.sigma = sigma;
            obj.gamma = gamma;
            if nargin<4
                obj.x0=mu;
            else
                obj.x0 = x0;
            end
        end

        function XS=sampleEM(obj,times,K)
            N=numel(times);
            XS=randn(N,K);
            XS(1,:)=obj.x0;
            dt=diff(times);
            S = sqrt((obj.sigma^2*2*obj.gamma)*dt);
            for n=2:N
                XS(n,:) = XS(n-1,:) + obj.gamma*(obj.mu-XS(n-1,:))*dt(n-1) + S(n-1)*XS(n,:);
            end
        end

        function XS=sampleGillespie(obj,times,K,x0)
            % [in]
            %  times: size:[N] vector of times to sample ate
            %  K: numer of sequences to simulate
            %  x0: [optional] starting point [default=obj.x0].  Can be vector size:[K] or scalar.
            % [out]
            %  XS: size: [N, K] - each column is a sample trajectory each row corresponds to times
            N=numel(times);
            XS=randn(N,K);
            if nargin<4
                XS(1,:)=obj.x0;
            else
                XS(1,:)=x0;
            end
            dt=diff(times);
            G = exp(-obj.gamma*dt);
            S = obj.sigma*sqrt(1-G.^2);
            for n=2:N
                XS(n,:) = (1-G(n-1))*obj.mu + G(n-1)*XS(n-1,:) + S(n-1)*XS(n,:);
            end
        end

        function XS=sampleGillespieEquilibrium(obj,times,K)
            % XS - size: [Ntimes, K] - each column is a sample
            N=numel(times);
            XS=randn(N,K);
            XS(1,:)=obj.mu+obj.sigma*randn(1,K);
            dt=diff(times);
            G = exp(-obj.gamma*dt);
            S = obj.sigma*sqrt(1-G.^2);
            for n=2:N
                XS(n,:) = (1-G(n-1))*obj.mu + G(n-1)*XS(n-1,:) + S(n-1)*XS(n,:);
            end
        end

        function llh = sampleLLH(obj,times,XS,varargin)
            % [in]
            %  times [N] vector of times
            %  XS [N x K] matrix of K column samples
            %  x0 - [optional] if samples start at x0 this should be set [default: x0=mu]
            %  t0 - [optional] if initial distribution starts at time t0 [default: t0=0]
            % [out]
            %  llh [N x K] matrix of llh of each time for each sample
            llh = obj.computeLLH(times,XS,obj.mu,obj.sigma,obj.gamma,varargin{:});
        end

        function llh = equilibriumLLH(obj,xs)
            % Probability of samples drawn indpenendenly from equillibrium distribution for process
            % The equilibrium distribution is normal arround mu with variance sigma^2
            % [in]
            %  xs - [Nx1] A vector of x values (iid)
            % [out]
            %  llh - [Nx1] Vector of log likelihoods
            V= obj.sigma^2;
            llh = -.5*(log(2*pi) + log(V) + (xs-obj.mu)./V);
        end

        function [M,V] = sampleGaussianDistribution(obj, times, XS, varargin)
            % [in]
            %  times [N] vector of times
            %  XS [N x K] matrix of K column samples
            %  x0 - [optional] if samples start at x0 this should be set [default: x0=mu]
            %  t0 - [optional] if initial distribution starts at time t0 [default: t0=0]
            % [out]
            %  M [N x K] matrix of distribution mean for each time for each sample
            %  V [N x K] matrix of distribution variance for each time for each sample
            [M,V] = obj.computeGaussianDistribution(times,XS,obj.mu,obj.sigma,obj.gamma,varargin{:});
        end

        
        function [mle_params, mle_llh] = estimateMLE(obj, times, XS)
            % [in]
            %  times [N] vector of times
            %  XS [N x K] matrix of K column samples
            % Params: [mu sigma gamma x0]
            function Z=llhFunc(P,Data)
                llh = OUProcess.computeLLH(Data.times,Data.XS,P(1),P(2),P(3),P(4));
                Z = sum(llh(:));
            end
            Data.times = times;
            Data.XS = XS;
            param_init_val = [mean(XS(:)),sqrt(var(XS(numel(XS)/2:end))),1,mean(XS(1,:))];
            mle = MLE(@llhFunc, Data,  param_init_val);
            mle.params.lowerBound =  [-1e6,1e-6,1e-6,-1e6];
            mle.params.upperBound = [1e6,1e6,1e6,1e6];
            mle.params.names = {'mu','sigma','gamma','x0'};
            mle.estimate();
            mle.viewProfileLikelihood();
            mle_params=mle.mle.param_est;
            mle_llh = mle.mle.llh;
        end

        function annalyzeEquilibrium(obj,XS,times)
            if nargin<3
                times = 1:size(XS,1);
            end
            figure();
            subplot(2,1,2);
            hold('on');
            K=size(XS,2);
            plot(times,std(XS,0,2),'DisplayName','Standard deviation of samples');
            plot([0,max(times)],obj.sigma*[1,1],'k:','DisplayName','Long-run std dev $\sigma$');
            plot(log(2)/obj.gamma*[1,1],ylim(),'r-','DisplayName','Relaxation Half-life $\ln{2}\gamma^{-1}$');
            H=legend('Location','Best');
            H.Interpreter='latex';
            xlabel('Time');
            ylabel('STD');
            title('Standard deviations');
            subplot(2,1,1);
            hold('on');
            for k=1:min(K,5)
                plot(times,XS(:,k),'DisplayName',sprintf('Sample %i',k));
            end
            plot([0,max(times)],obj.mu*[1,1],'k:','DisplayName','Long-run mean $\mu$');
            plot(log(2)/obj.gamma*[1,1],ylim(),'r-','DisplayName','Relaxation Half-life $\ln{2}\gamma^{-1}$');
            H=legend('Location','Best');
            H.Interpreter='latex';
            xlabel('time');
            ylabel('R');
            title('OU Sequences');
        end
    end

    methods
        function D = get.D(obj)
            D=sigma^2*2*gamma;
        end
    end

    methods (Static=true)
       function llh = computeLLH(times,XS,mu,sigma,gamma,varargin)
            % compute the sample LLH for a either an equillbrium start point or a 
            % fixed iniital point x0 at time t0
            % [in]
            %  times - size:[N] vector of times
            %  XS - size:[N x K] matrix of K column samples
            %  mu - long-run mean
            %  sigma - long-run standard deviation
            %  gamma - relaxationrate
            %  x0 - [optional] initial position defaults to equilibrium position mu
            %  t0 - [optional] initial time defaults to 0
            % [out]
            %  llh [N x K] matrix of llh of each time for each sample
            [M,V] = OUProcess.computeGaussianDistribution(times,XS,mu,sigma,gamma,varargin{:});
            llh=zeros(size(XS));
            if V(:,1)==0
                llh(1,:) = 0;
            else
                llh(1,:) = -.5*(log(2*pi) + log(V(1,:)) + (XS(1,:) - M(1,:)).^2./V(1,:));
            end
            llh(2:end,:) = -.5*(log(2*pi) + log(V(2:end,:)) + (XS(2:end,:) - M(2:end,:)).^2./V(2:end,:));
        end
    
        function [m,V]=singleStepGaussianDistribution(t,mu,sigma,gamma,x0)
            % [in]
            %  t - time
            %  mu - long-run mean
            %  sigma - long-run standard deviation
            %  gamma - relaxationrate
            %  x0 - [optional] initial position defaults to equilibrium position mu ncan be a vector or a
            %                   scalar.  If a fector should be size:[K]
            if nargin<5
                %Use equilibrium distribution
                m = mu;
                V = sigma^2;
            else
                %Observed x0 at time 0
                E=exp(-gamma*t);
                m = mu*(1-E)+x0*E;
                V = sigma^2*(1-E^2);
            end
        end

        function [M,V] = computeGaussianDistribution(times,XS,mu,sigma,gamma,x0,t0)
            % compute the sample gaussian distributions (mean and variance) an equillbrium start point or a 
            % fixed iniital point x0 at time t0
            % [in]
            %  times - size:[N] vector of times
            %  XS - size:[N x K] matrix of K column samples
            %  mu - long-run mean
            %  sigma - long-run standard deviation
            %  gamma - relaxationrate
            %  x0 - [optional] initial position defaults to equilibrium position mu ncan be a vector or a
            %                   scalar.  If a fector should be size:[K]
            %  t0 - [optional] initial time defaults to 0
            % [out]
            %  M [N x K] matrix of mean value of distribution for each time for each sample
            %  V  [N x K] matrix of variance of distribution for each time for each sample

            dt = diff(times);
            N=numel(times);
            K=size(XS,2);
            M = zeros(N,K);
            V = zeros(N,K);
            if nargin<6
                %Start from equilibrium distribution
                M(1,:)=mu;
                V(1,:)=sigma^2;
            else
                %Start from X(t0)=x0
                if nargin<7; t0=0; end  %t0 defaults to 0
                if times(1)<=0 
                    if  any(XS(1,:) ~= x0)
                        error('OUProcess:computeLLH','Invalid time sequence for a fixed starting point sequence.  Must have time(1)>0');
                    else
                        V(1,:) = 0;
                        M(1,:) = mu; % Essentially treat initial distribution as a delta, and try not to contribute to llh by adding 0.
                    end
                else
                    V(1,:) = sigma^2*(1-exp(-2*gamma*(times(1)-t0)));
                    E = exp(-gamma*(times(1)-t0));
                    M(1,:) = mu*(1-E) + x0*E;
                end
            end
            for n=2:N
                V(n,:) = sigma^2*(1-exp(-2*gamma*dt(n-1)));
                E = exp(-gamma*dt(n-1));
                M(n,:)= mu*(1-E)+XS(n-1,:)*E;
            end
        end
    end
end
