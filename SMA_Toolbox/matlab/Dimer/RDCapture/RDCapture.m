% File: RDcapture.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% 08/2016
% RDCapture 
%
% Computes survival probabilities for diffusion to capture problem for a circular trap at the origin.
%
% Uses matlab or a high-speed C++ implementation for computation

classdef RDCapture < IfaceMixin
    % RDCapture
    % Liklihood calculations and simulations for diffusion-to-capture with the Erban-Chapman model
    % 
    %
    properties
        D; %Mutual diffusion constant (um^2/s)
        rho;  %Binding radius (um)
        lambda; %Capture rate when within binding radius (1/s)
    end

    properties (Access=protected)
        GSorder=7; % Order of Gaver-Stehfest inverse laplace t-form approx
        GSfactors; % Qfactors for Gaver-Stehfest algorithm
        GSsigns;
    end
    
    properties (Access = protected, Transient=true)
        initialized=false; %True once the object is correctly initialized by the initializeTrack method
    end

    methods
        function obj=RDCapture(Dmutual, rho, lambda)
            % [in]
            %  Dmutual - Mutual diffusion constant (um^2/s)
            %  rho - Binding radius (um)
            %  lambda - Capture rate when within target (um^2/s)
            obj=obj@IfaceMixin(@RDCapture_Iface);
            obj.initialized=true;
            obj.D = Dmutual;
            obj.rho = rho;
            obj.lambda = lambda;
            obj.GSfactors=obj.generateGaverStefestFactors(obj.GSorder);
            obj.initialized = obj.openIface(obj.D, obj.rho, obj.lambda);
        end
        
        function solvePDE(obj)
            model= createpde();
            max_r = 3*obj.rho;
            %Geometry description matrix
            shapes = [1 0 0 obj.rho zeros(1,6); 3 4 -max_r -max_r max_r max_r -max_r max_r max_r -max_r]';
            shape_error_status = csgchk(shapes); % This checks the Geometry Description matrix
            assert(~any(shape_error_status));

            names = char('CaptureSite','Chamber')';
            set_formula = 'Chamber+CaptureSite';
            [geometry, shape_table] = decsg(shapes,set_formula,names);
            figure();
            pdegplot(geometry,'EdgeLabels','on')
            xlim([-max_r,max_r])
            axis equal
            model.geometryFromEdges(geometry);
            model.applyBoundaryCondition('neumann','edge',1:4,'g',0);
            model.specifyCoefficients('m',0,'d',1,'c',obj.D,'a',0,'f',0);
            model.generateMesh();
            initFun = @(p) 1/(2*pi*obj.rho/3)^2*exp(-.5*(p.x.^2+(p.y-3*obj.rho).^2)/(obj.rho/3));
            model.setInitialConditions(initFun);
            ts = linspace(0,1,100);
            R = model.solvepde(ts);
            colormap('hot');
            for k=1:numel(ts)
                pdeplot(model,'XYdata',R.NodalSolution(:,k),'ZData',R.NodalSolution(:,k),'Mesh','on');
                title(sprintf('time: %.3g',ts(k)));
                zlim([0,1]);
                drawnow();
            end
        end

        function plotSurivalProbability(obj,ts)
            ts=sort(ts);
            Nsimulate = 1e4;
            max_dt = min(ts(1), obj.rho^2/(2*obj.D))/10;
            fprintf('MaxSimulation dt:%.5g\n',max_dt);
            Nr0 = 100;
            Nt = numel(ts);
            max_r0 = 10*obj.rho;
            r0 = linspace(0,max_r0,Nr0);
            Q=zeros(Nt,Nr0);
            QC=zeros(Nt,Nr0);
            Qring=zeros(Nt,Nr0);
            Qpring=zeros(Nt,Nr0);
            Qinf=zeros(Nt,Nr0);
            for k = 1:Nt
                tsk = repmat(ts(k),1,Nr0);
                Q(k,:)=obj.computeSurvivalProb(r0, tsk);
                QC(k,:)=obj.computeSurvivalProb_C(r0, ts(k));
                Qring(k,:)=obj.computeSurvivalProb_AbsorbingRing(r0,tsk);
                Qpring(k,:)=obj.computeSurvivalProb_PermeableRing(r0,tsk);
                Qinf(k,:)=obj.computeSurvivalProb_InfiniteTrap(tsk);
            end
            Qsamp = obj.simulate(r0,ts,Nsimulate,max_dt);
%               Qring(:)=min(max(Qring(:),0),1);
%             Qpring(:)=min(max(Qpring(:),0),1);
            figure();
            title(sprintf('D: %.3g rho: %.3g lambda: %.4g',obj.D,obj.rho,obj.lambda));
            for k=1:Nt
                subplot(Nt,1,k);
                hold();
                plot(r0,Q(k,:),'k:','DisplayName','Disk');
                plot(r0,QC(k,:),'k-','DisplayName','Disk [C++]');
                plot(r0,Qring(k,:),'r-','DisplayName','Abs. Ring');
                plot(r0,Qpring(k,:),'r--','DisplayName','Per. Ring');
                plot(r0,Qinf(k,:),'b-','DisplayName','Infinte trap');
                plot(r0,Qsamp(k,:),'ko','MarkerSize',3,'DisplayName','Simulated');
                ylim([0,1]);
                legend('location','best');
                xlabel('Intial dist: r0');
                ylabel('Survival Prob: Q');
                title(sprintf('Time: %.3g',ts(k)));
            end
        end

        function plotSurivalProbability3D(obj)
            Nr0 = 100;
            Nt = 100;
            max_r0 = sqrt(10)*obj.rho;
            r0 = linspace(0,max_r0,Nr0);
            max_time = 5*max_r0^2/obj.D;
            t = logspace(-3,log(max_time), Nt);
            Q=zeros(Nt,Nr0);
            Qring=zeros(Nt,Nr0);
            Qpring=zeros(Nt,Nr0);
            Qinf=zeros(Nt,Nr0);
            for k = 1:Nr0
                r0k = repmat(r0(k),1,Nt);
                Q(:,k)=obj.computeSurvivalProb(r0k,t);
                Qring(:,k)=obj.computeSurvivalProb_AbsorbingRing(r0k,t);
                Qpring(:,k)=obj.computeSurvivalProb_PermeableRing(r0k,t);
                Qinf(:,k)=obj.computeSurvivalProb_InfiniteTrap(t);
            end
            Qring(:)=min(max(Qring(:),0),1);
            Qpring(:)=min(max(Qpring(:),0),1);
            
            figure();
            [X,Y] =meshgrid(r0,t);
            hold('on');
            h=surface(X,Y,Q,'DisplayName','Survival Probability Disk','EdgeColor','None');
            h=surface(X,Y,Qpring,'DisplayName','Survival Probability Permeable Ring','EdgeColor','None');
            h=surface(X,Y,Qring,'DisplayName','Survival Probability Ring','EdgeColor','None');
            h=surface(X,Y,Qinf,'DisplayName','Survival Probability Infinite','EdgeColor','None');
            colorbar();
            ax=gca();
            ax.TickDir='out';
            ax.Box='on';
            ax.BoxStyle='Full';
            grid('on');
            grid('minor');
            xlabel('Distance r0 (um)');
            ylabel('Time (s)');
            ax.YScale='log';
%             ax.YDir='reverse';
            ylim([min(Y(:)),max(Y(:))]);
            xlim([min(X(:)),max(X(:))]);
            title(sprintf('D:%.3g um^2/s rho:%.3g um lambda:%.3g 1/s',obj.D,obj.rho,obj.lambda));
%             legend('location','best');
            axis('vis3d');
        end

        function Q=computeSurvivalProb_C(obj, rs, ts)
            % Compute the value of survival probablity Q(t,r0) in C
            % [in]
            %   rs - size:K vector of times (or scalar) to evaluate at
            %   ts - size:K vector of distances (or scalar) to evaluate at
            % [out]
            %   Q - size:K vector of values of survival probability Q
            if (numel(ts)>1 && numel(rs)>1) && (numel(ts) ~= numel(rs))
                error('RDCapture:BadParameters','r0 and t must have the same length if both are vectors.');
            end
            Q = obj.call('survivalProb',rs,ts);
        end
        
        function Q=computeSurvivalProb(obj, rs, ts)
            % Compute the value of survival probablity Q(t,r0;D,rho,lambda)
            % [in]
            %   rs - size:K vector of times (or scalar) to evaluate at
            %   ts - size:K vector of distances (or scalar) to evaluate at
            %         If both rs and ts are non-scalar, they must have the same length
            % [out]
            %   mu - size:[K,L] vector of values of mu.  mu(i,j) is time i distance j.
            NR=numel(rs);
            NT=numel(ts);
            K=max(NR,NT);
            if NR>1 && NT>1 && NR~=NT
                error('RDCapture:compute_mu','If rs and ts are both vectors they must have the same length');
            elseif NR>1 && NT==1
                ts=repmat(ts,K,1);
            elseif NR==1 && NT>1
                rs=repmat(rs,K,1);
            end
            Q = ones(K,1);
            potential_capture = rs-obj.rho < 6*sqrt(2*obj.D*ts);
            Q(potential_capture) = arrayfun(@(r,t) obj.laplaceInverse(@(s) obj.computeLT_Q(r,s(:),obj.D,obj.rho,obj.lambda),t),...
                                            rs(potential_capture), ts(potential_capture));
        end

        function Q=computeSurvivalProb_AbsorbingRing(obj,rs,ts)
            function val = LT_Q(r,s) 
                omega = sqrt(s/obj.D);
%                 v1=besselk(0,omega*r);
%                 v2=besselk(0,omega*obj.rho);
%                 fprintf('Param: %.8g, v1:%.8g v2:%.8g Ratio: %.8g\n',omega(end)*obj.rho,v1(end),v2(end),v1(end)/v2(end));
                val=(1./s).*(1-besselk(0,omega*r)./besselk(0,omega*obj.rho));
            end
            Q = arrayfun(@(r,t) obj.laplaceInverse(@(s) LT_Q(r,s),t), rs,ts);
        end

        function Q=computeSurvivalProb_PermeableRing(obj,rs,ts)
            function Qt=LT_Q(r,s)
                omega = sqrt(s/obj.D);
                xi = 1/(2*pi*obj.D);
                Qt=(1./s).*(1-(xi*besseli(0,omega*obj.rho).*besselk(0,omega*r))./(1/obj.lambda + xi*besseli(0,omega*obj.rho).*besselk(0,omega*obj.rho)));
            end
            Q = arrayfun(@(r,t) obj.laplaceInverse(@(s) LT_Q(r,s),t),rs,ts);
        end

        function Q=computeSurvivalProb_InfiniteTrap(obj,ts)
            Q = exp(-obj.lambda*ts);
        end


        function Q=compute_Q(obj,rs, ts, D, rho,lambda)
            % Compute the value of survival probablity Q(t,r0;D,rho,lambda)
            % [in]
            %   rs - size:K vector of times (or scalar) to evaluate at
            %   ts - size:K vector of distances (or scalar) to evaluate at
            %         If both rs and ts are non-scalar, they must have the same length
            %   D - [optional] Mutual diffusion constant (um^2/s) [default = obj.D]
            %   rho - [optional] Capture radius (um)              [default = obj.rho]
            %   lambda - [optional] Capture rate (1/s)            [default = obj.lambda]
            % [out]
            %   mu - size:[K,L] vector of values of mu.  mu(i,j) is time i distance j.
            if nargin<4; D=obj.D; end
            if nargin<5; rho=obj.rho; end
            if nargin<6; lambda=obj.lambda; end
            NR=numel(rs);
            NT=numel(ts);
            K=max(NR,NT);
            if NR>1 && NT>1 && NR~=NT
                error('RDCapture:compute_mu','If rs and ts are both vectors they must have the same length');
            elseif NR>1 && NT==1
                ts=repmat(ts,K,1);
            elseif NR==1 && NT>1
                rs=repmat(rs,K,1);
            end
            Q=zeros(K,1);
            for k=1:K
                Q(k) = obj.laplaceInverse(@(s) obj.computeLT_Q(rs(k),s,D,rho,lambda), ts(k));
            end
        end

        function nu=compute_nu(obj, t,D,rho)
            % Compute the value of nu(t;D,rho)
            % [in]
            %   t - size:K vector of times to evaluate at
            %   D - [optional] Mutual diffusion constant (um^2/s) [default = obj.D]
            %   rho - [optional] Capture radius (um)              [default = obj.rho]
            % [out]
            %   nu - size:K vector of values of nu at each time
            if nargin<3; D=obj.D; end
            if nargin<4; rho=obj.rho; end
            nu = obj.laplaceInverse(@(s) obj.computeLT_nu(s(:),D,rho), t);
        end

        function nu = compute_nu_C(obj, ts)
            % Compute the value of nu(t) in C
            % [in]
            %   ts - size:K vector of distances (or scalar) to evaluate at
            % [out]
            %   nu - size:K vector of values of nu(t)            
            nu = obj.call('nu',ts);
        end

        function mu = compute_mu_C(obj, rs, ts)
            % Compute the value of mu(r0,t) in C
            % [in]
            %   rs - size:K vector of times (or scalar) to evaluate at
            %   ts - size:K vector of distances (or scalar) to evaluate at
            % [out]
            %   mu - size:K vector of values of mu(t)
            if (numel(ts)>1 && numel(rs)>1) && (numel(ts) ~= numel(rs))
                error('RDCapture:BadParameters','r0 and t must have the same length if both are vectors.');
            end
            mu = obj.call('mu',rs,ts);
        end

        function mu=compute_mu(obj, rs, ts, D, rho)
            % Compute the value of mu(t,r0;D,rho)
            % [in]
            %   rs - size:K vector of times (or scalar) to evaluate at
            %   ts - size:K vector of distances (or scalar) to evaluate at
            %         If both rs and ts are non-scalar, they must have the same length
            %   D - [optional] Mutual diffusion constant (um^2/s) [default = obj.D]
            %   rho - [optional] Capture radius (um)              [default = obj.rho]
            % [out]
            %   mu - size:K vector of values of mu(t)
            if nargin<4; D=obj.D; end
            if nargin<5; rho=obj.rho; end
            NR=numel(rs);
            NT=numel(ts);
            K=max(NR,NT);
            if NR>1 && NT>1 && NR~=NT
                error('RDCapture:compute_mu','If rs and ts are both vectors they must have the same length');
            elseif NR>1 && NT==1
                ts=repmat(ts,K,1);
            elseif NR==1 && NT>1
                rs=repmat(rs,K,1);
            end
            mu=zeros(K,1);
%             mu2=zeros(K,1);
            for k=1:K
                mu(k) = obj.laplaceInverse(@(s) obj.computeLT_mu(rs(k),s,D,rho), ts(k));
%                 mu2(k) = (1-marcumq(rs(k)./ sqrt(2*D*ts(k)),rho./sqrt(2*D*ts(k))));
            end
        end
        
        function survivalProb = simulate(obj, r0, ts, N, max_dt)
            % [in]
            %  r0 - initial radius
            %  ts - measurment times
            %  N - number of particles
            %  max_dt - maximum time step to take
            survivalProb = obj.call('simulate',r0, ts, uint32(N), max_dt);
        end

        function val = I0(obj, x)
            val = obj.call('I0',x);
        end
        function val = I1(obj, x)
            val = obj.call('I1',x);
        end
        function val = K0(obj, x)
            val = obj.call('K0',x);
        end
        function val = K1(obj, x)
            val = obj.call('K1',x);
        end
        function val = logI0(obj, x)
            val = obj.call('logI0',x);
        end
        function val = logI1(obj, x)
            val = obj.call('logI1',x);
        end
        function val = logK0(obj, x)
            val = obj.call('logK0',x);
        end
        function val = logK1(obj, x)
            val = obj.call('logK1',x);
        end
        function val = I0K0(obj, varargin)
            
            val = obj.call('I0K0',varargin{:});
        end
        function val = I0K1(obj, varargin)
            val = obj.call('I0K1',varargin{:});
        end
        function val = I1K0(obj, varargin)
            val = obj.call('I1K0',varargin{:});
        end
        function val = I1K1(obj, varargin)
            val = obj.call('I1K1',varargin{:});
        end
        function testBessel(obj)
            N=1e2;
            xs=logspace(-6,6,N);
            i0 = obj.I0(xs);
            i1 = obj.I1(xs);
            k0 = obj.K0(xs);
            k1 = obj.K1(xs);
            i0m = besseli(0,xs)';
            i1m = besseli(1,xs)';
            k0m = besselk(0,xs)';
            k1m = besselk(1,xs)';
            figure();
            ax=subplot(2,2,1);
            hold();
            ax.XScale='log';
            plot(xs,i0,'k-','DisplayName','I0');
            plot(xs,i1,'r-','DisplayName','I1');
            plot(xs,k0,'b-','DisplayName','K0');
            plot(xs,k1,'g-','DisplayName','K1');
            plot(xs,i0m,'ko','DisplayName','I0m');
            plot(xs,i1m,'ro','DisplayName','I1m');
            plot(xs,k0m,'bo','DisplayName','K0m');
            plot(xs,k1m,'go','DisplayName','K1m');
            legend('location','best');
            ylim([0,10]);
            
            ax=subplot(2,2,2);
            logi0 = obj.logI0(xs);
            logi1 = obj.logI1(xs);
            logk0 = obj.logK0(xs);
            logk1 = obj.logK1(xs);
            hold();
            ax.XScale='log';
            plot(xs,logi0,'k-','DisplayName','logI0');
            plot(xs,logi1,'r-','DisplayName','logI1');
            plot(xs,logk0,'b-','DisplayName','logK0');
            plot(xs,logk1,'g-','DisplayName','logK1');
            plot(xs,log(i0m),'ko','DisplayName','logI0m');
            plot(xs,log(i1m),'ro','DisplayName','logI1m');
            plot(xs,log(k0m),'bo','DisplayName','logK0m');
            plot(xs,log(k1m),'go','DisplayName','logK1m');
            legend('location','best');
            ax=subplot(2,2,3);
            hold();
            ax.XScale='log';
            i0k0 = obj.I0K0(xs);
            i0k1 = obj.I0K1(xs);
            i1k0 = obj.I1K0(xs);
            i1k1 = obj.I1K1(xs);

            plot(xs,i0k0,'k-','DisplayName','I0K0');
            plot(xs,i0k1,'r-','DisplayName','I0K1');
            plot(xs,i1k0,'b-','DisplayName','I1K0');
            plot(xs,i1k1,'g-','DisplayName','I1K1');
            plot(xs,i0m.*k0m,'ko','DisplayName','I0K)m');
            plot(xs,i0m.*k1m,'ro','DisplayName','I0K1m');
            plot(xs,i1m.*k0m,'bo','DisplayName','I1K0m');
            plot(xs,i1m.*k1m,'go','DisplayName','I1K1m');
            ylim([0,10]);
            legend('location','best');

            ax=subplot(2,2,4);
            hold();
            ys = ones(N,1)*1e1;
            i0k0 = obj.I0K0(xs,ys);
            i0k1 = obj.I0K1(xs,ys);
            i1k0 = obj.I1K0(xs,ys);
            i1k1 = obj.I1K1(xs,ys);
            k0y = besselk(0,ys);
            k1y = besselk(1,ys);
            ax.XScale='log';
            plot(xs,i0k0,'k-','DisplayName','I0K0');
            plot(xs,i0k1,'r-','DisplayName','I0K1');
            plot(xs,i1k0,'b-','DisplayName','I1K0');
            plot(xs,i1k1,'g-','DisplayName','I1K1');
            plot(xs,i0.*k0y,'ko','DisplayName','I0K0m');
            plot(xs,i0.*k1y,'ro','DisplayName','I0K1m');
            plot(xs,i1.*k0y,'bo','DisplayName','I1K0m');
            plot(xs,i1.*k1y,'go','DisplayName','I1K1m');
            legend('location','best');
            ylim([0,1e10]);
            ax.YScale='log';

        end
    end

    methods (Static=true)        
        function QLT = computeLT_Q(r0,s,D,rho,lambda)
            % Compute the Laplace Transform of survival probability Q(s,r0;D,rho,lambda)
            % [in]
            %   r0 - initial displacement (um)
            %   s - size:K vector of laplace-space coordinates
            %   D - Mutual diffusion constant (um^2/s)
            %   rho - Capture radius (um)
            %   lambda - Capture rate (1/s)
            % [out]
            %   QLT - size:K vector of values of nu laplace transform at each s value
            QLT = (1-RDCapture.computeLT_mu(r0,s,D,rho)./(1/lambda + RDCapture.computeLT_nu(s,D,rho)))./s;
        end

        function nuLT=computeLT_nu(s,D,rho)
            % Compute the Laplace Transform of nu(s;D,rho)
            % [in]
            %   s - size:K vector of laplace-space coordinates
            %   D - Mutual diffusion constant (um^2/s)
            %   rho - Capture radius (um)
            % [out]
            %   nuLT - size:K vector of values of nu laplace transform at each s value
            Z= rho*sqrt(s/D);
            nuLT=(1-2*besseli(1,Z).*besselk(1,Z))./s ;
        end

        function muLT=computeLT_mu(r0,s,D,rho)
            % Compute the Laplace Transform of mu(s,r0;D,rho)
            % [in]
            %   r0 - initial displacement (um)
            %   s - size:K vector of laplace-space coordinates
            %   D - Mutual diffusion constant (um^2/s)
            %   rho - Capture radius (um)
            % [out]
            %   muLT - size:K vector of values of mu laplace transform at each s value
            omega = sqrt(s/D);
            if r0<rho % start inside trap
                muLT = (1-rho*omega.*besseli(0,omega*r0).*besselk(1,omega*rho))./s;
            else
                muLT = (rho*omega.*besseli(1,omega*rho).*besselk(0,omega*r0))./s;
            end
        end
    end
    methods (Access=protected)
        function vals = laplaceInverse(obj,F, t)
            % Numerical laplace inverse with the Gaver-Stehfest algorithm
            % [in]
            %  F - a function handle that accepts vectors of s laplace space coordinates and retrurns
            %      vectors of the function value at each laplace coordinate
            %  ts - size:k array of times to evaluate at
            K=numel(t);
            vals=zeros(K,1);
            idxs = (1:2*obj.GSorder)';
            alphas = log(2)./t;
            for k=1:K
                vals(k) = alphas(k) * sum(obj.GSfactors.*F(alphas(k)*idxs));
            end
        end

    end
    methods (Static=true,Access=protected)
        function gsFactors = generateGaverStefestFactors(order)
            % Numerical laplace inverse with the Gaver-Stehfest algorithm
            gsFactors=zeros(2*order,order);
            for k=1:2*order
                for j = floor((k+1)/2):min(k,order)
                    gsFactors(k,j) = j^(order+1) / (factorial(j)^2*factorial(order-j)*factorial(k-j)*factorial(2*j-k)) * factorial(2*j);
                end
            end
            alternating_signs = ones(2*order,1);
            alternating_signs(mod(order,2)+1:2:end)=-1;
            gsFactors = alternating_signs.*sum(gsFactors,2);
            gsFactors = gsFactors(:);
        end
    end
end
