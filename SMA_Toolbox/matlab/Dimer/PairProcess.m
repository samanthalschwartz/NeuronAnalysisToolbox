classdef PairProcess
    properties (Constant=true)
        log2pi = log(2*pi);
    end
    
    properties
        params = struct( 'DA',0.1,... % [um^2/s]
                    'DB',0.1,... % [um^2/s]
                    'DAB',0.01,... % [um^2/s]
                    'Dtheta',2*pi*1e-1,... % [rad^2/s]
                    'rho',0.050,... % [um]
                    'sigma_rho',0.015,... % [um]
                    'gamma',1e-2... % [1/s]
                  );
        prior = struct( 'state',[.5, .5], ...
                        'posType','Gauss', ... Options: ['Rect','Circ','Gauss']
                        'posCenter',[0,0],... Gauss:center; Circ:center; Rect:LowerLeftCorner
                        'posExtent',[10 0; 0 10] ... 'Gauss':Covariance, 'Circ':radius, 'Rect':area
                );
    end

    methods
        function pos = samplePositionPrior(obj, N, prior)
            % Sample from the initial position prior
            % [in] N - number of samples
            % [in] prior - [optional] prior struct to override local obj.prior.  [default=obj.prior]
            % [out] pos - 2xN vector of initial positions
            if nargin<3
                prior=obj.prior; %#ok<*NASGU>
            end
            switch prior.posType
                case {'Rect','rect'}
                    pos = repmat(prior.posCenter(:),1,N) + [rand(1,N)*prior.posExtent(1); rand(1,N)*prior.posExtent(2)];
                case {'Circ','circ'}
                    rs = sqrt(prior.posExtent*rand(1,N));
                    thetas = rand(1,N)*2*pi;
                    pos = repmat(prior.posCenter(:),1,N) + [rs.*cos(thetas); rs.*sin(thetas)];
                case {'Gauss','gauss'}
                    pos = repmat(prior.posCenter(:),1,N) + chol(prior.posExtent)*randn(2,N);
            end
        end

        function [ys,bound] = sampleTrajectories(obj, ts, N,  params, prior)
            % [in]
            %  ts - Kx1 time scheudle
            %  N - number of tracks to simulate
            %  params - [optional] params struct [default: obj.params]
            %  prior - [optional] prior struct [default: obj.prior]
            % [out]
            %  ys - 4xKxN matrix of trajectories
            %  bound - 1xN matrix of states (0=Free 1=Bound)
            if nargin<4; params=obj.params; end
            if nargin<5; prior=obj.prior; end
            K=numel(ts);
            ys=zeros(4,K,N);

            [ys_init, bound] = obj.sampleInitialPosition(N, params, prior);
            Nbound = sum(bound);
            Nfree = N-Nbound;
            OU=OUProcess(params.rho,params.sigma_rho,params.gamma,params.rho);
            ys(:,1,:) = ys_init;
            ys(:,1,bound) = PairProcess.transV(squeeze(ys(:,1,bound)));
            for k=2:K
                dt = ts(k) - ts(k-1);
                %Simulate free
                ys(1:2,k,~bound) = randn(2,Nfree).*(2*params.DA*dt)+squeeze(ys(1:2,k-1,~bound));
                ys(3:4,k,~bound) = randn(2,Nfree).*(2*params.DB*dt)+squeeze(ys(3:4,k-1,~bound));
                %Simulate bound in V space
                ys(1:2,k,bound) = randn(2,Nbound).*(2*params.DAB*dt)+squeeze(ys(1:2,k-1,bound));          % c_x c_y
                ys(4,k,bound) = mod(randn(1,1,Nbound).*(2*params.Dtheta*dt)+ys(4,k-1,bound),2*pi); % theta
            end
            %OU model for r parameter for bound model
            ys(3,:,bound) = OU.sampleGillespie(ts,Nbound,ys(3,1,bound));
            for k=1:K %Transform bound to Y space
                ys(:,k,bound) = PairProcess.transVinv(squeeze(ys(:,k,bound))); % Transform back to Y space
            end
        end

        function [ys,bound] = sampleInitialPosition(obj, N, params, prior)
            % Sample N initial states according to the prior.state for the initial state and the
            % prior.posType for the position
            % [in]
            %  N - number of tracks to simulate
            %  params - [optional] params struct [default: obj.params]
            %  prior - [optional] prior struct [default: obj.prior]
            % [out]
            %  ys - 4xN matrix of initial positions
            %  bound -1xN matrix of states (0=Free 1=Bound)
            if nargin<3; params=obj.params; end
            if nargin<4; prior=obj.prior; end
            ys=zeros(4,N);
            bound= rand(1,N)>prior.state(1);
            Nbound = sum(bound);
            Nfree = N - Nbound;
            % Simulate bound initial positions
            %Simulate in V space
            ys(1:2, bound) = obj.samplePriorPos(Nbound,prior);         % c_x c_y
            ys(3,bound) = randn(1,Nbound)*params.sigma_rho+params.rho; % r
            ys(4,bound) = rand(1,Nbound)*2*pi;                         % theta
            ys(:,bound) = PairProcess.transVinv(ys(:,bound)); % Transform back to Y space
            % Simulate free initial positions
            ys(1:2, ~bound) = obj.samplePriorPos(Nfree,prior);
            ys(3:4, ~bound) = obj.samplePriorPos(Nfree,prior);
        end

        function plotSamples(obj, ts, ys, bound)
            ts=ts';
            figure();
            ax=axes();
            hold();
            colormap('lines');
            xmin = min(reshape(ys([1,3],:,:),1,[]));
            xmax = max(reshape(ys([1,3],:,:),1,[]));
            ymin = min(reshape(ys([2,4],:,:),1,[]));
            ymax = max(reshape(ys([2,4],:,:),1,[]));
            tmin = ts(1); tmax=ts(end);
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            zlim([tmin,tmax]);
            grid('on'); grid('minor');
            view(-45,30);
            asp(1) = (xmax-xmin)/(ymax-ymin);
            asp(2) = 1;
            asp(3) = mean([asp(1),asp(2)]);
            pbaspect(asp);
            xlabel('X (um)');
            ylabel('Y (um)');
            zlabel('T (s)');
            ax.TickDir='out';
            ax.Box='on';
            ax.BoxStyle='Full';
            for n=1:size(ys,3)
                Y = ys(:,:,n);
                N = size(Y,2);
                C = n * ones(N,1);
                lw = 1 + 2*~bound(n);
                h = surface(Y([1,1],:)',Y([2,2],:)',[ts,ts],[C,C],'LineWidth',lw,'EdgeColor','interp','CDataMapping','Direct');
                h = surface(Y([3,3],:)',Y([4,4],:)',[ts,ts],[C,C],'LineWidth',lw,'EdgeColor','interp','CDataMapping','Direct');
            end            
        end

        function [xs, xVar] = simulateObservation(obj, ys, params, prior)
            
        end

        function llh=initialObsLLH(obj, ys, params, prior)
            % LLH of the intial observations y0 in each column of y0s under our prior distribution
            % [in]
            %  ys - 4xN or 4xKxN vector of positions in y-space for each of N pairs of length K
            %  params - [optional] params struct [default: obj.params]
            %  prior - [optional] prior struct [default: obj.prior]
            % [out]
            %   llh  - 2xN first row=freeLLH second=boundLLH. 
            if nargin<3; params=obj.params; end
            if nargin<4; prior=obj.prior; end
            N=size(ys,ndims(ys));
            if ndims(ys)==3
                v0s = PairProcess.transV(squeeze(ys(:,1,:)));
            else
                v0s = PairProcess.transV(ys);
            end
            rs=v0s(3,:);
            llh = zeros(2,N);
            switch prior.posType
                case {'Rect','rect'}
                    llh(1,:) = log(PairProcess.pairR_PDF_rect(rs,prior.posExtent(1),prior.posExtent(2)));
                case {'Circ','circ'}
                    llh(1,:) = log(PairProcess.pairR_PDF_circ(rs,prior.posExtent));
                case {'Gauss','gauss'}
                    llh(1,:) = log(PairProcess.pairR_PDF_2Dgauss(rs,prior.posExtent));
            end
            llh(2,:) = max(-1e9,log(normpdf(rs,params.rho, params.sigma_rho)));
        end

        function llh=trajectoryLLH(obj, ts, ys, params, prior)
            % This computes the model likelihoods for the free and bound models for
            % a given trajectory in Y-space where rows of ys are [ax ay bx by]'
            % [in]
            %  ts - 1xK vector of times
            %  ys - 4xKxN vector of positions in y-space for each of N pairs of length K
            %  params - [optional] params struct [default: obj.params]
            %  prior - [optional] prior struct [default: obj.prior]
            % [out]
            %   llh  - 2xKxN first row=freeLLH second=boundLLH for each time point
            if nargin<4; params=obj.params; end
            if nargin<5; prior=obj.prior; end
            N=size(ys,3);
            K=size(ys,2);
            vs = reshape(PairProcess.transV(reshape(ys,4,[])),size(ys));
            llh = zeros(2,K,N);
            llh(:,1,:) = obj.initialObsLLH(ys, params, prior);
            for k=2:K % for each time
                dt = ts(k)-ts(k-1);
                %Calculate free LLH
                vA = 2*params.DA*dt;
                vB = 2*params.DB*dt;
                sigmaF = diag([vA,vA,vB,vB]);
                llh(1,k,:) = sum(PairProcess.multiGaussianLLH(squeeze(ys(:,k,:)-ys(:,k-1,:)), sigmaF));
                %Calculate bound LLH
                vAB = 2*params.DAB*dt;
                vTheta = 2*params.Dtheta*dt;
                sigmaB = diag([vAB,vAB,vTheta]);
                llh(2,k,:) = sum(PairProcess.multiGaussianLLH( squeeze(vs([1,2,4],k,:)-vs([1,2,4],k-1,:)), sigmaB));
                muR = params.rho + (vs(3,k-1,:)-params.rho)*exp(-params.gamma*dt);
                sR = sqrt(params.sigma_rho^2/(2*params.gamma)*(1-exp(-2*params.gamma*dt)));
                llh(2,k,:) = llh(2,k,:) + log(normpdf(vs(3,k,:), muR, sR))-log(vs(3,k,:));
            end
        end

        function plotTrajectoryLLH(obj, ts, ys, bound, params, prior)
            if nargin<5; params=obj.params; end
            if nargin<6; prior=obj.prior; end
            llh=trajectoryLLH(obj, ts, ys, params, prior);
            figure();
            K=size(llh,3);
            N=size(llh,2);
            ax=subplot(2,1,1);
            hold(ax);
            for k=1:K
                if bound(k)
                    style='b-';
                else
                    style='r-';
                end
                plot(ts,llh(1,:,k),style)
            end
            title('Free Model');
            xlabel('Observation time');

            ax=subplot(2,1,2);
            hold(ax);
            for k=1:K
                if bound(k)
                    style='b-';
                else
                    style='r-';
                end
                plot(ts,llh(2,:,k),style)
            end
            ylim([-1000,100]);
            title('Bound Model');
            xlabel('Observation time');

        end

        function [state, BF] = classify(obj, ts, Yvec)

        end

        function estimateBayesError(obj, ts, params, prior)
        end

        function plotBayesError(obj, tmax, params, prior)
        end
    end
    
    methods (Static=true)
        function  H=plot_CDF_match(samp,xs,cdfs)
            H=figure();
            [F,X] = ecdf(samp);
            plot(X,F,'b-','DisplayName','Sample empircal cdf');
            hold();
            plot(xs,cdfs,'r-','DisplayName','Computed cdf');
            ylabel('Probability');
            legend('Location','Best');
        end


        function  H=plot_PDF_match(samp,xs,pdfs)
            H=figure();
%             [F,X] = ecdf(samp);
%             [V,C] = ecdfhist(F,X);
%             plot(C,V,'b-','DisplayName','Sample empircal pdf');
            hold();
            plot(xs,pdfs,'r-','DisplayName','Computed pdf');
            H=histogram(samp,'Normalization','pdf');
            H.DisplayName='Sample pdf';
            ylabel('Probability');
            legend('Location','Best');
        end

        function plotR_CDF_line(a)
            N=1e5;
            rs = linspace(0,a,N);
            p_rs = PairProcess.pairR_CDF_line(rs,a);
            dist = abs(rand(1,N) - rand(1,N))*a;
            PairProcess.plot_CDF_match(dist, rs, p_rs);
            xlabel('PairDistance: r');
            title(sprintf('Pair Distance - Uniform Line [w:%.3g]',a));
        end

        function plotR_PDF_line(a)
            N=1e6;
            rs = linspace(0,a,N);
            p_rs = PairProcess.pairR_PDF_line(rs,a);
            dist = abs(rand(1,N) - rand(1,N))*a;
            PairProcess.plot_PDF_match(dist, rs, p_rs);
            xlabel('PairDistance: r');
            title(sprintf('Pair Distance - Uniform Line [w:%.3g]',a));
        end


        function plotR_PDF_rect(w,h)
            N=1e6+1;
            rs = linspace(0,sqrt(w^2+h^2),N);
            samp = PairProcess.sampleRect(N,0,0,w,h);
            p_rs = PairProcess.pairR_PDF_rect(rs,w,h);
            dist = sqrt(sum((samp(:,1:floor(N/2)) - samp(:,ceil(N/2)+1:end)).^2));
            PairProcess.plot_PDF_match(dist, rs, p_rs);
            xlabel('PairDistance: r');
            title(sprintf('Pair Distance - Uniform Rect [w:%.3g, h:%.3g]',w,h));
        end

        function plotR_PDF_circ(rho)
            N=1e6;
            rs = linspace(0,2*rho,N);
            samp = PairProcess.sampleCirc(N,0,0,rho);
            p_rs = PairProcess.pairR_PDF_circ(rs,rho);
            dist = sqrt(sum((samp(:,1:floor(N/2)) - samp(:,ceil(N/2)+1:end)).^2));
            PairProcess.plot_PDF_match(dist, rs, p_rs);
            xlabel('PairDistance: r');
            title(sprintf('Pair Distance - Uniform Circ [radius:%.3g]',rho));
        end


        function plotR_PDF_gauss(Sigma)
            N=1e6;
            samp = PairProcess.sampleGauss(N,[0,0],Sigma);
            dist = sqrt(sum((samp(:,1:floor(N/2)) - samp(:,ceil(N/2)+1:end)).^2));            
            rs = linspace(0,max(dist)*1.1,N);
            p_rs = PairProcess.pairR_PDF_2Dgauss(rs,Sigma);
            PairProcess.plot_PDF_match(dist, rs, p_rs);
            xlabel('PairDistance: r');
            title(sprintf('Pair Distance - Bivariate Gauss [%.3g, %.3g; %.3g, %.3g]',Sigma(1,1),Sigma(1,2),Sigma(2,1),Sigma(2,2)));
            figure();
            plot(samp(1,:),samp(2,:),'k.','MarkerSize',2);
            axis('equal');
        end

        function plot_norm_PDF_2Dgauss(mu, Sigma)
            N=1e6;
            samp = PairProcess.sampleGauss(N,mu,Sigma);
            dist = sqrt(sum(samp.^2));            
            rs = linspace(0,max(dist)*1.1,N);
            pdf_rs = PairProcess.norm_PDF_2Dgauss(rs,mu,Sigma);
            PairProcess.plot_PDF_match(dist, rs, pdf_rs);
            xlabel('PairDistance: r');
            title(sprintf('Bivariate Gauss Norm: \\mu=[%.3g; %.3g] \\Sigma=[%.3g, %.3g; %.3g, %.3g]',mu(1),mu(2),Sigma(1,1),Sigma(1,2),Sigma(2,1),Sigma(2,2)));
            figure();
            plot(samp(1,:),samp(2,:),'k.','MarkerSize',2);
            axis('equal');
        end


        function plotR_InitialModelComparison()
            Sigma = .01*diag([1,1]);
            w=2; h=2;
            radius=1;
            rho=0.050;
            sigma_rho=0.015;
            N=1000;
            rs=logspace(-3,2,N);
            p_rect = PairProcess.pairR_PDF_rect(rs,w,h);
            p_circ = PairProcess.pairR_PDF_circ(rs,radius);
            p_gauss = PairProcess.pairR_PDF_2Dgauss(rs,Sigma);
            p_bound = normpdf(rs,rho,sigma_rho);
            figure();
            subplot(2,1,1);
            plot(rs,p_bound,'k-','DisplayName','Bound');
            hold();
            plot(rs,p_rect,'r-','DisplayName','Rect');
            plot(rs,p_circ,'m-','DisplayName','Disk');
            plot(rs,p_gauss,'b-','DisplayName','Gauss');
            ax=gca();
            ax.YScale='log';
            ax.YLim = [1e-3,3e1];
            ax.XScale='log';
            H=legend('location','best');
%             H.interpreter='latex';
            xlabel('Relative distance $r$');
            ylabel('Probability');
            subplot(2,1,2);
            hold('on');
            p_R_given_z0_gauss = [p_bound',p_gauss'];
            p_R_and_z0_gauss = p_R_given_z0_gauss .* repmat(.5,N,2);
            p_z0_gauss = p_R_and_z0_gauss ./ repmat(sum(p_R_and_z0_gauss,2),1,2);
            
            p_R_given_z0_rect = [p_bound',p_rect'];
            p_R_and_z0_rect = p_R_given_z0_rect .* repmat(.5,N,2);
            p_z0_rect = p_R_and_z0_rect ./ repmat(sum(p_R_and_z0_rect,2),1,2);

            p_R_given_z0_circ = [p_bound',p_circ'];
            p_R_and_z0_circ = p_R_given_z0_circ .* repmat([.9,.1],N,1);
            p_z0_circ = p_R_and_z0_circ ./ repmat(sum(p_R_and_z0_circ,2),1,2);
            plot(rs,p_z0_rect(:,1),'r-','DisplayName','p(z0=Bound | prior=rect)');
            plot(rs,p_z0_rect(:,2),'r--','DisplayName','p(z0=Free | prior=rect)');
            plot(rs,p_z0_circ(:,1),'m-','DisplayName','p(z0=Bound | prior=disk)');
            plot(rs,p_z0_circ(:,2),'m--','DisplayName','p(z0=Free | prior=disk)');
            plot(rs,p_z0_gauss(:,1),'b-','DisplayName','p(z0=Bound | prior=gauss)');
            plot(rs,p_z0_gauss(:,2),'b--','DisplayName','p(z0=Free | prior=gauss)');
            ax=gca();
            ax.XScale='log';
            ax.YScale='log';
            ax.YLim = [1e-4,1];
            plot([rho,rho],ax.YLim,'r:','DisplayName','Mean binding radius rho');
            ylabel('p(Bound|prior_model)');
            xlabel('r0 [um]');
            legend('location','best');
        end
        
%         function plotR_PDF_rect(w,h)
%             N=1e6;
%             rs = linspace(0,sqrt(max(w,h)),N);
%             samp = PairProcess.sampleRect(N,0,0,w,h);
%             figure();
%             histogram(pdist(samp'),'Normalization','pdf');
%             hold();
%             p_rs = PairProcess.pairR_PDF_rect(rs,w,h);
%             plot(rs,p_rs,'r-');
%         end

        function pos = samplePriorPos(N, prior)
            switch prior.posType
                case {'Gauss','gauss'}
                    pos = PairProcess.sampleGauss(N,prior.posCenter,prior.posExtent);
                case {'Circ','circ','disk'}
                    pos = PairProcess.sampleCirc(N,prior.posCenter(1),prior.posCenter(2),prior.posExtent);
                case {'Rect','rect'}
                    pos = PairProcess.sampleRectr(N,prior.posCenter(1),prior.posCenter(2),prior.posExtent(1),prior.posExtent(2));
            end
        end

        function pos = sampleRect(N,l,b,w,h)
            pos = [rand(1,N)*w+l; rand(1,N)*h+b];
        end

        function pos = sampleGauss(N,mu,Sigma)
            s=chol(Sigma);
            n = size(Sigma,1);
            pos = repmat(mu(:),1,N) + (randn(n,N)'*s)';
        end

        function pos = sampleCirc(N,x,y,radius)
            r = sqrt(radius^2*rand(1,N));
            theta = 2*pi*rand(1,N);
            pos = [x+r.*cos(theta); y+r.*sin(theta)];
        end

        function p_rs = pairR_PDF_line(rs,a)
            filter = rs<a;
            p_rs=zeros(size(rs));
            p_rs(filter) = rs(filter).*(2*(1/a^2)*(a./rs(filter) - 1));
        end


        function p_rs = pairR_CDF_line(rs,a)
            ss = rs.^2;
            filter = rs<=a;
            p_rs=zeros(size(rs));
            p_rs(filter) = (1/a^2)*(2*a*rs(filter) - ss);
        end

        function p_rs = pairR_PDF_rect(rs,w,h)
            % rs - list of distances
            % w - width of rect
            % h - height of rect
            % [out]
            %  p_rs - pdf for distance r between 2 points drawn from rect with sizes: [w,h]
            a=min(w,h);
            b=max(w,h);
            ss=rs.^2; % squares
            A=a*b; % area
            filter = rs <= a;
            p_rs=zeros(size(rs));
            p_rs(filter) = (2*rs(filter)).* ((1/A^2)*(A*pi+ss(filter)-2*(a+b)*rs(filter)));
            
            filter = (a < rs) & (rs <= b);
            p_rs(filter) = (2*rs(filter)).* ((1/A^2)*(-a^2 + 2*A*asin(a./rs(filter)) + 2*b*sqrt(ss(filter)-a^2) - 2*b*rs(filter)));
            filter = (b < rs);
            p_rs(filter) = (2*rs(filter)).* ((1/A^2)*(-a^2 -b^2 - A*pi -ss(filter) +...
                                                      2*A*(asin(a./rs(filter))+asin(b./rs(filter))) + ...
                                                      2*b*sqrt(ss(filter)-a^2) + 2*a*sqrt(ss(filter)-b^2) ));
        end
        
        function p_r = pairR_PDF_circ(r,a)
            theta = 2*acos(r/(2*a));
            p_r = (2*r).*(theta - sin(theta))/(pi*a^2);
        end

        
        function p_r = pairR_PDF_2Dgauss(r,Sigma)
            % [in]
            %   r - vector Nx1 of distances to evaluate pdf at
            %   Sigma - 2x2 positive definite symmetric matrix
            % [out]
            %  p_r - vector Nx1 of pdf vals for each r value
            [~,p]=chol(Sigma);
            if p>0 || ~issymmetric(Sigma) % Make sure we are positive definite and symmetric
                error('PairProcess:pairR_PDF_2Dgauss','Invalid covariance matrix');
            end
            lambda = 2*eig(Sigma);
            omega = sum(lambda);
            q2=lambda(2)/lambda(1);
            q=sqrt(q2);
            q2inv = 1/q2;
            q4=q2*q2;
            p_r = (1/(q*omega))*(1+q2)*r.*exp(-(1+q2)^2/(4*q2*omega)*r.^2).*besseli(0,(1-q4)/(4*q2*omega)*r.^2);
        end


%        function p_r = norm_PDF_2Dgauss(r,mu,Sigma)
%             % [in]
%             %   r - vector Nx1 of distances to evaluate pdf at
%             %   mu - mean vector 2x1
%             %   Sigma - 2x2 positive definite symmetric matrix
%             % [out]
%             %  p_r - vector Nx1 of pdf vals for each r value
%             [~,p]=chol(Sigma);
%             if p>0 || ~issymmetric(Sigma) % Make sure we are positive definite and symmetric
%                 error('PairProcess:pairR_PDF_2Dgauss','Invalid covariance matrix');
%             end
%             lambda = eig(Sigma);
%             rp = norm(mu);
%             p_r = sqrt(prod(lambda))* (r.*exp(-(r.^2+rp^2)/2).*besseli(0,rp*r));
%         end



        function cdfs= gaussianCDF(ds,V)
            % ds - size[1,N] distances
            % V - variance
            cdfs = .5*(1+erf(ds./sqrt(2*V)));
        end

        function pdfs=gaussianPDF_interval(ds,V)
            % ds - size[2,N] ds(1,:) interval start, ds(2,:) interval end Note: all(ds(1,:) < ds(2,:))
            % V - variance scalar
            pdfs = gaussianCDF(ds(2,:),V)-gaussianCDF(ds(1,:),V);
        end

        function llhs=multiGaussianLLH(ds,Sigma)
            % [in]
            % ds: size:[k,N] each row is a distance vector in R^k
            % SigmaL size:[k,k]
            % [out]
            % llhs: size:[k,N] Log-likelihoods along each dimension
            llhs=-.5*(size(ds,1)*PairProcess.log2pi + log(det(Sigma)) + ds.*(Sigma\ds));
        end

        function vs = transV(ys)
            % Implement the Y->V transform
            % [in]
            %  ys - [4xN] rows [ax,ay,bx,by]
            % [out]
            %  vs - [4xN] rows [cx,cy,r,theta]
            vs = [ .5*(ys(1:2,:) + ys(3:4,:)); ...
                   sqrt(sum((ys(1:2,:) - ys(3:4,:)).^2)); ...
                   mod(atan2(ys(2,:) - ys(4,:),ys(1,:) - ys(3,:)),2*pi)]; 
        end

        function ys = transVinv(vs)
            % Implement the V->Y transform
            % [in]
            %  vs - [4xN] rows [cx,cy,r,theta]
            % [out]
            %  ys - [4xN] rows [ax,ay,bx,by]
            ys = [vs(1,:) + .5*vs(3,:).*cos(vs(4,:));...
                  vs(2,:) + .5*vs(3,:).*sin(vs(4,:));...
                  vs(1,:) - .5*vs(3,:).*cos(vs(4,:));...
                  vs(2,:) - .5*vs(3,:).*sin(vs(4,:))];
        end

        function us = transU(ys)
            % Implement the Y->U transform
            % [in]
            %  ys - [4xN] rows [ax,ay,bx,by]
            % [out]
            %  us - [4xN] rows [cx,cy,d_x,d_y]
            us = [.5*(ys(1:2,:) + ys(3:4,:)); ys(1:2,:) - ys(3:4,:)];
        end

        function ys = transUinv(us)
            % Implement the U->Y transform
            % [in]
            %  us - [4xN] rows [cx,cy,d_x,d_y]
            % [out]
            %  ys - [4xN] rows [ax,ay,bx,by]
            ys = [us(1:2,:) + .5*us(3:4,:); us(1:2,:) - .5*us(3:4,:)];
        end
    end

end
