
classdef UnscentedTransform < handle
    properties
        kappa=0;
        alpha=0.5;
        beta=2;
    end
    methods
        function obj = UnscentedTransform(kappa_,alpha_,beta_)
            if nargin>=1
                obj.kappa=kappa_;
            end
            if nargin>=1
                obj.alpha=alpha_;
            end
            if nargin>=1
                obj.beta=beta_;
            end
        end

        function [est_mean,est_cov, sigmaP, sigmaP_tform] = estimate(obj, func, orig_mean, orig_cov)
            [sigmaP,meanW, covW] = obj.sigmaPoints(orig_mean, orig_cov);
            Nx = size(sigmaP,1);
            Nsigma = size(sigmaP,2);
            sigmaP_tform = func(sigmaP);
            est_mean = sum(sigmaP_tform.*repmat(meanW,Nx,1),2);
            delta = sigmaP_tform - repmat(est_mean,1,Nsigma);
            est_cov = delta*diag(covW)*delta';
%             figure();
%             subplot(2,1,1);
%             hold();
%             plot(sigmaP(1,:),sigmaP(2,:),'ro','DisplayName','sigmaP a');
%             plot(sigmaP(3,:),sigmaP(4,:),'bo','DisplayName','sigmaP b');
%             axis('equal');
%             legend('location','best');
%             subplot(2,1,2);
%             hold();
%             plot(sigmaP_tform(1,:),sigmaP_tform(2,:),'ro','DisplayName','sigmaP tform a');
%             plot(sigmaP_tform(3,:),sigmaP_tform(4,:),'bo','DisplayName','sigmaP tform b');
%             axis('equal');
%             legend('location','best');

        end

        function [sigmaP, meanW, covW] = sigmaPoints(obj, x_mean, x_cov)
            Nx = size(x_mean,1);
            lambda = obj.alpha^2*(Nx+obj.kappa)-Nx;
            Nsigma = Nx*2+1; % number of sigma points
            V = chol((Nx+lambda)*x_cov);
            sigmaP = repmat(x_mean,1,Nsigma);
            sigmaP(:,2:Nx+1) = sigmaP(:,2:Nx+1)+V;
            sigmaP(:,Nx+2:end) = sigmaP(:,Nx+2:end)-V;
            centerMeanWeight = lambda/(Nx+lambda);
            centerCovWeight = lambda/(Nx+lambda)+(1-obj.alpha^2+obj.beta);
            generalWeight = 1/(2*(Nx+lambda));
            meanW =[centerMeanWeight, repmat(generalWeight,1,Nsigma-1)];
            covW =[centerCovWeight, repmat(generalWeight,1,Nsigma-1)];
        end

        function testEstimate2D(obj, f, o_mean, o_cov)
            if nargin==1
                o_mean = [0.050; pi/3];
                o_cov = diag([0.010^2, (pi/8)^2]);
                f = @(xs) [xs(1,:).*cos(xs(2,:)); xs(1,:).*sin(xs(2,:))]; % Convert from radial coordinates
            end
            [e_mean, e_cov, sigmaP, f_sigmaP] = obj.estimate(f, o_mean, o_cov);
            o_cov_chol = chol(o_cov);
            o_ellipse = obj.covEllipse(o_mean, o_cov_chol);
            o3_ellipse = obj.covEllipse(o_mean, o_cov_chol,3);
            Nsample=1e4;
            o_sample = obj.sampleGaussian(o_mean,o_cov_chol,Nsample);
            f_mean = f(o_mean);
            f_ellipse = f(o_ellipse);
            f3_ellipse = f(o3_ellipse);
            f_sample = f(o_sample);
            
            e_cov_chol = chol(e_cov);
            e_ellipse = obj.covEllipse(e_mean, e_cov_chol);
            e3_ellipse = obj.covEllipse(e_mean, e_cov_chol,3);

                        
            figure();
            hold('on');
            plot(o_sample(1,:), o_sample(2,:),'r.','MarkerSize',1,'DisplayName','Orig Sample');
            plot(o_mean(1),o_mean(2),'s','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean');
            plot(o_ellipse(1,:),o_ellipse(2,:),'r-','LineWidth',3,'DisplayName','Orig 1-sigma Cov');
            plot(o3_ellipse(1,:),o3_ellipse(2,:),'r-','LineWidth',0.75,'DisplayName','Orig 3-sigma Cov');
            plot(sigmaP(1,:),sigmaP(2,:),'d','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean');

            plot(f_sample(1,:), f_sample(2,:),'b.','MarkerSize',1,'DisplayName','Transformed Sample');
            plot(f_mean(1),f_mean(2),'s','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','Transformed mean');
            plot(f_ellipse(1,:),f_ellipse(2,:),'b-','LineWidth',3,'DisplayName','Transformed 1-sigma Cov');
            plot(f3_ellipse(1,:),f3_ellipse(2,:),'b-','LineWidth',0.75,'DisplayName','Transformed 3-sigma Cov');

            plot(e_mean(1),e_mean(2),'s','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','UT est mean');
            plot(e_ellipse(1,:),e_ellipse(2,:),'m-','LineWidth',3,'DisplayName','UT est 1-sigma Cov');
            plot(e3_ellipse(1,:),e3_ellipse(2,:),'m-','LineWidth',0.75,'DisplayName','UT est 3-sigma Cov');
            plot(f_sigmaP(1,:),f_sigmaP(2,:),'d','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean');
            axis('equal');
            legend('location','best');
        end

        function testEstimate4D(obj, f, o_mean, o_cov)
            [e_mean, e_cov, sigmaP, f_sigmaP] = obj.estimate(f, o_mean, o_cov);
            o_cov_chol = chol(o_cov);
            oa_ellipse = obj.covEllipse(o_mean(1:2), o_cov_chol(1:2,1:2),1,1e4);
            oa3_ellipse = obj.covEllipse(o_mean(1:2), o_cov_chol(1:2,1:2),3,1e4);
            ob_ellipse = obj.covEllipse(o_mean(3:4), o_cov_chol(3:4,3:4),1,1e4);
            ob3_ellipse = obj.covEllipse(o_mean(3:4), o_cov_chol(3:4,3:4),3,1e4);
            Nsample=1e4;
            o_sample = obj.sampleGaussian(o_mean,o_cov_chol,Nsample);
            f_mean = f(o_mean);
            temp = f([oa_ellipse;ob_ellipse]);
            fa_ellipse = temp(1:2,:);
            fb_ellipse = temp(3:4,:);
            temp = f([oa3_ellipse;ob3_ellipse]);
            fa3_ellipse = temp(1:2,:);
            fb3_ellipse = temp(3:4,:);
            f_sample = f(o_sample);
            
            e_cov_chol = chol(e_cov);
            Pa = [1,0;0,1;0,0;0,0]';
            Pb = [0,0;0,0;1,0;0,1]';
            eCovA = Pa*e_cov_chol*Pa';
            eCovB = Pb*e_cov_chol*Pb';
            ea_ellipse = obj.covEllipse(e_mean(1:2), eCovA);
            ea3_ellipse = obj.covEllipse(e_mean(1:2), eCovA,3);
            eb_ellipse = obj.covEllipse(e_mean(3:4), eCovB);
            eb3_ellipse = obj.covEllipse(e_mean(3:4), eCovB,3);
            e_sample = obj.sampleGaussian(e_mean,e_cov_chol,Nsample);
                        
            figure();
            subplot(2,1,1);
            hold('on');
            plot(o_sample(1,:), o_sample(2,:),'r.','MarkerSize',1,'DisplayName','Orig Sample (A)');
            plot(o_mean(1),o_mean(2),'s','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean (A)');
            plot(oa_ellipse(1,:),oa_ellipse(2,:),'r-','LineWidth',3,'DisplayName','Orig 1-sigma Cov (A)');
            plot(oa3_ellipse(1,:),oa3_ellipse(2,:),'r-','LineWidth',0.75,'DisplayName','Orig 3-sigma Cov (A)');
            plot(sigmaP(1,:),sigmaP(2,:),'d','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean (A)');

            plot(o_sample(3,:), o_sample(4,:),'b.','MarkerSize',1,'DisplayName','Orig Sample (B)');
            plot(o_mean(3),o_mean(4),'s','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean (B)');
            plot(ob_ellipse(1,:),ob_ellipse(2,:),'b-','LineWidth',3,'DisplayName','Orig 1-sigma Cov (B)');
            plot(ob3_ellipse(1,:),ob3_ellipse(2,:),'b-','LineWidth',0.75,'DisplayName','Orig 3-sigma Cov (B)');
            plot(sigmaP(3,:),sigmaP(4,:),'d','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','Orig mean (B)');


            plot(f_sample(1,:), f_sample(2,:),'m.','MarkerSize',1,'DisplayName','Transformed Sample (A)');
            plot(f_mean(1),f_mean(2),'s','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','Transformed mean (A)');
            plot(fa_ellipse(1,:),fa_ellipse(2,:),'m-','LineWidth',3,'DisplayName','Transformed 1-sigma Cov (A)');
            plot(fa3_ellipse(1,:),fa3_ellipse(2,:),'m-','LineWidth',0.75,'DisplayName','Transformed 3-sigma Cov (A)');

            plot(f_sample(3,:), f_sample(4,:),'g.','MarkerSize',1,'DisplayName','Transformed Sample (B)');
            plot(f_mean(3),f_mean(4),'s','MarkerFaceColor',[0,1,0],'MarkerEdgeColor',[0,0,0],'DisplayName','Transformed mean (B)');
            plot(fb_ellipse(1,:),fb_ellipse(2,:),'g-','LineWidth',3,'DisplayName','Transformed 1-sigma Cov (B)');
            plot(fb3_ellipse(1,:),fb3_ellipse(2,:),'g-','LineWidth',0.75,'DisplayName','Transformed 3-sigma Cov (B)');

            plot(e_mean(1),e_mean(2),'s','MarkerFaceColor',[.5,.5,0],'MarkerEdgeColor',[0,0,0],'DisplayName','UT est mean (A)');
            plot(ea_ellipse(1,:),ea_ellipse(2,:),'-','Color',[.5,.5,0], 'LineWidth',3,'DisplayName','UT est 1-sigma Cov (A)');
            plot(ea3_ellipse(1,:),ea3_ellipse(2,:),'-','Color',[.5,.5,0],'LineWidth',0.75,'DisplayName','UT est 3-sigma Cov (A)');
            plot(e_sample(1,:),e_sample(2,:),'o','markersize',1,'MarkerFaceColor',[.5,.5,0],'MarkerEdgeColor','none','DisplayName','Estimated sample (A)');
            plot(f_sigmaP(1,:),f_sigmaP(2,:),'d','MarkerFaceColor',[.5,.5,0],'MarkerEdgeColor',[0,0,0],'DisplayName','Transformed sigma points (A)');

            plot(e_mean(3),e_mean(4),'s','MarkerFaceColor',[0,.5,.5],'MarkerEdgeColor',[0,0,0],'DisplayName','UT est mean (B)');
            plot(eb_ellipse(1,:),eb_ellipse(2,:),'-','Color',[0,.5,.5], 'LineWidth',3,'DisplayName','UT est 1-sigma Cov (B)');
            plot(eb3_ellipse(1,:),eb3_ellipse(2,:),'-','Color',[0,.5,.5],'LineWidth',0.75,'DisplayName','UT est 3-sigma Cov (B)');
            plot(e_sample(3,:),e_sample(4,:),'o','markersize',1,'MarkerFaceColor',[0,.5,.5],'MarkerEdgeColor','none','DisplayName','Estimated sample (B)');
            plot(f_sigmaP(3,:),f_sigmaP(4,:),'d','MarkerFaceColor',[0,.5,.5],'MarkerEdgeColor',[0,0,0],'DisplayName','Transformed sigma points (B)');
            axis('equal');
            legend('location','best');
            
            subplot(2,1,2);
            d = sqrt(sum((e_sample(1:2,:) -e_sample(3:4,:)).^2));
            plot(1:Nsample, d,'k-');
        end


        function cross_cov = estimateCrossModelCov(obj,y_sigmaPts, y_mean, x_mean, x_cov)
            Nx = size(y_sigmaPts,1);
            Nsigma = size(y_sigmaPts,2);
            [x_sigmaPts, ~, centerCovWeight, generalWeight] = obj.sigmaPoints(x_mean, x_cov);
            x_delta =  repmat(x_mean,1,Nsigma) - x_sigmaPts;
            y_delta = repmat(y_mean,1,Nsigma) - y_sigmaPts;
            
            cross_cov = y_delta*diag(w)*x_delta';
        end
    end
    methods (Static=true)
        function sample = sampleGaussian(x_mean, x_cov_chol, N)
            % [in]
            %   x_mean : size:[Nx,1]
            %   x_cov : size:[Nx,Nx]
            %   N : [optional] number of samples to return >=1 [default=1000]
            % [out]
            %  sample : size[Nx,N] each column is an Nx dimensional sample
            if nargin<3
                N=1000;
            end
            Nx = size(x_mean,1);
            sample = repmat(x_mean,1,N)+x_cov_chol*randn(Nx,N);
        end

        function ellipse = covEllipse(x_mean, x_cov_chol, sigmaFactor, Ncirc)
            % [in]
            %   x_mean : size:[Nx,1]
            %   x_cov : size:[Nx,Nx]
            %   sigmaFactor : [optional] scalar>0 number of standard deviations outward to draw ellipse [default=1]
            %   Npts : [optional] number of points to display in circle [default=100]
            % [out]
            %  sample : size[Nx,N] each column is an Nx dimensional sample
            if nargin<3
                sigmaFactor=1;
            end
            if nargin<4
                Ncirc = 100; % number of points in the circle
            end
            alpha = linspace(0,2*pi,Ncirc);
            circle = [cos(alpha); sin(alpha)];
            ellipse = repmat(x_mean(:),1,Ncirc) + sigmaFactor*x_cov_chol*circle;
        end

    end
end


