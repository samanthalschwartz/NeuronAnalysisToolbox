classdef PairInteractionMCMCModel < handle
    properties (Constant=true)
        LOG2PI = log(2*pi);

        MODELS={'DistanceMarkov','Diffusive'};
        DEFAULT_RHO = 0.050; %microns - bond distance
        DEFAULT_SIGMA_REG = 0.005; %microns - 
        DefaultParams=struct('DistanceMarkov', struct(... %Distance Markov model
                                        'D1',0.050,... %(um^2/s) monomer diffusion constant 
                                        'D2',0.050,... %(um^2/s) dimer diffusion constant 
                                        'rho',0.050,... %(um) dimer bond distance between fluorophores
                                        'kon',1e-2,... %(1/s) (on) F->B rate (Fudge)
                                        'koff',1e0...%(1/s) (off) B->F rate.  (True kinetic rate)
                                        ),...
                             'Diffusive', struct(... %ErbanChapman with Diffusion
                                        'DA',0.150,... %(um^2/s) particle A diffusion constant 
                                        'DB',0.150,... %(um^2/s) particle B diffusion constant 
                                        'DAB',0.035,... %(um^2/s) dimer diffusion constant 
                                        'Dphi',5e1,... %(rad^2/s) dimer angle diffusion rate
                                        'rho',0.050,... %(um) dimer bond distance between fluorophores
                                        'sigma_rho',0.010,... %(um) dimer bond distance standard deviation
                                        'gamma',1e2,... %(1/s) dimer bond distance relaxation rate
                                        'rho_bind',0.050,... %(um) binding distance.  Binding is possible when distance < rho_bind
                                        'rho_unbind',0.050,... %(um) unbinding distance.  
                                        'lambda_bind',1e2,... % [1/s] Rate of binding when within rho_bind distance
                                        'koff',1e0...%(1/s) (off) B->F rate.  (True kinetic rate)
                                        )...
                            );

        %The Prior structure type describes the prior knowledge about the initial state and positions for 
        % a prior experiment (cell)
        DefaultPrior = struct( 'state',[.5, .5], ...
                               'posShape','Gaussian', ... Options: ['Rectangular','Circular','Gaussian']
                               'posPoint',[0,0],... Gauss:center; Circ:center; Rect:LowerLeftCorner [left,bot]
                               'posExtent',[1,0;0,1] ... 'Gauss':Covariance, 'Circ':radius, 'Rect':[w,h]
                              );

        DefaultData=struct('experimentId',[],... % Inteteger index indicating from which "cell" this data is taken.  Priors on positions are per-cell.
                    'Nobs',[],... % number of observations
                    'obs',[],...   % size:[Nobs, 9] matrix of observed positions [t ax ay bx by sax say sbx sby]
                    'particleTrue',[]... % size:[Nobs, 5] matrix of true values for hidden variables [state ax ay bx by]
                    );

        DefaultExperiment = struct('name',[],...
                                   'dateTime',[],...
                                   'frameRate',[],...
                                   'simulation',false,... % boolean - true if this is simulated data
                                   'simulationInteraction',[],... % Simulation interaction type
                                   'simulationParams',[],...      % Simulation parameters used
                                   'registrationSigma',0.005,... %(um) The standard deviation of the registration error                                        
                                   'prior',[]... % A Strucuture of shape DefaultPrior
                                  );

        DefaultResults = struct('model',[],...
                                'params',[],...
                                'particles',[],... % [Nobs, 5, K] cube of sampled particles
                                'weights',[],...   % [K,1] vector of particle weights
                                'llhParticles',[],...
                                'llhTrue',[],... % log-likelihood of true particle (if availible)
                                'llhTruePrior',[],...
                                'llhTrueAll',[]...
                                );
    end

    properties
        Data; %Struct-array of the data (pairs) and any known true data from simulation
        Experiment; %Struct-array of the experiments "cells" from which the various data are taken
        Results; %Struct-array of results under a particular model and parameter set

        simDefaultParams;
    end

    properties (Dependent=true)
       Ndata; 
    end

    methods
        
        function obj = PairInteractionMCMCModel()
            obj.simDefaultParams = obj.DefaultParams.Diffusive;
        end

        % Forward declare this method in browseParticles.m
        browseParticles(obj, exper_idx, Nparticles, ProposalType);

        function resetAll(obj)
            obj.Data=[];
            obj.Experiment=[];
            obj.Results=[];
        end

        function sim = simulateTrajectories(obj,N,nSteps,Tmax,interaction_type,simParams)
            % [in]
            %  N - number of pairs to simulate
            %  nSteps - number of steps to simulate
            %  Tmax - maximum time (s)
            %  interaction_type - [optional] options  'Bound' [default='Fixed']
            %    {'Bound', ... % Always bound
            %     'Flyby', ... %Non-interacting nearby routes at random time (mimics expexted "free" pair)
            %     'FixedFlyby', ... %Non-interacting crossing routes at time T/2 exactly.
            %     'Fixed', ... % A single fixed duration
            %     'Exponential', ... %Random single binding period at random position in tracjetory
            %     'ReactionSim'}; % Brownian reaction simulation selected interacting pairs
            %  simParams - A parameter structure for the Diffusive model with which to simulate
            if nargin<6
                simParams = obj.simDefaultParams;
            end
            if nargin<5
                interaction_type = 'Fixed';
            end
            nSimSteps = 10*nSteps;
            simDeltaT = Tmax/nSimSteps;
            frameDeltaT = Tmax/nSteps;
            sim = DimerSimulator(nSimSteps,simDeltaT);
            sim.DA = simParams.DA;
            sim.DB = simParams.DB;
            sim.DAB = simParams.DAB;
            sim.dimerRigidLength = simParams.rho;
            sim.dimerLengthStddev = simParams.sigma_rho;
            sim.dimerRelaxationRate = simParams.gamma;
            sim.dimerOrientationDriftRate = simParams.Dphi;
            sim.interactionKoffRate = simParams.koff;
            sim.interactionFixedTime = 1/simParams.koff;
            sim.reactionSimBindRadius = simParams.rho_bind;
            sim.reactionSimUnbindRadius = simParams.rho_unbind;
            sim.reactionSimBindProb = simParams.lambda_bind;
            sim.interactionType = interaction_type;

            Ts = sim.simulateTrajectories(N);
            Os = sim.observeTrajectories(Ts,frameDeltaT);
            if isempty(obj.Experiment)
                exper.name = 'Simulation1';
            else
                exper.name = findUniqueName([obj.Experiment(:).name],'Simulation%i');
            end
            exper.frameRate=frameDeltaT;
            exper.registrationSigma=sim.overlayErrorStddev;
            exper.simulation = true;
            exper.simulationInteraction = interaction_type;
            exper.simulationParams = simParams;
            exper.dateTime = datetime();
            obj.loadSimulatedData(Os, exper);
        end

        function sim = simulateTrajectoriesSMC(obj,Nparticles,Nsteps,Tmax,params)
            % [in]
            if nargin<5
                params = obj.simDefaultParams;
                params.wmean=0.025;
                params.wsigma=0.005;
            end
            
            
            smc = PairInteractionSMC(params);
            smc.setPrior(obj.DefaultPrior);
            [Os,Ps,llhAll, simParams] = smc.simulate(Nparticles,Nsteps,Tmax,'Ancestral',params);
            
            if isempty(obj.Experiment)
                exper.name = 'SMCSimulation1';
            else
                exper.name = findUniqueName([obj.Experiment(:).name],'SMCSimulation%i');
            end
            exper.frameRate=simParams.dt;
            exper.registrationSigma=0;
            exper.simulation = true;
            exper.simulationInteraction = 'Ancestral';
            exper.simulationParams = simParams;
            exper.dateTime = datetime();
            Ts=cell(Nparticles,1);
            for n=1:Nparticles
                Ts{n} = [Os(:,:,n)', Ps(:,:,n)'];
            end
            obj.loadSimulatedData(Ts, exper);
        end

        function loadSimulatedData(obj, observedTrajs, experiment)
            % [in]
            %  observedTrajs - Cell array of trajectory data where each cell is a matrix with 14 columns:
            %       [t obs_ax obs_ay obs_bx obs_by Sax Say Sbx Sby trueState true_ax true_ay true_bx true_by]
            exp_id = obj.addExperiment(experiment);
            for k=1:numel(observedTrajs)
                obs=observedTrajs{k};
                data.experimentId = exp_id;
                data.Nobs = size(obs,1);
                data.obs = obs(:,1:9);
                data.particleTrue = obs(:,10:14);
                obj.addData(data);
            end
            obj.estimatePositionPrior(exp_id);
        end

        function data_idxs = getExperimentDataIdxs(obj, exper_id)
            data_idxs = find([obj.Data(:).experimentId]==exper_id);
        end

        function [trueLLH,freeLLH] = compareTrueLLH(obj, data_idx)
            if nargin<2
                Ds = obj.Data;
            else
                Ds=obj.Data(data_idx);
            end
            N = numel(Ds);
            cLLH = zeros(N,1);
            cParallelLLH = zeros(N,1);
            trueLLH = zeros(N,1);
            freeLLH = zeros(N,1);
%             smc = PairInteractionSMC();
            for n = 1:N
                D = Ds(n);
                exper = obj.Experiment(D.experimentId);
                trueLLH(n) = obj.computeSingleParticleLLH_Distance(D.obs, D.particleTrue, exper.simulationParams, exper.prior, true);
%                 smc.setParams(exper.simulationParams);
%                 cLLH(n) = smc.computeLLH(D.obs',D.particleTrue');
%                 fprintf('True particle [%i] log-likelihood: Matlab:%.9g C:%.9g\n',n,trueLLH,cLLH);
                free_particle = D.particleTrue;
                free_particle(:,1) = 0;
                freeLLH(n) = obj.computeSingleParticleLLH_Distance(D.obs, free_particle, exper.simulationParams, exper.prior, true);
            end
%             cParallelLLH = smc.computeLLH(cellmap(@transp,[obj.Data(:).obs]ell(obj.Data(:).particleTrue)));
        end

        function debugLLH(obj,data_idx)
            if nargin<2
                data_idx=1;
            end
            D = obj.Data(data_idx);
            E = obj.Experiment(D.experimentId);
            params = E.simulationParams;
            prior = E.prior;
            smc = PairInteractionSMC(params);
            smc.setPrior(prior);
            [~, llhRef] = obj.computeSingleParticleLLH_Distance(D.obs,D.particleTrue, params,prior, false);
            llh=llhRef';

            llhC = smc.computeLLH_debug(D.obs', D.particleTrue');
            llhC_sum = smc.computeLLH(D.obs', D.particleTrue');
            llhCdebug_sum = sum(llhC(:));
            llhMatlab = sum(llh(:));
            d = abs(llhC-llhRef');
            d(abs(d)<=sqrt(eps)) = 0;
%             assert(max(d(:))==0);
            Nobs = size(llh,2); % ncols
            figure();
            hold('on');
            plot(1:Nobs,llh(1,:),'k-','DisplayName','state (free=0/bound=1)');
            plot(1:Nobs,llh(2,:),'-','Color',[1  0   0],'DisplayName','y ax / v cx');
            plot(1:Nobs,llh(3,:),'-','Color',[.7 0  0],'DisplayName','y ay / v cy');
            plot(1:Nobs,llh(4,:),'-','Color',[.4 0  0],'DisplayName','y bx / v r + logAbsDet(V^-1)');
            plot(1:Nobs,llh(5,:),'-','Color',[.2  0  0],'DisplayName','y by / v phi');
            plot(1:Nobs,llh(6,:),'-','Color',[0 .5   1],'DisplayName','x ax');
            plot(1:Nobs,llh(7,:),'-','Color',[0 .5  .7],'DisplayName','x ay');
            plot(1:Nobs,llh(8,:),'-','Color',[0 .8  .4],'DisplayName','x_bx');
            plot(1:Nobs,llh(9,:),'-','Color',[0 .5  .2],'DisplayName','x_by');
            plot(1:Nobs,llhC(1,:),'k--','DisplayName','[C] state (free=0/bound=1)');
            plot(1:Nobs,llhC(2,:),'--','Color',[1  0   0],'DisplayName','[C] y ax / v cx');
            plot(1:Nobs,llhC(3,:),'--','Color',[.7 0  0],'DisplayName','[C] y ay / v cy');
            plot(1:Nobs,llhC(4,:),'--','Color',[.4 0  0],'DisplayName','[C] y bx / v r + logAbsDet(V^-1)');
            plot(1:Nobs,llhC(5,:),'--','Color',[.2  0  0],'DisplayName','[C] y by / v phi');
            plot(1:Nobs,llhC(6,:),'--','Color',[0 .5   1],'DisplayName','[C] x ax');
            plot(1:Nobs,llhC(7,:),'--','Color',[0 .5  .7],'DisplayName','[C] x ay');
            plot(1:Nobs,llhC(8,:),'--','Color',[0 .8  .4],'DisplayName','[C] x_bx');
            plot(1:Nobs,llhC(9,:),'--','Color',[0 .5  .2],'DisplayName','[C] x_by');
            xlabel('observation idx');
            ylabel('llh');
            legend('location','best');
        end
    
        function [smc,D] = makeExperimentSMC(obj, exper_idx)
            E = obj.Experiment(exper_idx);
            prior = E.prior;
            params = E.simulationParams;
            data_idx = [obj.Data(:).experimentId]==exper_idx;
            D = obj.Data(data_idx);
            
            smc = PairInteractionSMC(params);
            smc.setPrior(prior);
            data = cellmap(@transpose, {D(:).obs});
            smc.addData(data);
        end


        function [Pmean, Pvar] = estimateMCMC(obj, data_idx, Nsamples, Nparticles, proposal_method)
            if nargin<5
                proposal_method='Transition';
            end
            if nargin<4
                Nparticles=200;
            end
            D = obj.Data(data_idx);
            E = obj.Experiment(D.experimentId);
            params = E.simulationParams;
            prior = E.prior;
            llh = zeros(Nsamples,1);
            Nobs = size(D.obs,1);
            P = zeros(5,Nobs,Nsamples);
            smc = PairInteractionSMC(params);
            smc.setPrior(prior);
            smc.addData(D.obs');
        
            smc.runParticleFilter(Nparticles, proposal_method);
            [newP, llh(1)] = smc.sampleParticle();
            P(:,:,1) = newP{1};
            llh_true = smc.computeLLH(D.particleTrue');
            fprintf('OptimalLLH: %.6g\n',llh_true);
            naccept=0;
            nreject=0;
            for n=2:Nsamples
                smc.runParticleFilter(Nparticles, proposal_method);
                [newP, new_llh] = smc.sampleParticle();
                Z = exp(new_llh-llh(n-1));  %Acceptance probability
                if rand() < Z %accept;
                    llh(n) = new_llh;
                    P(:,:,n) = newP{1};
                    naccept = naccept+1;
                else %reject
                    llh(n) = llh(n-1);
                    P(:,:,n) = P(:,:,n-1);
                    nreject = nreject+1;
                end
                accept_ratio = naccept/n*100;
                fprintf('Step:[%i/%i] #acc:%i #rej:%i acc%%:%.3f%% AccProb:%.3f old_llh:%.8g new_llh:%.8g \n',...
                        n,Nsamples,naccept,nreject,accept_ratio, min(1,Z), sum(llh(n-1,:)),sum(new_llh));
            end
            Pmean = mean(P,3);
            PT = D.particleTrue;
            ts = D.obs(:,1);
            dt = diff(ts);
            figure();
            subplot(2,2,1);
            hold();
            title('Binding State');
            vs_mean = obj.transV(Pmean(2:5,:)');
            vs_true = obj.transV(PT(:,2:5));
            plot(ts,vs_true(:,3),'b--','DisplayName','True distance');
            plot(ts,vs_mean(:,3),'b-','DisplayName','Estimated distance');
            yl=ylim();
            stairs([0;ts(1:end-1)], 0.99*yl(2)*PT(:,1), 'k-','DisplayName','True state estimate');
            plot(ts, 0.99*yl(2)*Pmean(1,:), 'r-','DisplayName','State estimate');
            for t=2:Nobs
                if PT(t,1)
                    xs=ts([t-1,t,t,t-1]);
                    ys=[0,0, yl(2), yl(2)];
                    h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                    h.Annotation.LegendInformation.IconDisplayStyle='off';
                end
            end
            legend('location','best');

            bound= logical(PT(:,1));
            subplot(2,2,2);
            title('Error in cx,cy,r when Z=Bound');
            hold();
            plot(ts(bound), abs(vs_mean(bound,1)-vs_true(bound,1)), 'r-', 'DisplayName','cx error');
            plot(ts(bound), abs(vs_mean(bound,2)-vs_true(bound,2)), 'b-', 'DisplayName','cy error');
            plot(ts(bound), abs(vs_mean(bound,3)-vs_true(bound,3)), 'k-', 'DisplayName','r error');
            legend('location','best');

            subplot(2,2,3);
            title('Error in phi when Z=Bound');
            hold();
            plot(ts(bound), mod(abs(vs_mean(bound,4)-vs_true(bound,4)),pi), 'k-', 'DisplayName','phi error');
            legend('location','best');

            subplot(2,2,4);
            hold();
            title('Particle position error');
            plot(ts,abs(PT(:,2)-Pmean(2,:)'),'r-','DisplayName','ax error');
            plot(ts,abs(PT(:,3)-Pmean(3,:)'),'g-','DisplayName','ay error');
            plot(ts,abs(PT(:,4)-Pmean(4,:)'),'b-','DisplayName','bx error');
            plot(ts,abs(PT(:,5)-Pmean(5,:)'),'m-','DisplayName','by error');
            plot(ts,mean(D.obs(:,6:9),2),'k-','DisplayName','mean obs std');
            plot(ts,3*mean(D.obs(:,6:9),2),'k--','DisplayName','mean obs 3 x std');
            legend('location','best');
            yl=ylim();
            for t=2:Nobs
                if PT(t,1)
                    xs=ts([t-1,t,t,t-1]);
                    ys=[yl(1),yl(1), yl(2), yl(2)];
                    h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                    h.Annotation.LegendInformation.IconDisplayStyle='off';
                end
            end
        end


        function [P,w_norm,llhP,pobs] = runParticleFilter(obj, data_idx, K, plotOn)
            if nargin<2
                Ds = obj.Data;
            else
                Ds=obj.Data(data_idx);
            end
            if nargin <4
                plotOn=false;
            end
            N = numel(Ds);
            P = cell(N,1);
            w_norm = cell(N,1);
            llhP = cell(N,1);
            pobs = cell(N,1);
            smc = PairInteractionSMC();
            for n=1:N
                D=Ds(n);
                exper = obj.Experiment(D.experimentId);
                [P{n}, w_norm{n}, llhP{n}, pobs{n}] = obj.condensationParticleFilter(D.obs, K, exper.simulationParams, exper.prior, plotOn);
%                 smc.setParams(exper.simulationParams);
%                 smc.setPrior(exper.prior);
%                 ps = smc.getParams();
%                 smc.clearData();
%                 smc.addData({D.obs'});
%                 smc.runParticleFilter(100);
%                 [particles,llh] = smc.sampleParticle();
            end
            if N==1
                P=P{1};
                w_norm=w_norm{1};
                llhP=llhP{1};
                pobs=pobs{1};
            end
        end

        function [P,w_norm,llhP,llh_obs] = runStateParticleFilter(obj, data_idx, K, params, prior)
            if nargin<2
                Ds = obj.Data;
            else
                Ds=obj.Data(data_idx);
            end
            N = numel(Ds);
            P = cell(N,1);
            w_norm = cell(N,1);
            llhP = cell(N,1);
            llh_obs = zeros(N,1);
            for n=1:N
                D=Ds(n);
                if nargin<5
                    exper = obj.Experiment(D.experimentId);
                    if nargin<4
                        params = exper.simulationParams;
                    end
                    prior = exper.prior;
                end
                ys = D.particleTrue(:,2:5);
                ts = D.obs(:,1);
                [P{n}, w_norm{n}, llhP{n}, llh_obs(n)] = obj.condensationStatesParticleFilter(ts,ys, K, params, prior);
            end
            if N==1
                P=P{1};
                w_norm=w_norm{1};
                llhP=llhP{1};
            end
        end

        function [ys_mean, ys_cov, llh_xs] = runKalmanFilter(obj, data_idx, state, plotOn)
            if nargin<2
                Ds = obj.Data;
            else
                Ds=obj.Data(data_idx);
            end
            N = numel(Ds);
            if state % bound
                vs_mean = cell(N,1);
                vs_cov = cell(N,1);
            else  
                ys_mean = cell(N,1);
                ys_cov = cell(N,1);
            end
            llh_xs = zeros(N,1);
            for n=1:N
                D=Ds(n);
                exper = obj.Experiment(D.experimentId);
                params = exper.simulationParams;
                prior = exper.prior;
                if state % Bound
                    [vs_mean{n}, vs_cov{n}, llh_xs(n)] = obj.computeUnscentedKalmanFilter_Bound(D.obs, params, prior);
                    if nargin>3 && plotOn
                        figure();
                        subplot(2,2,1);
                        hold('on');
                        ts = D.obs(:,1);
                        vs_true = PaitInteractionMCMCModel.transV(D.particleTrue(:,2:5)); 
                        vs_est = vs_mean{n}';
                        plot(ts,vs_true(:,1),'r:','DisplayName','cx true');
                        plot(ts,vs_true(:,2),'g:','DisplayName','cy true');
                        plot(ts,vs_true(:,3),'b:','DisplayName','r true');
                        plot(ts,vs_true(:,4),'m:','DisplayName','phi true');
                        plot(ts,vs_est(:,1),'r-','DisplayName','cx est');
                        plot(ts,vs_est(:,2),'g-','DisplayName','cy est');
                        plot(ts,vs_est(:,3),'b-','DisplayName','r est');
                        plot(ts,vs_est(:,4),'m-','DisplayName','phi est');
                        legend('location','best');
                        xlabel('time');
                        ylabel('position');

                        subplot(2,2,2);
                        rmse = sqrt(sum((vs_true-vs_est).^2,2));
                        hold('on');
                        plot(ts,rmse,'r-','DisplayName','UKF RMSE');
                        legend('location','best');
                        xlabel('time');
                        ylabel('position error');
                    end
                else%free
                    [ys_mean{n}, ys_cov{n}, llh_xs(n)] = obj.computeKalmanFilter_Free(D.obs, params, prior);
                    [ys_sm_mean, ys_sm_cov] = obj.computeKalmanSmoother_Free(D.obs, params, prior);
                    if nargin>3 && plotOn
                        figure();
                        subplot(2,2,1);
                        hold('on');
                        ts = D.obs(:,1);
                        ys_true = D.particleTrue(:,2:5); 
                        ys_sm = ys_sm_mean';
                        ys_est = ys_mean{n}';
                        plot(ts,ys_true(:,1),'r:','DisplayName','ax true');
                        plot(ts,ys_true(:,2),'g:','DisplayName','ay true');
                        plot(ts,ys_true(:,3),'b:','DisplayName','bx true');
                        plot(ts,ys_true(:,4),'m:','DisplayName','by true');
                        plot(ts,ys_sm(:,1),'r--','DisplayName','ax sm');
                        plot(ts,ys_sm(:,2),'g--','DisplayName','ay sm');
                        plot(ts,ys_sm(:,3),'b--','DisplayName','bx sm');
                        plot(ts,ys_sm(:,4),'m--','DisplayName','by sm');
                        plot(ts,ys_est(:,1),'r-','DisplayName','ax est');
                        plot(ts,ys_est(:,2),'g-','DisplayName','ay est');
                        plot(ts,ys_est(:,3),'b-','DisplayName','bx est');
                        plot(ts,ys_est(:,4),'m-','DisplayName','by est');
                        legend('location','best');
                        xlabel('time');
                        ylabel('position');

                        subplot(2,2,2);
                        res = sqrt(sum((ys_true-ys_est).^2,2));
                        sm_res = sqrt(sum((ys_true-ys_sm).^2,2));
                        hold('on');
                        plot(ts,res,'r-','DisplayName','Filter RMSE');
                        plot(ts,sm_res,'b-','DisplayName','Smother RMSE');
                        legend('location','best');
                        xlabel('time');
                        ylabel('position error');
                    end
                end
                
            end
            if N==1
                ys_mean = ys_mean{1};
                ys_cov = ys_cov{1};
            end

        end
        


        function koff = statesKoffMCMC(obj, exper_idx, Nsamples, Nparticles, proposal_method)
            % use the states particle filter to estimate the koff from known positions ys
            % N number of generated samples
            if nargin<5
                proposal_method='Transition';
            end
            if nargin<4
                Nparticles=200;
            end
            E = obj.Experiment(exper_idx);
            prior = E.prior;
            params = E.simulationParams;
            data_idx = [obj.Data(:).experimentId]==exper_idx;
            D = obj.Data(data_idx);
            NData = numel(D);
            koff_sigma_factor = 0.015;
            koff = zeros(Nsamples,1);
            llh = zeros(Nsamples,NData);
            true_koff = params.koff;
            koff(1) = true_koff;
%             [~,w_norm,llhP] =  obj.runStateParticleFilter(data_idx, K, params, prior);
% %             llh(1,:) = cellfun(@(w,L) L(randsample(K,1,true,w)), w_norm, llhP);
%             for m=1:ND
%                 [~,w_norm,llhP,llhObs] = obj.condensationParticleFilter(D(m).obs, K, params, prior,false);
%                 idx = randsample(numel(w_norm),1,true,w_norm);  %random weighted sample of N from N with replacement!                
%                 llh(1,m)=llhP(idx);
%             end
            smc = PairInteractionSMC(params);
            smc.setPrior(prior);
            data = cellmap(@transpose, {D(:).obs});
            smc.addData(data);
            smc.runParticleFilter(Nparticles,proposal_method);
            [ps,llh(1,:)] = smc.sampleParticle();
            llh_true = smc.computeLLH(cellmap(@transpose, {D(:).particleTrue}));
            fprintf('OptimalLLH: %.6g\n',sum(llh_true));
%             llhRef = smc.computeLLH(ps);
%             fprintf('Maximum LLH discrepency: %.6g\n',max(abs(llh(1,:)-llhRef')));
            naccept=0;
            nreject=0;
            for n=2:Nsamples
                sigma_koff = koff(n-1)*koff_sigma_factor;
                new_koff = koff(n-1) + randn(1,1)*sigma_koff;
                params.koff=new_koff;
                smc.runParticleFilter(Nparticles,proposal_method);
                [ps,new_llh] = smc.sampleParticle();
%                 llhRef = smc.computeLLH(ps);
%                 fprintf('Maximum LLH discrepency: %.6g',max(abs(new_llh-llhRef)));

                q_log_ratio = obj.gaussianLLH(new_koff-koff(n-1), new_koff*koff_sigma_factor) - ...
                              obj.gaussianLLH(new_koff-koff(n-1), koff(n-1)*koff_sigma_factor);
                p_log_ratio = sum(new_llh) - sum(llh(n-1,:));
                Z = exp(p_log_ratio + q_log_ratio);
                if rand() < Z %accept
                    koff(n) = new_koff;
                    llh(n,:) = new_llh;
                    naccept = naccept+1;
                else %reject
                    koff(n) = koff(n-1);
                    llh(n,:) = llh(n-1,:);
                    nreject = nreject+1;
                end
                accept_ratio = naccept/n*100;
                fprintf('Step:[%i/%i] #acc:%i #rej:%i acc%%:%.3f%% meankoff:%.4g koff:%.4g new_koff:%.4g AccProb:%.3f old_llh:%.8g new_llh:%.8g \n',...
                        n,Nsamples,naccept,nreject,accept_ratio, mean(koff(1:n)), koff(n),new_koff,min(1,Z),...
                        sum(llh(n-1,:)),sum(new_llh));
            end
            mkoff = mean(koff);
            figure();
            histogram(koff,'DisplayName','MCMC est of p(koff|data)');
            hold();
            yl=ylim();
            plot(true_koff*[1,1],yl,'r:','DisplayName','True Koff'); 
            plot(mkoff*[1,1],yl,'r-','DisplayName','MCMC esitmated Koff'); 
            xlabel('koff [1/sec]');
            ylabel('probability');
            legend('location','best');
        end

        function plotParticles(obj, data_idx,K)
            [P,w_norm,llhP,pobs] = obj.runParticleFilter(data_idx, K);
            D=obj.Data(data_idx);
            ts = D.obs(:,1);
            figure();
            ax=subplot(3,2,1);
            hold(ax);
            stairs(ts, D.particleTrue(:,1),'r-');
            lag = min(diff(ts))/10;
%             for k=1:K
%                 stairs(ts+lag,P(:,1,k),'k-');
%             end
            meanZs = squeeze(P(:,1,:))*w_norm';
            stairs(ts+lag,meanZs,'b-');
            subplot(3,2,2);
            hold('on');
            plot(ts,D.particleTrue(:,2),'r-','DisplayName','ax');
            plot(ts,D.particleTrue(:,3),'g-','DisplayName','ay');
            plot(ts,D.particleTrue(:,4),'b-','DisplayName','bx');
            plot(ts,D.particleTrue(:,5),'m-','DisplayName','by');
            meanAx = squeeze(P(:,2,:))*w_norm';
            meanAy = squeeze(P(:,3,:))*w_norm';
            meanBx = squeeze(P(:,4,:))*w_norm';
            meanBy = squeeze(P(:,5,:))*w_norm';
            plot(ts,meanAx,'r:','DisplayName','est_ax');
            plot(ts,meanAy,'g:','DisplayName','est_ay');
            plot(ts,meanBx,'b:','DisplayName','est_bx');
            plot(ts,meanBy,'m:','DisplayName','est_by');
            xlabel('time (t)');
            ylabel('position (um)');
            legend('location','best');
        end
    end
       

    methods (Static=true, Access=public)
        %% Prior position distribution helpers
        function llh = llh_R_prior(rs, prior)
            % Prior probability of seperation rs
            % [in]
            %  rs - size:[N,1] distance samples
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  llh - size:[N,1] log-likelihood of rs(n)
            ext =  prior.posExtent;
            switch lower(prior.posShape)
                case {'gaussian','gauss'}
                    llh = log(PairInteractionMCMCModel.pairR_PDF_2Dgauss(rs,ext));
                case {'circular','circ'}
                    llh = log(PairInteractionMCMCModel.pairR_PDF_circ(rs, ext));
                case {'rectangular','rect'}
                    llh = log(PairInteractionMCMCModel.pairR_PDF_2Dgauss(rs,ext(1),ext(2)));                
            end
        end

        function llh = llh_pos_prior(pos, prior)
            % Sample the position prior
            % [in]
            %  pos - size:[N,2] position samples
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  llh - size:[N,1] total (sum) log-likelihood of pos(n)
            ext =  prior.posExtent;
            pt = prior.posPoint;
            N = size(pos,1);
            switch lower(prior.posShape)
                case {'gaussian','gauss'}
                    llh = PairInteractionMCMCModel.multiGaussianLLH((pos-repmat(pt,N,1))',ext);
                case {'circular','circ'}
                    llh = repmat(-log(pi*ext^2),N,1);
                case {'rectangular','rect'}
                    llh = repmat(-sum(log(ext)),N,1);
            end
        end

        function [m,V] = pos_prior_distribution_free(prior)
            % Get the mean and covariance for a free pair of particles given a gaussian prior
            % [in]
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  m - size:[4,1] mean positions rows:[ax ay bx by]
            %  V - size:[4,4] covariance
            ext =  prior.posExtent;
            pt = prior.posPoint;
            switch lower(prior.posShape)
                case {'gaussian','gauss'}
                    m = [pt(:); pt(:)];
                    V = blkdiag(ext,ext);
                otherwise
                    error('PairInteractionMCMCModel:pos_prior_distribution_free','Non gaussian prior');
            end
        end

        function [v0_mean,v0_cov] = pos_prior_distribution_bound(params,prior)
            % Get the mean and covariance for a bound pair of particles
            % in v-space
            %
            % For theta:
            % Use mean of 0 and variance of (2*pi)^2 to simulate uniform distribution over [0,2*pi) when taken mod(2*pi) 
            % we'll use vonMises to compute the pdf, but for a gaussian representation, 
            % this will be just as good.
            % [in]
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  v0_mean - size:[4,1] mean positions rows:[cx cy r theta]
            %  v0_cov - size:[4,4] covariance
            ext =  prior.posExtent;
            pt = prior.posPoint;
            switch lower(prior.posShape)
                case {'gaussian','gauss'}
                    v0_mean = [pt(:); params.rho; 0];
                    v0_cov = blkdiag(ext,params.sigma_rho^2,(2*pi)^2);
                otherwise
                    error('PairInteractionMCMCModel:pos_prior_distribution_free','Non gaussian prior');
            end
        end


        function pos = samplePosPrior(N, prior)
            % Sample the position prior
            % [in]
            %  N - numer of samples
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  pos - size:[N,2] position samples
            ext =  prior.posExtent;
            pt = prior.posPoint;
            switch lower(prior.posShape)
                case {'gaussian','gauss'}
                    pos = repmat(prior.posPoint,N,1) + randn(N,2)*chol(ext);
                case {'circular','circ'}
                    r = sqrt(ext^2*rand(N,1));
                    theta = 2*pi*rand(N,1);
                    pos = [pt(1)+r.*cos(theta); pt(2)+r.*sin(theta)];
                case {'rectangular','rect'}
                    pos = rand(N,2).*repmat(ext,N,1)+repmat(pt,N,1);
            end
        end

        %% Prior sampling distribution q()
        


        function [z0,y0,llh,llh_all] = sample_Q0(x0,w0,N,params,prior)
            % Sample the prior proposal distribution for state given true positions
            %   (z0,y0) ~ q(z0, y0 | x0 w0) = q(z0 | y0) q(y0 | x0 w0).
            % Also report the log-likelihood of the sample.
            % [in]
            %  x0 - size:[N,4] observations at time t
            %  w0 - size:[N,4] observations sigmas at time t
            %  N - number of samples
            %  params - params struct
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  y0 - size:[N,4] initial y(0) positions
            %  llh - size:[N,1] total (sum) log-likelihood of (z0, y0) 
            %  llh_all - [optional] size:[N,5] log-likelihood of (z0, y0) by component
            %                       columns: [z ax ay bx by]
            if size(x0,1)==1 && N>1
                x0 = repmat(x0(:)',N,1);
            end
            if size(w0,1)==1 && N>1
                w0 = repmat(w0(:)',N,1);
            end
            llh_all = zeros(N,5);
            y0 = x0 + randn(N,4).*w0;
            llh_all(:,2:5) = PairInteractionMCMCModel.gaussianLLH(x0-y0,w0.^2);            
            [z0,pz0] = PairInteractionMCMCModel.sample_Q0_Z(y0,params,prior);
            llh_all(:,1) = log(pz0);
            llh = sum(llh_all,2);
        end

X
        function [z0,pz0] = sample_Q0_Z(y0,params,prior)
            % Sample the prior proposal distribution for Z given Y
            %   (z0) ~ q(z0 | y0).
            % Also report the log-likelihood of the sample.
            % [in]
            %  y0 - size:[N,4] true positions [ax ay bx by] at time t
            %  params - params struct
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  z0 - size:[1,N] initial z(0)states 0=free 1=bound
            %  pz0 - [optional] size:[1,N] probability of each z0
            N = size(y0,1);
            llh_free =  PairInteractionMCMCModel.llhG0_Y(y0,zeros(N,1),params,prior)';
            llh_bound =  PairInteractionMCMCModel.llhG0_Y(y0,ones(N,1),params,prior)';
            llh = [llh_free, llh_bound]; % cols are: [free, bound];
            p_y0_given_z0 = exp(llh); % p(y0 | z0)
            p_y0_and_z0 = p_y0_given_z0 .* repmat(prior.state,N,1);  % p(y0 z0) = p(y0|z0)p(z0)
            % Bayes rule p(z0 | y0) = p(y0 | z0) p(z0) / p(y0), where p(y0)=sum_x p(y0|z0)p(z0)
            p_z0 = p_y0_and_z0 ./ repmat(sum(p_y0_and_z0,2),1,2); 
            z0 = rand(N,1) > p_z0(:,1);
            pz0 = p_z0(:,1); %set the output probability to the prob for the free state
            pz0(z0) = p_z0(z0,2); % then fix the bound probability for those that are not free
            z0=z0';
            pz0=pz0'; %return as a row vector
        end

        function [zt,pzt] = sample_Q_Z(yt,z0,y0,t,params,RD)
            % This is the proposal distribution for state given true positions and previous positions and state
            % Sample zt ~ q(zt | yt, z0, y0, theta) 
            % for a set of independent particles with known yt,z0,y0 and given t, and params theta
            % [in]
            %  yt - size:[N,4] true positions at time t
            %  z0 - size:[N] initial states 0=free 1=bound
            %  y0 - size:[N,4] initial positions
            %  t - scalar or vector of times.  if vector must be size:[N] and z(i) is sampled at time t(i)
            %  params - params struct
            %  RD - RDCapture object 
            % [out]
            %  zt - size:[N,1] sampled z(t) values according to data-dependent sampling propogator q().
            %  pzt - size:[N,1] probability of each zt
            N=numel(z0);
            llh_yt_free = PairInteractionMCMCModel.llhG_Y(yt,zeros(N,1),t,z0,y0,params); %log(p(yt|zt=F,z0,y0))
            llh_yt_bound = PairInteractionMCMCModel.llhG_Y(yt,ones(N,1),t,z0,y0,params); %log(p(yt|zt=B,z0,y0))
            p_zt_free = PairInteractionMCMCModel.pG_Z(zeros(N,1),t,z0,y0,params,RD); %P(zt=F) = 1-P(zt=B)
            p_zt_and_yt = [exp(llh_yt_free).*p_zt_free, exp(llh_yt_bound).*(1-p_zt_free)];
            % Bayes rule p(zt | yt z0 y0) = p(yt | zt z0 y0) p(zt | z0 y0) / p(yt | z0 y0), 
            %   where p(yt | z0 y0)=sum_zt p(yt | zt z0 y0) p(zt | z0 y0).
            p_zt_given_yt = p_zt_and_yt ./ repmat(sum(p_zt_and_yt,2),1,2); 
            zt = rand(N,1) > p_zt_given_yt(:,1); %is bound if greater than the free prob.
            pzt = p_zt_given_yt(:,1);
            pzt(zt) = p_zt_given_yt(zt,2);
        end


        function [zt,pzt] = sample_Q_Z_X(xt,wt,z0,y0,t,params,RD)
            % This is the forward state sampling distribution for state given true positions and previous positions and state
            % Sample zt ~ q(zt | xt, wt, z0, y0, theta) 
            % for a set of independent particles with known yt,z0,y0 and given t, and params theta
            % [in]
            %  xt - size:[N,4] observed positions at time t
            %  wt - size:[N,4] observed position variances at time t
            %  z0 - size:[N] initial states 0=free 1=bound
            %  y0 - size:[N,4] initial positions
            %  t - scalar or vector of times.  if vector must be size:[N] and z(i) is sampled at time t(i)
            %  params - params struct
            %   RD - RDCapture object 
            % [out]
            %  zt - size:[N,1] sampled z(t) values according to data-dependent sampling propogator q().
            %  pzt - size:[N,1] probability of each zt
            N=numel(z0);
            if numel(t)==1
                t = t*ones(N,1);
            end

            %log(p(xt|wt,zt=F,z0,y0))
            Q = 2*[params.DA, params.DA, params.DB, params.DB];
            xt_cov_free = repmat(Q,N,1).*repmat(t,1,4)+wt.^2;
            llh_xt_free = sum(PairInteractionMCMCModel.gaussianLLH(xt-y0,xt_cov_free),2);
            %log(p(xt|wt,zt=B,z0,y0))
            ut=UnscentedTransform(0,.5,2);
            v0 = PairInteractionMCMCModel.transV(y0);
            vxt = PairInteractionMCMCModel.transV(xt);
            vmean = v0;
            E=exp(-params.gamma*t);
            vmean(:,3) = E.*vmean(:,3) + (1-E).*params.rho;
            vvar = [2*params.DAB*t, 2*params.DAB*t, params.sigma_rho^2*(1-E.^2), min(2*pi^2,2*params.Dphi*t)];
            func = @(vs) PairInteractionMCMCModel.transVinv(vs')';
            llh_xt_bound = zeros(N,1);
            
            plotOn=false;
            if plotOn
                figure();
            end
            for n=1:N
                dmin = params.rho+5*params.sigma_rho+3*2*t(n)*(params.DA+params.DB);
                if vxt(3)<dmin
                    [est_y_mean,est_y_cov] = ut.estimate(func, vmean(n,:)', diag(vvar(n,:)));
%                     ymean = PairInteractionMCMCModel.transVinv(vmean(n,:));
                    % Use UT estimate to infer xt distribution
                    est_xt_mean = est_y_mean; 
                    est_xt_cov = est_y_cov+diag(wt(n,:).^2);
                    [pred_yt_mean, pred_yt_cov] = PairInteractionMCMCModel.gaussianBayes(est_y_mean,est_y_cov, xt(n,:), wt(n,:).^2);
                    llh_xt_bound(n) = PairInteractionMCMCModel.multiGaussianLLH(est_xt_mean-xt(n,:)',est_xt_cov);
                    if plotOn
                        subplot(1,1,1);
                        cla();
                        hold('on');
                        Npts=1000;
                        vts = repmat(vmean(n,:),Npts,1) + randn(Npts,4)*chol(diag(vvar(n,:)));
                        yvts = PairInteractionMCMCModel.transVinv(vts);
                        
                        plot(yvts(:,1),yvts(:,2),'.','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'DisplayName','ya(t) sample drawn from v(t) dist');
                        plot(yvts(:,3),yvts(:,4),'.','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'DisplayName','yb(t) sample drawn from v(t) dist');
                        a_cov_chol = diag(wt(n,1:2));
                        a1_ellipse = PairInteractionMCMCModel.covEllipse(xt(n,1:2), a_cov_chol, 1);
                        a3_ellipse = PairInteractionMCMCModel.covEllipse(xt(n,1:2), a_cov_chol, 3);
                        b_cov_chol = diag(wt(n,3:4));
                        b1_ellipse = PairInteractionMCMCModel.covEllipse(xt(n,3:4), b_cov_chol, 1);
                        b3_ellipse = PairInteractionMCMCModel.covEllipse(xt(n,3:4), b_cov_chol, 3);

                        est_y_cov_chol = chol(est_y_cov);
                        est_a1_ellipse = PairInteractionMCMCModel.covEllipse(est_y_mean(1:2), est_y_cov_chol(1:2,1:2), 1);
                        est_a3_ellipse = PairInteractionMCMCModel.covEllipse(est_y_mean(1:2), est_y_cov_chol(1:2,1:2), 3);
                        est_b1_ellipse = PairInteractionMCMCModel.covEllipse(est_y_mean(3:4), est_y_cov_chol(3:4,3:4), 1);
                        est_b3_ellipse = PairInteractionMCMCModel.covEllipse(est_y_mean(3:4), est_y_cov_chol(3:4,3:4), 3);

                        
                        pred_yt_cov_chol = chol(pred_yt_cov);
                        pred_a1_ellipse = PairInteractionMCMCModel.covEllipse(pred_yt_mean(1:2), pred_yt_cov_chol(1:2,1:2), 1);
                        pred_a3_ellipse = PairInteractionMCMCModel.covEllipse(pred_yt_mean(1:2), pred_yt_cov_chol(1:2,1:2), 3);
                        pred_b1_ellipse = PairInteractionMCMCModel.covEllipse(pred_yt_mean(3:4), pred_yt_cov_chol(3:4,3:4), 1);
                        pred_b3_ellipse = PairInteractionMCMCModel.covEllipse(pred_yt_mean(3:4), pred_yt_cov_chol(3:4,3:4), 3);

                        est_xt_cov_chol = chol(est_xt_cov);
                        est_axt1_ellipse = PairInteractionMCMCModel.covEllipse(est_xt_mean(1:2), est_xt_cov_chol(1:2,1:2), 1);
                        est_axt3_ellipse = PairInteractionMCMCModel.covEllipse(est_xt_mean(1:2), est_xt_cov_chol(1:2,1:2), 3);
                        est_xbt1_ellipse = PairInteractionMCMCModel.covEllipse(est_xt_mean(3:4), est_xt_cov_chol(3:4,3:4), 1);
                        est_xbt3_ellipse = PairInteractionMCMCModel.covEllipse(est_xt_mean(3:4), est_xt_cov_chol(3:4,3:4), 3);

                        free_xt_cov_chol = diag(sqrt(xt_cov_free(n,:)));
                        free_axt1_ellipse = PairInteractionMCMCModel.covEllipse(y0(n,1:2), free_xt_cov_chol(1:2,1:2), 1);
                        free_axt3_ellipse = PairInteractionMCMCModel.covEllipse(y0(n,1:2), free_xt_cov_chol(1:2,1:2), 3);
                        free_xbt1_ellipse = PairInteractionMCMCModel.covEllipse(y0(n,3:4), free_xt_cov_chol(3:4,3:4), 1);
                        free_xbt3_ellipse = PairInteractionMCMCModel.covEllipse(y0(n,3:4), free_xt_cov_chol(3:4,3:4), 3);


                        plot(est_a1_ellipse(1,:),est_a1_ellipse(2,:),'--','Color',[.5,0,0],'LineWidth',3,'DisplayName','est ya(t) cov 1-sigma');
%                         plot(est_a3_ellipse(1,:),est_a3_ellipse(2,:),'--','Color',[.5,0,0],'LineWidth',1,'DisplayName','est obs (a) cov 3-sigma');
                        plot(est_b1_ellipse(1,:),est_b1_ellipse(2,:),'--','Color',[0,0,.5],'LineWidth',3,'DisplayName','est obs (b) cov 1-sigma');
%                         plot(est_b3_ellipse(1,:),est_b3_ellipse(2,:),'--','Color',[0,0,.5],'LineWidth',1,'DisplayName','est obs (b) cov 3-sigma');
                        
                        plot(pred_a1_ellipse(1,:),pred_a1_ellipse(2,:),'r-.','LineWidth',3,'DisplayName','UT inference for p(at|xt,y0,zt=B) [1-sigma]');
%                         plot(pred_a3_ellipse(1,:),pred_a3_ellipse(2,:),'r-.','LineWidth',2,'DisplayName','UT inference for p(at|xt,y0,zt=B) [3-sigma]');
                        plot(pred_b1_ellipse(1,:),pred_b1_ellipse(2,:),'b-.','LineWidth',3,'DisplayName','UT inference for p(bt|xt,y0,zt=B) [1-sigma]');
%                         plot(pred_b3_ellipse(1,:),pred_b3_ellipse(2,:),'b-.','LineWidth',2,'DisplayName','UT inference for p(bt|xt,y0,zt=B) [3-sigma]');
                        
                        plot(est_axt3_ellipse(1,:),est_axt3_ellipse(2,:),'--','Color',[.8 .5 0],'LineWidth',1,'DisplayName','UT inference for p(axt|wt,y0,zt=B) [3-sigma]');
                        plot(est_xbt1_ellipse(1,:),est_xbt1_ellipse(2,:),'--','Color',[0 .5 .8],'LineWidth',2,'DisplayName','UT inference for p(bxt|wt,y0,zt=B) [1-sigma]');
%                         plot(est_xbt3_ellipse(1,:),est_xbt3_ellipse(2,:),'--','Color',[0 .5 .8],'LineWidth',1,'DisplayName','UT inference for p(bxt|wt,y0,zt=B) [3-sigma]');

                        plot(free_axt1_ellipse(1,:),free_axt1_ellipse(2,:),'-','Color',[.8 .5 0],'LineWidth',2,'DisplayName','p(axt|wt,y0,zt=F) [1-sigma]');
%                         plot(free_axt3_ellipse(1,:),free_axt3_ellipse(2,:),'-','Color',[.8 .5 0],'LineWidth',1,'DisplayName','UT inference for p(axt|wt,y0,zt=F) [3-sigma]');
                        plot(free_xbt1_ellipse(1,:),free_xbt1_ellipse(2,:),'-','Color',[0 .5 .8],'LineWidth',2,'DisplayName','p(bxt|wt,y0,zt=F) [1-sigma]');
%                         plot(free_xbt3_ellipse(1,:),free_xbt3_ellipse(2,:),'-','Color',[0 .5 .8],'LineWidth',1,'DisplayName','UT inference for p(bxt|wt,y0,zt=F) [3-sigma]');
                        
                        plot(y0(n,1),y0(n,2),'rd','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','ya(0) - prior position');
                        plot(y0(n,3),y0(n,4),'bd','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','yb(0) - prior position');
                        plot(est_y_mean(1),est_y_mean(2),'ro','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','UT-est ya(t)');
                        plot(est_y_mean(3),est_y_mean(4),'bo','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','UT-est yb(t)');
                        
                        plot(xt(n,1),xt(n,2),'rs','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','obs a(t)');
                        plot(xt(n,3),xt(n,4),'bs','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','obs b(t)');                        
                        plot(a1_ellipse(1,:),a1_ellipse(2,:),'-','Color',[1,0,0],'LineWidth',3,'DisplayName','obs xa(t) cov 1-sigma');
%                         plot(a3_ellipse(1,:),a3_ellipse(2,:),'-','Color',[1,0,0],'LineWidth',0.75,'DisplayName','obs xa cov 3-sigma');
                        plot(b1_ellipse(1,:),b1_ellipse(2,:),'-','Color',[0,0,1],'LineWidth',3,'DisplayName','obs xb(t) cov 1-sigma');
%                         plot(b3_ellipse(1,:),b3_ellipse(2,:),'-','Color',[0,0,1],'LineWidth',0.75,'DisplayName','obs xb cov 3-sigma');

                        legend('location','best');
                        axis('equal');
                        drawnow();
                    end
                end
            end


            p_zt_free = PairInteractionMCMCModel.pG_Z(zeros(N,1),t,z0,y0,params,RD); %P(zt=F) = 1-P(zt=B)
            p_zt_and_xt = [exp(llh_xt_free).*p_zt_free, exp(llh_xt_bound).*(1-p_zt_free)];
            % Bayes rule p(zt | xt wt z0 y0) = p(xt | zt wt z0 y0) p(zt | z0 y0) / p(xt | wt z0 y0), 
            %   where p(xt|w z0 y0)=sum_zt p(xt | zt wt z0 y0) p(zt | z0 y0).
            p_zt_given_xt = p_zt_and_xt ./ repmat(sum(p_zt_and_xt,2),1,2); 
            zt = rand(N,1) > p_zt_given_xt(:,1); %bound if greater than the free prob.
            pzt = p_zt_given_xt(:,1);
            pzt(zt) = p_zt_given_xt(zt,2);
%             sum(zt)
         end

         function [yt,llh,yt_mean,yt_cov]=sample_Q_YFree(xt,wt,t,y0,params)
            % Sample the prior proposal distribution for the free model using the Bayes's rule for 
            % linear gaussian systems.  Given that the prdictive distribution 
            % g(yt | zt=free, y0)=N(yt|y0,Q_t), Q_t=2*t*diag([D_AB,D_AB]); is gaussian 
            % and the measurement p(xt | yt, wt)=N(xt|yt,R), R=diag(wt.^2); is gaussian, the
            % posterior distribution p( yt | xt, wt, y0, zt=Free) = N(yt | yt_mean, yt_cov) is also gaussian,
            % And we can compute the posterior exactly.  So this is the optimal proposal distribution.
            % Sample:
            %  yt ~ p( yt | xt, wt, y0, zt=Free) 
            % Also report the log-likelihood of the sample, and the posterior distribution mean and cov for yt.
            % This method assumes zt=0 (free)
            % [in]
            %  xt - size:[N,4] observed positions at time t: [oax oay obx oby]
            %  wt - size:[N,4] observations sigmas at time t
            %   t - scalar t>0 time elapsed (i.e., \Delta t).
            %  y0 - size:[N,4] positions at time 0 [ax ay bx by]
            %  params - params struct
            % [out]
            %  yt - size:[N,4] y(t) samples
            %  llh - size:[N,1] covEllipse(x_mean, x_cov_chol, sigmaFactor, Ncirc)log-likelihood of yt
            %  yt_mean - size:[4,N] mean of posterior distribution p( yt | xt, wt, y0, zt=Free)
            %  yt_cov - size:[4,4,N] cov of posterior distribution p( yt | xt, wt, y0, zt=Free)
            %  NOTE: last two parameters are in a better order for column-major formatting which is transposed from the
            %        normal way of reporting yt vars in this class
            Q = 2*diag([params.DA, params.DA, params.DB, params.DB]); %predictive distribution covariance
            N = size(xt,1);
            if numel(t)==1
                t = t*ones(N,1);
            end
            yt = zeros(N,4);
            llh = zeros(N,1);
            yt_mean = zeros(4,N);
            yt_cov = zeros(4,4,N);
            for n=1:N
                Qt = t(n)*Q;
                Rt = diag(wt(n,:).^2);
                yt_cov(:,:,n) = inv(inv(Rt)+inv(Qt));
                yt_mean(:,n) = yt_cov(:,:,n)*(Rt\xt(n,:)' + Qt\y0(n,:)');
                yt(n,:) = (yt_mean(:,n) + chol(yt_cov(:,:,n))*randn(4,1))';
                llh(n)=PairInteractionMCMCModel.multiGaussianLLH(yt(n,:)'-yt_mean(:,n),yt_cov(:,:,n));
            end
        end

        function [yt,llh]=sample_Q_YBound(xt,wt,t,z0,y0,params)
            % Sample the prior proposal distribution for the bound model.
            % For the center position, we using the Bayes's rule for 
            % linear gaussian systems, simmilar to the free model, but with the appropriate distribution
            % For the distance r and angle theta, we sample from the proposal distribution, ignoring the estimate
            % of the 

            % Sample:
            %  yt ~ p( yt | xt, wt, y0, zt=Bound) 
            % Also report the log-likelihood of the sample, and the posterior distribution mean and cov for yt.
            % This method assumes zt=1 (bound)
            % [in]
            %  xt - size:[N,4] observed positions at time t: [oax oay obx oby]
            %  wt - size:[N,4] observations sigmas at time t
            %   t - scalar t>0 time elapsed (i.e., \Delta t).
            %  z0 - size:[N,1] state at time 0 (0=free, 1=bound)
            %  y0 - size:[N,4] positions at time 0 [ax ay bx by]
            %  params - params struct
            % [out]
            %  yt - size:[N,4] yt samples  (converted from vt)
            %  llh - size:[N,1] log-likelihood of yt (using -log(rt) as logabsdet(Vinv) as correction
            %  vt_mean - size:[4,N] mean of posterior distribution p( vt | xt, wt, y0, z0, zt=Bound)
            %  vt_cov - size:[4,4,N] cov of posterior distribution p( vt | xt, wt, y0, z0, zt=Bound)
            %  NOTE: last two parameters are in a better order for column-major formatting which is transposed from the
            %        normal way of reporting yt vars in this class
            [ct,llh_ct] = PairInteractionMCMCModel.sample_Q_Vc_Bound(xt,wt,t,y0,params);
            [rt,llh_rt] = PairInteractionMCMCModel.sample_Q_Vr_Bound(xt,wt,t,z0,y0,params);
            [phi_t,llh_phi_t] = PairInteractionMCMCModel.sample_Q_Vphi_Bound(xt,wt,t,z0,y0,params);
            vt=[ct, rt, phi_t];
            llh = llh_ct+llh_rt+llh_phi_t-log(rt); %correct for log(abs(det(V^{-1})))
            assert(all(isreal(llh)))
%             vt_mean = [ct_mean; rt_mean; phi_t_mean];
%             vt_cov = zeros(4,4,N);
%             for n=1:N
%                 vt_cov(:,:,n) = blkdiag(ct_cov(:,:,n), rt_var(n), phi_t_var(n));
%             end
            yt = PairInteractionMCMCModel.transVinv(vt);
            assert(isvector(llh));
        end

        function [ct,llh,ct_mean,ct_cov]=sample_Q_Vc_Bound(xt,wt,t,y0,params)
            % Sample the center position for the bound model based on the current measurement
             % [in]
            %  xt - size:[N,4] observed positions at time t: [oax oay obx oby]
            %  wt - size:[N,4] observations sigmas at time t
            %   t - scalar t>0 time elapsed (i.e., \Delta t).
            %  y0 - size:[N,4] positions at time 0 [ax ay bx by]
            %  params - params struct
            % [out]
            %  ct - size:[N,2] ct samples  
            %  llh - size:[N,1] log-likelihood of log(p(ct | xt, wt, y0, zt=Bound) 
            %  ct_mean - size:[2,N] mean of posterior distribution p( vt | xt, wt, y0, zt=Bound)
            %  ct_cov - size:[2,2,N] cov of posterior distribution p( vt | xt, wt, y0, zt=Bound)
            Qc = 2*diag([params.DAB, params.DAB]); %predictive distribution covariance
            C=[.5 0 .5 0; 0 .5 0 .5];
            N=size(xt,1);
            if numel(t)==1
                t = t*ones(N,1);
            end
            ct = zeros(N,2);
            llh = zeros(N,1);
            ct_mean = zeros(2,N);
            ct_cov = zeros(2,2,N);
            for n=1:N
                Qct = t(n)*Qc; %Covaraince for propogator
                Rct = C*diag(wt(n,:).^2)*C';
                xct = C*xt(n,:)';
                c0 = C*y0(n,:)';
                ct_cov(:,:,n) = inv(inv(Rct) + inv(Qct));
                ct_mean(:,n) = ct_cov(:,:,n)*(Rct\xct + Qct\c0);
                ct(n,:) = (ct_mean(:,n) + chol(ct_cov(:,:,n))*randn(2,1))';
                llh(n) = PairInteractionMCMCModel.multiGaussianLLH(ct(n,:)'-ct_mean(:,n),ct_cov(:,:,n));
            end
        end

        function [rt,llh,rt_mean,rt_var]=sample_Q_Vr_Bound(xt,wt,t,z0,y0,params)
            N=numel(z0);
            if numel(t)==1
                t = t*ones(N,1);
            end
            E = exp(-params.gamma*t)';
            rt_mean = zeros(1,N);
            rt_var = zeros(1,N);
            bound_0 = z0==1;
            Nbound_0 = sum(bound_0);
            Nfree_0 = N - Nbound_0;
            r0 = PairInteractionMCMCModel.transVr(y0)';
            if Nbound_0
                rt_mean(bound_0) = E(bound_0).*r0(bound_0) + (1-E(bound_0))*params.rho; 
                rt_var(bound_0) = params.sigma_rho^2*(1-E(bound_0).^2);
            end
            if Nfree_0
                rt_mean(~bound_0) = params.rho; 
                rt_var(~bound_0) = params.sigma_rho^2;
            end
            rt = abs(rt_mean' + randn(N,1).*sqrt(rt_var'));
            llh = PairInteractionMCMCModel.gaussianLLH(rt'-rt_mean,rt_var)';
        end

        function [phi_t, llh, mu, kappa]=sample_Q_Vphi_Bound(xt,wt,t,z0,y0,params)          
            N = numel(z0);
            [mu, kappa] = PairInteractionMCMCModel.estimateVonMisesKappa(xt(1,1:2)-xt(1,3:4), diag(wt(1,1:2)+wt(1,3:4)), 10);
            [phi_t,llh] = PairInteractionMCMCModel.sampleVonMises(mu,kappa,N);
%             figure();
%             histogram(phi_t,10,'Normalization','pdf','EdgeColor','none');
%             hold('on');
%             as=linspace(0,2*pi,1000);
%             apdfs = exp(PairInteractionMCMCModel.vonMisesLLH(mu-as,kappa));
%             plot(as,apdfs,'r-');
%             [phi_t2,llh2] = PairInteractionMCMCModel.sampleVonMises(mu,kappa,5000);
%             histogram(phi_t2,500,'Normalization','pdf','EdgeColor','none');
%             plot(phi_t2,exp(llh2),'o','MarkerFaceColor',[0,1,0],'MarkerEdgeColor',[0,0,0],'MarkerSize',3);
%             plot(phi_t,exp(llh),'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0]);
            phi_t=phi_t(:);
            llh=llh(:);
        end


        function [zt,yt,llh,llh_all] = sample_Q(xt,wt,t,z0,y0,params,RD)
            % Sample the prior proposal distribution
            %   (z0,y0) ~ q(z0, y0 | x0 w0) = q(z0 | y0) q(y0 | x0 w0).
            % Also report the log-likelihood of the sample.
            % [in]
            %  xt - size:[N,4] observed positions at time t: [oax oay obx oby]
            %  wt - size:[N,4] observations sigmas at time t
            %   t - scalar t>0 time elapsed (i.e., \Delta t).
            %  z0 - size:[N,1] state at time 0 (0=free, 1=bound)
            %  y0 - size:[N,4] positions at time 0 [ax ay bx by]
            %  params - params struct
            %  RD - RDCapture object 
            % [out]
            %  zt - size:[N] initial z(0)states 0=free 1=bound
            %  yt - size:[N,4] initial y(0) positions
            %  llh - size:[N,1] total (sum) log-likelihood of (z0(t), y0(t)) 
            N = numel(z0);
            if numel(t)==1
                t = t*ones(N,1);
            end
            if size(xt,1)==1 && N>1
                xt = repmat(xt(:)',N,1);
            end
            if size(wt,1)==1 && N>1
                wt = repmat(wt(:)',N,1);
            end
            llh_all = zeros(N,5);

            [zt,pzt] = PairInteractionMCMCModel.sample_Q_Z_X(xt,wt,z0,y0,t,params,RD); 
            bound = zt==1;
            Nbound = sum(bound);
            Nfree = N - Nbound;
            
            yt = zeros(N,4);
            llh = zeros(N,1);
            if Nfree
                [yt(~bound,:),llh(~bound)] = PairInteractionMCMCModel.sample_Q_YFree(xt(~bound,:),wt(~bound,:),...
                                                t(~bound),y0(~bound,:),params);
            end
            if Nbound
               [yt(bound,:),llh(bound)] = PairInteractionMCMCModel.sample_Q_YBound(xt(bound,:),wt(bound,:),...
                                                t(bound),z0(bound),y0(bound,:),params);
            end
            llh = llh + log(pzt);
        end


        %% Propogator g() sampling and llh computation
        function [z0,y0,llh,llh_all] = sample_G0(N,params,prior)
            % Sample the prior (z0,y0) ~ g(z(0), y(0))
            % [in]
            %  N - numer of samples
            %  params - params struct
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  y0 - size:[N,4] initial y(0) positions
            %  llh - [optional] size:[N,1] total (sum) log-likelihood of z0(t) y0(t)
            %  llh_all - [optional] size:[N,5] columns are log-likelihood of 
            %             cols are [z0, a0, -, b0, -] if z0=free
            %             cols are [z0, c0, -, r0+log(det(V^-1)), phi0] if z0=free
            z0 = PairInteractionMCMCModel.sample_G0_Z(N,prior);
            y0 = PairInteractionMCMCModel.sample_G0_Y(z0,params,prior);
            if nargout>2
                [llh,llh_all] = PairInteractionMCMCModel.llhG0(z0,y0,params,prior);
            end
        end

        function z0 = sample_G0_Z(N,prior)
            % Sample the prior z0 ~ g(z(0))
            % [in]
            %  N - numer of samples
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            z0 = rand(N,1) > prior.state(1);
        end


        function y0 = sample_G0_Y(z0,N,params,prior)
            % Sample the prior y0 ~ g(y(0) | z(0))
            % [in]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  params - params struct
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  y0 - size:[N,4] initial y(0) positions
            bound = z0==1;
            Nbound = sum(bound);
            Nfree = N - Nbound;
            if Nfree
                y0(~bound,1:2) = PairInteractionMCMCModel.samplePosPrior(Nfree, prior);
                y0(~bound,3:4) = PairInteractionMCMCModel.samplePosPrior(Nfree, prior);
            end
            if Nbound
                y0(bound,1:2) = PairInteractionMCMCModel.samplePosPrior(Nbound, prior);
                y0(bound,3) = params.rho + randn(Nbound,1).*params.sigma_rho;
                y0(bound,4) = 2*pi*rand(N,1);
            end
        end
        
        function [llh, llh_all] = llhG0(z0,y0,params,prior)
            % This is prior log-likelihood: log(g(z0,y0))
            % for a set of independent particles with known yt and params theta
            % [in]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  y0 - size:[N,4] initial y(0) positions
            %  params - params struct
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  llh - size:[N,1] total (sum) log-likelihood of z0(t) y0(t)
            %  llh_all - [optional] size:[N,5] columns are log-likelihood of 
            %             cols are [z0, a0, -, b0, -] if z0=free
            %             cols are [z0, c0, -, r0+log(det(V^-1)), phi0] if z0=free
            N = numel(z0);
            llh_all = zeros(N,5);
            llh_all(:,1) = PairInteractionMCMCModel.llhG0_Z(z0,prior);
            [~, llh_all(:,2:5)] = PairInteractionMCMCModel.llhG0_Y(y0,z0,params,prior);
            llh = sum(llh_all,2);
        end

        function llh = llhG0_Z(z0,prior)
            % This is prior log-likelihood: log(g(z0,y0))
            % for a set of independent particles with known yt and params theta
            % [in]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  llh - size:[N,1] total (sum) log-likelihood of z(t)
            llh = log(PairInteractionMCMCModel.pG0_Z(z0,prior));
        end

        function pzt = pG0_Z(z0,prior)
            % Probability of the state propogator g() for state zt
            % g0(z0) for the initial hidden state z0 using the prior.states  
            % [in]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  pzt - size:[N,1] probability of given zt
            pzt = prior.state(z0+1);
        end


        function [llh, llh_all] = llhG0_Y(y0,z0,params,prior)
            % This is prior log-likelihood: log(g(z0,y0))
            % for a set of independent particles with known yt and params theta
            % [in]
            %  y0 - size:[N,4] initial y(0) positions if size is [1,4] we expand it with N copies to [N,4]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  params - params struct
            %  prior - prior structure as defined by DefaultPrior
            % [out]
            %  llh - size:[N,1] total (sum) log-likelihood of z0(t) y0(t)
            %  llh_all - [optional] size:[N,5] columns are log-likelihood of 
            %             cols are [a0, -, b0, -] if z0=free
            %             cols are [c0, -, r0+log(det(V^-1)), phi0] if z0=free
            N = numel(z0);
            if numel(y0) == 4
                y0 = repmat(y0(:)',N,1);%if size is [1,4] (or [4,1]) we expand it with N copies to [N,4]
            end
            llh_all = zeros(N,4);
            
            bound = z0==1;
            Nbound = sum(bound);
            Nfree = N - Nbound;

            if Nfree
                llh_all(~bound,1) = .5*PairInteractionMCMCModel.llh_pos_prior(y0(~bound,1:2),prior);
                llh_all(~bound,2) = llh_all(~bound,1);
                llh_all(~bound,3) = .5*PairInteractionMCMCModel.llh_pos_prior(y0(~bound,3:4),prior);
                llh_all(~bound,4) = llh_all(~bound,3);
            end
            if Nbound
                v0 = PairInteractionMCMCModel.transV(y0(bound,:));
                r0 = v0(:,3);
                llh_all(bound,1) = .5*PairInteractionMCMCModel.llh_pos_prior(v0(1:2),prior);
                llh_all(bound,2) = llh_all(bound,1);
                llh_all(bound,3) = -log(r0) + ... %log(det(V^{-1}))
                                   PairInteractionMCMCModel.gaussianLLH(r0-params.rho, params.sigma_rho^2);
                llh_all(bound,4) = -PairInteractionMCMCModel.LOG2PI;
            end
            llh = sum(llh_all,2)';
        end

        function [zt,yt,llh,llh_all] = sample_G(z0,y0,t,params)
            % This is the forward state propogator.
            % Sample (z(t), y(t)) ~ g(z(t), y(t) | z(0) y(0) theta)
            % for a set of independent particles with known z0,y0,t, and params theta
            % [in]
            %  z0 - size:[N] initial z(0)states 0=free 1=bound
            %  y0 - size:[N,4] initial y(0) positions
            %  t - scalar or vector of times.  if vector must be size:[N] and z(i) is sampled at time t(i)
            %  params - params struct
            % [out]
            %  zt - size:[N] sampled z(t) states 0=free 1=bound
            %  yt - size:[N,4] sampled y(t) values according to propogator
            %  llh - [optional] size:[N,1] total (sum) log-likelihood of (z(t), y(t))
            %  llh_all - [optional] size:[N,5] individual llhs for each param columns are log-likelihood
            %            are [llh(zt), llh(axt), llh(ayt), llh(bxt), llh(byt)] if zt=free,
            %            and [llh(zt), llh(cxt), llh(cyt), llh(r)+logdet(V^-{1}), llh(theta)] if zt=bound
            zt = PairInteractionMCMCModel.sampleG_Z(z0,y0,t,params);
            yt = PairInteractionMCMCModel.sampleG_Y(zt,z0,y0,t,params);
            if nargout>2
                [llh,llh_all] = PairInteractionMCMCModel.llhG(zt,yt,t,z0,y0,params);
            end
        end

        function [llh,llh_all] = llhG(zt,yt,t,z0,y0,params,RD)
            % Log-likelihood of the state propogator g().
            % log(g(z(t),y(t) | z(0), y(0))) for the hidden states z(t) and y(t).  
            % [AndrieuDoucetHolenstein 2010 notation]
            % [in]
            %   zt - size:[N,1] state at time t 0=free 1=bound 
            %   yt - size:[N,4] Position at time t: [ax ay bx by] 
            %   t - time >0.  (This can also be thought of as dt=t_i-t_{i-1} for a sequence)
            %   z0 - size:[N,1] state at time 0: 0=free 1=bound 
            %   y0 - size:[N,4] Position at time 0: [ax ay bx by] 
            %   params - params struct
            %   RD - RDCapture object 
            % [out]
            %  llh - size:[N,1] total (sum) log-likelihood of (z(t), y(t))
            %  llh_all - [optional] size:[N,5] individual llhs for each param columns are log-likelihood
            %            are [llh(zt), llh(axt), llh(ayt), llh(bxt), llh(byt)] if zt=free,
            %            and [llh(zt), llh(cxt), llh(cyt), llh(r)+logdet(V^-{1}), llh(theta)] if zt=bound
            N = numel(zt);
            llh_all = zeros(N,5);
            llh_all(:,1) = PairInteractionMCMCModel.llhG_Z(zt,t,z0,y0,params,RD);
            [~, llh_all(:,2:5)] = PairInteractionMCMCModel.llhG_Y(yt,zt,t,z0,y0,params);
%             assert(all(isfinite(llh_all(:))));
            llh = sum(llh_all,2);
        end

        function [zt,pzt] = sampleG_Z(z0,y0,t,params,RD)
            % This is the forward state propogator for Z.
            % Sample z(t) ~ g(z(t) | z(0), y(0), theta) 
            % for a set of independent particles with known z0,y0 and given t, and params
            % [in]
            %  z0 - size:[N] initial states 0=free 1=bound
            %  y0 - size:[N,4] initial positions
            %  t - scalar or vector of times.  if vector must be size:[N] and z(i) is sampled at time t(i)
            %  params - params struct
            %  RD - RDCapture object 
            % [out]
            %  zt - size:[N,1] sampled zt values according to propogator zt ~ p( . | z0, y0)
            %  pzt - size:[N,1] probability: p(zt | z0, y0)
            N=numel(z0);
            zt = zeros(N,1);
            pzt = zeros(N,1);
            if numel(t)==1
                t=t*ones(N,1);
            end
            bound = z0==1;
            Nbound = sum(bound);
            Nfree = N - Nbound;
            if Nfree
                r0 = sqrt(sum((y0(~bound,1:2)-y0(~bound,3:4)).^2,2)); %distance
                pSurvive = RD.computeSurvivalProb(t(~bound),r0); %P to remain free
                zt(~bound) = rand(Nfree,1) > pSurvive; % Bind if the unit random variable > survival prob. 
                pzt(~bound) = pSurvive;
                captured = ~bound & zt; % true if was free and is now bound
                pzt(captured) = 1 - pzt(captured); % Flip state for captured
            end
            if Nbound
                pSurvive = exp(-params.koff*t(bound));
                zt(bound) = rand(Nbound,1) < pSurvive; % Was bound.  Will remain bound if...
                pzt(bound) = pSurvive;
                released = bound & ~zt; % true if was bound and now are free
                pzt(released) = 1-pzt(released); % flip the state
            end
        end

        function llh = llhG_Z(zt,t,z0,y0,params,RD)
            % Log-likelihood of the state propogator g() for state z
            % log(g(z(t) | z(0), y(0))) for the hidden states z(t) and y(t).  
            % [in]
            %   zt - size:[N,1] state at time t 0=free 1=bound 
            %   t - time >0.  (This can also be thought of as dt=t_i-t_{i-1} for a sequence)
            %   z0 - size:[N,1] state at time 0: 0=free 1=bound 
            %   y0 - size:[N,4] Position at time 0: [ax ay bx by] 
            %   params - params struct
            %   RD - RDCapture object 
            % [out]
            %  llh - size:[N,1] log-likelihood of g(z(t) | z0,y0,theta)
            llh = log(PairInteractionMCMCModel.pG_Z(zt,t,z0,y0,params,RD));
        end

        function pzt = pG_Z(zt,t,z0,y0,params,RD)
            % Probability of the state propogator g() for state zt
            % g(z(t) | z(0), y(0)) for the hidden states z(t).  
            % [in]
            %   zt - size:[N,1] state at time t 0=free 1=bound 
            %   t - time >0.  (This can also be thought of as dt=t_i-t_{i-1} for a sequence)
            %   z0 - size:[N,1] state at time 0: 0=free 1=bound 
            %   y0 - size:[N,4] Position at time 0: [ax ay bx by] 
            %   params - params struct
            %   RD - RDCapture object 
            % [out]
            %  pzt - size:[N,1] probability of given zt
            N=numel(zt);
            pzt = zeros(N,1);
            if numel(t)==1
                t=t*ones(N,1);
            end
            if numel(y0) == 4
                y0 = repmat(y0(:)',N,1);%if size is [1,4] (or [4,1]) we expand it with N copies to [N,4]
            end
            bound_0 = z0(:)==1;
            bound_t  = zt(:)==1;
            Nbound_0 = sum(bound_0);
            Nfree_0 = N - Nbound_0;
            if Nfree_0 % previously free so next state is determined by diffusion-to-capture
                r0 = sqrt(sum((y0(~bound_0,1:2)-y0(~bound_0,3:4)).^2,2)); %distance
                pzt(~bound_0) = RD.computeSurvivalProb(t(~bound_0),r0); % survival prob for each previously free particle
                pzt(~bound_0 & bound_t) = 1-pzt(~bound_0 & bound_t); % free->bound prob is 1- (free->free)
            end
            if Nbound_0 % previously bound so use exponential cdf
                pzt(bound_0) = exp(-params.koff*t(bound_0));
                pzt(bound_0 & ~bound_t) = 1-pzt(bound_0 & ~bound_t);% bound->free prob is 1- (bound->bound)
            end
        end


        function yt = sampleG_Y(zt,z0,y0,t,params)
            % This is the forward state propogator for Y given Z.
            % Sample y(t) ~ g( y(t) | z(t), z(0) y(0) theta)
            % for a set of independent particles with known z,z0,y0 and given t, and params
            % [in]
            %  zt - size:[N] current states 0=free 1=bound
            %  z0 - size:[N] initial states 0=free 1=bound
            %  y0 - size:[N,4] initial positions
            %  t - scalar or vector of times.  if vector must be size:[N] and z(i) is sampled at time t(i)
            %  params - params struct
            % [out]
            %  yt - size:[N,4] sampled y values according to propogator
            N = numel(z0);
            yt = zeros(N,4);
            sqrt_t = sqrt(t);
            if isscalar(sqrt_t)
                t=t*ones(N,1);
                sqrt_t=sqrt_t*ones(N,1);
            end
            bound_0 = z0==1;
            bound_t  = zt==1;
            Nbound_t = sum(bound_t);
            Nbound_0 = sum(bound_0);
            Nfree_t = N - Nbound_t;
            if Nfree_t %For currently free
                filt = ~bound_t;
                Nfilt = sum(filt);
                sigmaF = repmat(sqrt_t(filt),1,4).*repmat(sqrt(2*[params.DA,params.DA,params.DB,params.DB]),Nfilt,1);
                yt(filt,:) = y0(filt,:)+randn(Nfilt,4).*sigmaF;
            end
            if Nbound_t
                v0 = PairInteractionMCMCModel.transV(y0);
                if Nbound_0 %Bound -> Bound transition
                    filt = bound_t & bound_0;
                    Nfilt = sum(filt);
                    sigmaC = repmat(sqrt(2*params.DAB)*sqrt_t(filt),1,2);
                    yt(filt,1:2) = v0(filt,1:2) + randn(Nfilt,2).*sigmaC;  
                    E=exp(-params.gamma*t(filt));
                    sigmaR = params.sigma_rho.*sqrt(1-E.^2);
                    yt(filt,3) = params.rho*(1-E) + v0(filt,3).*E + randn(Nfilt,1).*sigmaR;
                    sigmaPhi = sqrt(2*params.Dphi)*sqrt_t(filt);
                    yt(filt,4) = v0(filt,4) + randn(Nfilt,1).*sigmaPhi;
                else % Free -> bound transition
                    filt = bound_t & ~bound_0;
                    Nfilt = sum(filt);
                    sigmaC = repmat(sqrt(2*params.DAB)*sqrt_t(filt),1,2);
                    yt(filt,1:2) = v0(filt,1:2) + randn(Nfilt,2).*sigmaC; 
                    yt(filt,3) = params.rho + randn(Nfilt,1)*params.sigma_rho; %At equilibrium
                    yt(filt,4) = rand(Nfilt,1)*2*pi; %uniform over [0,2pi)
                end
                yt(bound_t,:) = PairInteractionMCMCModel.transVinv(yt(bound_t,:)); %Transform bound V space back to Y space
            end                               
        end

        function [llh, llh_all] = llhG_Y(yt,zt,t,z0,y0,params)
            % NOTE: The order of parameters can be confusing.  we put yt first since the other
            %       variables zt, z0, y0 are conditioned upon.
            % Log-likelihood of the state propogator g() for state y(t)
            % log(g(y(t) | z(t), z(0), y(0))) for the hidden states y(t).  
            % [in]
            %   yt - size:[N,4] state at time t 0=free 1=bound 
            %   zt - size:[N,1] state at time t 0=free 1=bound 
            %   t - time >0.  (This can also be thought of as dt=t_i-t_{i-1} for a sequence)
            %   z0 - size:[N,1] state at time 0: 0=free 1=bound 
            %   y0 - size:[N,4] Position at time 0: [ax ay bx by] 
            %   params - params struct
            % [out]
            %  llh - size:[N,1] log-likelihood of g(z(t) | z0,y0, hiddetheta)
            %  llh_all - [optional] size:[N,4] individual llhs for each param columns are log-likelihood
            %            are [llh(zt), llh(axt), llh(ayt), llh(bxt), llh(byt)] if zt=free,
            %            and [llh(zt), llh(cxt), llh(cyt), llh(r)+logdet(V^-{1}), llh(theta)] if zt=bound
            N=numel(zt);
            llh_all = zeros(N,4);
            if numel(t)==1
                t=t*ones(N,1);
            end

            if numel(y0) == 4
                y0 = repmat(y0(:)',N,1);%if size is [1,4] (or [4,1]) we expand it with N copies to [N,4]
            end
            if numel(yt) == 4
                yt = repmat(yt(:)',N,1);%if size is [1,4] (or [4,1]) we expand it with N copies to [N,4]
            end
            bound_0 = z0(:)==1;
            bound_t  = zt(:)==1;
            Nbound_t = sum(bound_t);
            Nbound_0 = sum(bound_0);
            Nfree_t = N - Nbound_t;
            if Nfree_t %For currently free
                filt = ~bound_t;
                Nfilt = sum(filt);
                varF = repmat(t(filt),1,4).*repmat(2*[params.DA,params.DA,params.DB,params.DB],Nfilt,1);
                llh_all(filt,:) = PairInteractionMCMCModel.gaussianLLH(y0(filt,:)-yt(filt,:),varF);
            end
            if Nbound_t
                v0 = PairInteractionMCMCModel.transV(y0);
                vt = PairInteractionMCMCModel.transV(yt);
                if Nbound_0 %Bound -> Bound transition
                    filt = bound_t & bound_0;
                    varC = repmat(2*params.DAB*t(filt),1,2);
                    llh_all(filt,1:2) = PairInteractionMCMCModel.gaussianLLH(v0(filt,1:2)-vt(filt,1:2),varC);  
                    E=exp(-params.gamma*t(filt));
                    varR = params.sigma_rho^2.*(1-E.^2);
                    muR = params.rho*(1-E) + v0(filt,3).*E;
                    llh_all(filt,3) =  PairInteractionMCMCModel.gaussianLLH(muR - vt(filt,3),varR);
                    kappa = 1./(2*t(filt)*params.Dphi);
                    llh_all(filt,4) = PairInteractionMCMCModel.vonMisesLLH(v0(filt,4) - vt(filt,4), kappa); %note factor of 2 in diffusion const is accounted for in vonMises defn.
                else % Free -> bound transition
                    filt = bound_t & ~bound_0;
                    varC = repmat(2*params.DAB*t(filt),1,2);
                    llh_all(filt,1:2) = PairInteractionMCMCModel.gaussianLLH(v0(filt,1:2)-vt(filt,1:2),varC); 
                    llh_all(filt,3) = PairInteractionMCMCModel.gaussianLLH(params.rho -vt(filt,3), params.sigma_rho^2); %At equilibrium
                    found = find(filt);
                    for ni=1:numel(found)
                        n = found(ni);
                        center = v0(n,1:2);
                        cholC = eye(2)*sqrt(2*t(n)*(params.DA+params.DB));
                        [muPhi, kappaPhi] = PairInteractionMCMCModel.estimateVonMisesKappa(center, cholC);
                        llh_all(filt,4) = PairInteractionMCMCModel.vonMisesLLH(vt(n,4)-muPhi,kappaPhi);
                    end
                end
                llh_all(bound_t,3) = llh_all(bound_t,3) - log(vt(bound_t,3)); % log(det(V^{-1}))
            end
            llh = sum(llh_all,2);
        end

        %% Measurement h() sampling and llh computation
        function xt = sampleH(wt,yt)
            % This is the observation distribution observed pos X(t) given true pos y(t) and 
            %  observed localization std sigma w(t).
            % We include overlay error along X and Y directions.
            % Sample x(t) ~ h( x(t) | w(t), y(t), theta)
            % for a set of independent particles with known yt and params theta
            % [in]
            %  wt - size:[N,4] standard deviation of observerd postion [sax, say, sbx, sby]
            %                  (including registration sigma)
            %  yt - size:[N,4] true positions [ax ay bx by]
            % [out]
            %  xt - size:[N,4] sampled x(t) according to observation
            N = size(wt,1);
            xt = yt + randn(N,4).*repmat(wt,N,1);
        end

        function [llh, llh_all] = llhH(xt,wt,yt)
            % This is the observation distribution observed pos X(t) given true pos y(t) and 
            %  observed localization std sigma w(t).
            % We include overlay error along X and Y directions.
            % Return log(h(x(t) | w(t), y(t))
            % for a set of independent particles with known yt and params theta
            % [in]
            %  xt - size:[N,4] observed positions [oax oay obx oby]
            %  wt - size:[N,4] observed standard deviation of postion [sax, say, sbx, sby]
            %                  (including registration sigma)
            %  yt - size:[N,4] true positions [ax ay bx by]
            % [out]
            %  llh - size:[N,1] total (sum) log-likelihood of x(t)
            %  llh_all - [optional] size:[N,4] log-likelihood of x(t) by component
            N = size(yt,1);
            if size(xt,1)==1 && N>1
                xt = repmat(xt(:)',N,1);
            end
            if size(wt,1)==1 && N>1
                wt = repmat(wt(:)',N,1);
            end
            llh_all = PairInteractionMCMCModel.gaussianLLH(xt-yt,wt.^2);
            llh = sum(llh_all,2);
        end


        function [P,w_norm,llhP,llh_obs] = condensationParticleFilter(obs, N, params, prior,plotOn)
            % This is a very simple 1-data point particle filter that uses the
            % condensation particle filter (i.e. using the state propogators as the proposal
            % distribution)
            % [in]
            %  obs - data observation matrix [t ax ay bx by sax say sbx sby]
            %  N - number of particles to sample
            %  params - params struct to use for model
            %  prior - prior struct to use
            % [out]
            %  P - size:[T,5,N] sampled particles columsn are [z ax ay bx by]
            %  w_norm - size:[N] normalized importance weights.
            %  llhP - size:[N] log-likelihood of each particle given obs
            %  llh_obs - size:[1] probability of the observation p(x_{1:T}, w_{1:T}) estimated be weighting
            if nargin<5
                plotOn=false;
            end
            T = size(obs,1);
            times = obs(:,1);
            xs = obs(:,2:5);
            vxs = PairInteractionMCMCModel.transV(xs);
            ws = obs(:,6:9);
            dt = diff(times);
            RD = RDCapture(params.DA+params.DB, params.rho_bind, params.lambda_bind);
            
            P = zeros(T,5,N);
            llhP = zeros(T,N); %likelihood contribution from each time t for each particle N.  sum(llhP)=overall llh for particles 1:N            
            llhQ = zeros(T,N); %likelihood contribution from each time t for each particle N.  sum(llhP)=overall llh for particles 1:N            
            log_alpha = zeros(T,N); %incrimental importance for time t of paritcle N

            %Initialize prior
            [z0,y0,llhQ(1,:)] = PairInteractionMCMCModel.sample_Q0(xs(1,:),ws(1,:),N,params,prior);
            P(1,1,:) = z0;
            P(1,2:5,:) = y0';
            llhG = PairInteractionMCMCModel.llhG0(z0,y0,params,prior);
            llhH = PairInteractionMCMCModel.llhH(xs(1,:),ws(1,:),y0);
            llhP(1,:) = llhG+llhH;
            log_alpha(1,:) =llhP(1,:)-llhQ(1,:);
            logw=log_alpha(1,:); % Log weights
            w_norm = logNormalize(logw); %normalized weights
            if plotOn
                figure('Position',[10,10,1200,1000]);
            end
            for t=2:T
                if plotOn
                    subplot(4,2,1);
                    cla();
                    hold('on');                    
                    for tt=max(1,t-2):t-1
                        if tt==t-1
                            color_a=[1,0,0];
                            color_b=[0,0,1];
                        else
                            color_a=[.5,0,0];
                            color_b=[0,0,.5];
                        end
                        plot(squeeze(P(tt,2,:)),squeeze(P(tt,3,:)),'.','MarkerEdgeColor',color_a,'DisplayName','a0 particles');
                        plot(squeeze(P(tt,4,:)),squeeze(P(tt,5,:)),'.','MarkerEdgeColor',color_b,'DisplayName','b0 particles');
                        if tt>2
                            plot(xs(tt-1:tt,1),xs(tt-1:tt,2),'-','Color',color_a,'LineWidth',0.75);
                            plot(xs(tt-1:tt,3),xs(tt-1:tt,4),'-','Color',color_b,'LineWidth',0.75);
                        end
                        plot(xs(tt,1),xs(tt,2),'s','MarkerFaceColor',color_a,'MarkerEdgeColor',[0,0,0],'DisplayName','obs_a0');
                        plot(xs(tt,3),xs(tt,4),'s','MarkerFaceColor',color_b,'MarkerEdgeColor',[0,0,0],'DisplayName','obs_b0');
                        a_cov_chol = chol(diag(ws(tt,1:2).^2));
                        a1_ellipse = PairInteractionMCMCModel.covEllipse(xs(tt,1:2), a_cov_chol, 1);
                        a3_ellipse = PairInteractionMCMCModel.covEllipse(xs(tt,1:2), a_cov_chol, 3);
                        b_cov_chol = chol(diag(ws(tt,3:4).^2));
                        b1_ellipse = PairInteractionMCMCModel.covEllipse(xs(tt,3:4), b_cov_chol, 1);
                        b3_ellipse = PairInteractionMCMCModel.covEllipse(xs(tt,3:4), b_cov_chol, 3);
                        plot(a1_ellipse(1,:),a1_ellipse(2,:),'-','LineWidth',3,'Color',color_a,'DisplayName','obs_a0 cov 1-sigma');
                        plot(a3_ellipse(1,:),a3_ellipse(2,:),'-','LineWidth',0.75,'Color',color_a,'DisplayName','obs_a0 cov 3-sigma');
                        plot(b1_ellipse(1,:),b1_ellipse(2,:),'-','LineWidth',3,'Color',color_b,'DisplayName','obs_b0 cov 1-sigma');
                        plot(b3_ellipse(1,:),b3_ellipse(2,:),'-','LineWidth',0.75,'Color',color_b,'DisplayName','obs_b0 cov 3-sigma');
                    end
                    axis('equal');
                    subplot(4,2,2);
                    hold('off');
                    cla();
                    a_cov_chol = chol(diag(ws(t-1,1:2).^2));
                    a1_ellipse = PairInteractionMCMCModel.covEllipse(xs(t-1,1:2), a_cov_chol, 1);
                    a3_ellipse = PairInteractionMCMCModel.covEllipse(xs(t-1,1:2), a_cov_chol, 3);
                    b_cov_chol = chol(diag(ws(t-1,3:4).^2));
                    b1_ellipse = PairInteractionMCMCModel.covEllipse(xs(t-1,3:4), b_cov_chol, 1);
                    b3_ellipse = PairInteractionMCMCModel.covEllipse(xs(t-1,3:4), b_cov_chol, 3);

                    hold('on');
                    for n=1:N
                        if P(t-1,1,n)==1
                            h=surface([P(t-1,[2 2],n); P(t-1,[4 4],n)],[P(t-1,[3 3],n); P(t-1,[5 5],n)],times([1 1;1 1]),'EdgeColor',[0,1,0],'EdgeAlpha',0.75*sqrt(w_norm(n))+0.25); 
                            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        end
                    end
                    plot(squeeze(P(t-1,2,:)),squeeze(P(t-1,3,:)),'r.','DisplayName','a0 particles');
                    plot(squeeze(P(t-1,4,:)),squeeze(P(t-1,5,:)),'b.','DisplayName','b0 particles');
                    plot(xs(t-1,1),xs(t-1,2),'s','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0],'DisplayName','obs_a0');
                    plot(xs(t-1,3),xs(t-1,4),'s','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,0],'DisplayName','obs_b0');
                    plot(a1_ellipse(1,:),a1_ellipse(2,:),'r-','LineWidth',3,'DisplayName','obs_a0 cov 1-sigma');
                    plot(a3_ellipse(1,:),a3_ellipse(2,:),'r-','LineWidth',0.75,'DisplayName','obs_a0 cov 3-sigma');
                    plot(b1_ellipse(1,:),b1_ellipse(2,:),'b-','LineWidth',3,'DisplayName','obs_b0 cov 1-sigma');
                    plot(b3_ellipse(1,:),b3_ellipse(2,:),'b-','LineWidth',0.75,'DisplayName','obs_b0 cov 3-sigma');
                    axis('equal');

                    subplot(4,2,3);
                    hold('off');
                    cla();
                    hold('on');
                    for p=1:N
                        if P(t-1,1,p)==1
                            color=[1,0,0];
                        else
                            color=[0,0,1];
                        end
                        plot(1:t-1,llhP(1:t-1,p),'-','Color',color);                    
                    end
                    title('llhP');

                    subplot(4,2,4);
                    hold('off');
                    cla();
                    hold('on');
                    for p=1:N
                        if P(t-1,1,p)==1
                            color=[1,0,0];
                        else
                            color=[0,0,1];
                        end
                        plot(1:t-1,llhQ(1:t-1,p),'-','Color',color);     
                    end
                    title('llhQ');

                    subplot(4,2,5);
                    hold('off');
                    cla();
                    hold('on');
                    for p=1:N
                        if P(t-1,1,p)==1
                            color=[1,0,0];
                        else
                            color=[0,0,1];
                        end
                        plot(1:t-1,log_alpha(1:t-1,p),'-','Color',color);
                        if w_norm(p)>0.01
                            plot(t-1,log_alpha(t-1,p),'o','MarkerFaceColor',color,'MarkerEdgeColor','None','MarkerSize',sqrt(w_norm(p))*10);
                        end
                    end
                    title('log_alpha');

                    subplot(4,2,6);
                    cla();
                    hold('on');
                    if t>2
                        ds = sqrt((squeeze(P(1:t-1,2,:))-squeeze(P(1:t-1,4,:))).^2+(squeeze(P(1:t-1,3,:))-squeeze(P(1:t-1,5,:))).^2);
                        for p=1:N                    
                            plot(1:t-1,ds(:,p)','-');
                        end
                    end
                    plot(1:t-1,vxs(1:t-1,3),'k-','linewidth',3);
                    plot([1,t-1],params.rho*[1,1],'b--','linewidth',1);
%                     yl=ylim();
%                     
%                     ylim([0,yl(2)]);
                    ylim('auto');
                    title('Distance: r');

                    subplot(4,2,7);
                    cla();
                    hold('on');
                    if t>2
                        stairs(1:t-1,sum(repmat(w_norm,t-1,1).*squeeze(P(1:t-1,1,:)),2),'-');
                    end
                    ylim([0,1.1]);
                    title('state');
                    subplot(4,2,8);
                    cla();
                    hold('on');
                    histogram(w_norm,N);

                    drawnow();
                end

                %Resample
                ess = 1/sum(w_norm.^2);
                if ess<max(3,N/20)
                    if plotOn
                        fprintf('T:%i N:%i ESS:%.3g ',t,N,ess);
                        fprintf('Resampling!\n');
                    end
                    idxs = randsample(N,N,true,w_norm); %random weighted sample of N from N with replacement!                
                    temp_P = zeros(t-1,5,N);
                    temp_llhP = zeros(t-1,N);
                    temp_llhQ = zeros(t-1,N);
                    temp_log_alpha = zeros(t-1,N);
                    for n=1:N
                        temp_P(:,:,n) = P(1:t-1,:,idxs(n));
                        temp_llhP(:,n) = llhP(1:t-1,idxs(n));
                        temp_llhQ(:,n) = llhQ(1:t-1,idxs(n));
                        temp_log_alpha(:,n) = log_alpha(1:t-1,idxs(n));                    
                    end
                    P(1:t-1,:,:) = temp_P;
                    llhP(1:t-1,:) = temp_llhP;
                    llhQ(1:t-1,:) = temp_llhQ;
                    log_alpha(1:t-1,:) = temp_log_alpha;
                    logw = logw(idxs);
                end
                            
                %Compute next steps
                z0 = squeeze(P(t-1,1,:));
                y0 = squeeze(P(t-1,2:5,:))';
                [zt,yt,llhQ(t,:)] = PairInteractionMCMCModel.sample_Q(xs(t,:),ws(t,:), dt(t-1),z0,y0,params,RD);
                
                P(t,1,:) = zt;
                P(t,2:5,:) = yt';
                llhG = PairInteractionMCMCModel.llhG(zt,yt,dt(t-1),z0,y0,params,RD);
                llhH = PairInteractionMCMCModel.llhH(xs(t,:),ws(t,:),yt);
                llhP(t,:) = llhG+llhH;

                log_alpha(t,:) = llhP(t,:)-llhQ(t,:);
                logw = logw + log_alpha(t,:);
                logw(~isfinite(logw)) = -inf; %Remove any NAN or +inf errors
                w_norm = logNormalize(logw); %normalized weights
                s_w_norm = sum(w_norm);
                
%                 if max(log_alpha(t,:))>max(log_alpha(t-1,:))
%                     [~,idx] = max(log_alpha(t,:));
%                     fprintf('id:%i [State:%i->%i] llhP:%.6g llhQ:%.6g logw:%.6g w_norm:%.6g \n',idx,P(t-1,1,idx),P(t,1,idx),llhP(t-1,idx),llhQ(t-1,idx),logw(idx), w_norm(idx));
%                 end
                
                assert(abs(1-s_w_norm) < sqrt(eps));
                assert(isreal(w_norm));
            end
            llhP = sum(llhP,1);
            llh_obs = logSum(logw)-log(N); % llh_obs(t) = log(p(x(1:T))).  it is the average weight
        end


        function [P,w_norm,llhP,llh_obs] = condensationStatesParticleFilter(ts, ys, N, params, prior, plotOn)
            % This is a very simple 1-data point particle filter that uses the
            % condensation particle filter assuming the Ys the true positions are known exactly.
            % [in]
            %  ts - size:[T,1] times
            %  ys - size:[T,4] true positionms [ax ay bx by]
            %  N - number of particles to sample
            %  params - params struct to use for model
            %  prior - prior struct to use
            % [out]
            %  P - size:[T,N] sampled particles state 0=free 1=bound at time T for particle N
            %  w_norm - size:[N] normalized importance weights.
            %  llhP - size:[T,N] log-likelihood of each particle z(t) given observed y(t) at time T for particle N
            %  pobs - size:[T] probability of the observation p(y_{1:t}) for each estimated as the unnormalized weight
            
            T = numel(ts);
            RD = RDCapture(params.DA+params.DB, params.rho_bind, params.lambda_bind);
                
            dt = diff(ts);
            P = zeros(T,N);
            llhZ = zeros(T,N); % for debugging
            llhY = zeros(T,N); % for debugging
            llhQ = zeros(T,N); % for debugging
            llhP = zeros(T,N); %likelihood contribution from each time t for each particle N.  sum(llhP)=overall llh for particles 1:N            
            log_alpha = zeros(T,N); %incrimental importance for time t of paritcle N

            %Initialize prior
            y0=repmat(ys(1,:),N,1);
            [z0,p_q0_z0] = PairInteractionMCMCModel.sample_Q0_Z(y0,params,prior);
            P(1,:) = z0;
            p_g0_z = PairInteractionMCMCModel.pG0_Z(z0,prior); %prior on state
            
            llhY0 = PairInteractionMCMCModel.llhG0_Y(ys(1,:),z0,params,prior);
            llhQ(1,:) = log(p_q0_z0);
            llhY(1,:) = llhY0;
            llhZ(1,:) = log(p_g0_z);
            llhP(1,:) = log(p_g0_z)+llhY0; %prior for p(z0,y0) = p(y0|z0)p(z0) 
            log_alpha(1,:) = llhP(1,:)-log(p_q0_z0); % incrimental importance weight
            logw=log_alpha(1,:); % Log weights
            w_norm = exp(logw-logSum(logw)); %normalized weights
            for t=2:T
%                 if nargin>5 && plotOn && t==T
                  if t==T
                    figure();
                    ds=sqrt(sum((ys(:,1:2)-ys(:,3:4)).^2,2));
                    trng = 1:t-1;
                    subplot(4,2,1);
                    hold('off');
                    plot(trng,ds(trng),'r-','DisplayName','distance');
                    hold('on');
                    yl=ylim();
                    plot(trng,yl(2)*mean(P(trng,:),2),'-b','DisplayName','unweighted average particle state');
                    plot(trng,yl(2)*sum(P(trng,:).*repmat(w_norm(:)',t-1,1),2),'-k','DisplayName','mean particle state');
                    legend('Location','Best');
                    ylabel('dist');
                    xlabel('time');
                    subplot(4,2,2);
                    hold('off');
                    for n=1:N
                        plot(trng,llhZ(trng,n),'-');
                        hold('on');
                    end
                    ylabel('llhZ');
                    xlabel('time');
                    subplot(4,2,3);
                    hold('off');
                    for n=1:N
                        plot(trng,llhY(trng,n),'-');
                        hold('on');
                    end
                    ylabel('llhY');
                    xlabel('time');
                    subplot(4,2,4);
                    hold('off');
                    for n=1:N
                        plot(trng,llhP(trng,n),'-');
                        hold('on');
                    end
                    ylabel('llhP');
                    xlabel('time');

                    subplot(4,2,5);
                    hold('off');
                    for n=1:N
                        plot(trng,llhQ(trng,n),'-');
                        hold('on');
                    end
                    ylabel('llhQ');
                    xlabel('time');

                    subplot(4,2,6);
                    hold('off');
                    for n=1:N
                        plot(trng,log_alpha(trng,n),'-');
                        hold('on');
                    end
                    ylabel('log_alpha (inc. importance)');
                    xlabel('time');

                    subplot(4,2,7);
                    hold('off');
                    histogram(log(w_norm),N,'DisplayName',sprintf('Normalized log-importance weights step:%i var:%.6g',t-1,var(w_norm)));
                    xlabel('normalized weights');
                    legend('Location','Best');
                end

                %Resample
                ess = 1/sum(w_norm.^2);
%                 fprintf('T:%i N:%i ESS:%.3g ',t,N,ess);
                if ess<max(3,N/20)
%                     fprintf('Resampling!\n');
                    idxs = randsample(N,N,true,w_norm); %random weighted sample of N from N with replacement!                
                    temp_P = zeros(t-1,N);
                    temp_llhQ = zeros(t-1,N);
                    temp_llhZ = zeros(t-1,N);
                    temp_llhY = zeros(t-1,N);
                    temp_llhP = zeros(t-1,N);
                    temp_log_alpha = zeros(t-1,N);
                    for n=1:N
                        temp_P(:,n) = P(1:t-1,idxs(n));
                        temp_llhQ(:,n) = llhQ(1:t-1,idxs(n));
                        temp_llhZ(:,n) = llhZ(1:t-1,idxs(n));
                        temp_llhY(:,n) = llhY(1:t-1,idxs(n));
                        temp_llhP(:,n) = llhP(1:t-1,idxs(n));
                        temp_log_alpha(:,n) = log_alpha(1:t-1,idxs(n));                    
                    end
                    P(1:t-1,:) = temp_P;
                    llhQ(1:t-1,:) = temp_llhQ;
                    llhZ(1:t-1,:) = temp_llhZ;
                    llhY(1:t-1,:) = temp_llhY;
                    llhP(1:t-1,:) = temp_llhP;
                    log_alpha(1:t-1,:) = temp_log_alpha;
                    logw = logw(idxs);
                else
%                     fprintf('<Continue>\n');
                end

                %Compute next steps
                z0 = P(t-1,:);
                y0 = ys(t-1,:);
                yt = ys(t,:);
                [zt,p_q_zt] = PairInteractionMCMCModel.sample_Q_Z(yt,z0,y0,dt(t-1),params,RD);
                p_g_zt = PairInteractionMCMCModel.pG_Z(zt,dt(t-1),z0,y0,params,RD);
                llh_Y = PairInteractionMCMCModel.llhG_Y(yt,zt,dt(t-1),z0,y0,params);
                P(t,:) = zt;         
                llhQ(t,:) = log(p_q_zt);
                llhZ(t,:) = log(p_g_zt);
                llhY(t,:) = llh_Y;
                llhP(t,:) = llh_Y+log(p_g_zt);
                log_alpha(t,:) = llhP(t,:)-log(p_q_zt);
                logw = logw + log_alpha(t,:);
                logw(~isfinite(logw)) = -inf; %Make sure any NaN's are given 0 weight
                w_norm = exp(logw-logSum(logw)); %normalized weights
                assert(abs(1-sum(w_norm))<sqrt(eps));                
            end
            llhP = sum(llhP,1);
            llh_obs = logSum(logw)-log(N); % llh_obs(t) = log(p(x(1:T))).  it is the average weight
        end

        function [marginal_ys_mean, marginal_ys_cov, llh_xs, pred_ys_mean, pred_ys_cov] = computeKalmanFilter_Free(obs,params,prior,plotOn)
            % Run the kalman filter for the free model predicting the sequence of marginal distributions
            % p(x_t | y_{1:t}, w_{1:t}) as well as an estimate of the likelihood p(x_{1:T} | w_{1:T})
            % [out]
            %
            %
            %  llh_xs = log(p(x_{1:T} | w_{1:T}));
            N = size(obs,1);
            ts = obs(:,1);
            dt = diff(ts);
            xs = obs(:,2:5);
            ws = obs(:,6:9);
            marginal_ys_mean=zeros(4,N); % marginal in time p(y(t) |  x_(1:t), w_(1:t))
            marginal_ys_cov=zeros(4,4,N);
            pred_ys_mean = zeros(4,N); % predictive distribution  p(y(t) |  x_(1:t-1), w_(1:t-1))
            pred_ys_cov = zeros(4,4,N);
            llh_xs = 0;
            Q = 2*diag([params.DA, params.DA, params.DB, params.DB]);
            [pred_ys_mean(:,1),pred_ys_cov(:,:,1)] = PairInteractionMCMCModel.pos_prior_distribution_free(prior);
            for n=1:N
                S = pred_ys_cov(:,:,n) + diag(ws(n,:).^2);
                K = pred_ys_cov(:,:,n)/S; % same as pred_yt_cov*inv(K)
                pred_xs_mean = pred_ys_mean(:,n); %This is pedantic to remind me how this works.  We simply have no input or nonlinearity here so equality is exact.  Remove later for speed.
                residual = xs(n,:)' - pred_xs_mean; % Difference between predicted and actual data
                marginal_ys_mean(:,n) = pred_ys_mean(:,n) + K*residual;
                marginal_ys_cov(:,:,n) = (eye(4) - K)*pred_ys_cov(:,:,n);
                llh_xs = llh_xs + PairInteractionMCMCModel.multiGaussianLLH(pred_ys_mean(:,n) - xs(n,:)', S);
                %Compute the predicted values for the next iter using this iter's values skip on last iter
                if n<N
                    pred_ys_mean(:,n+1) = marginal_ys_mean(:,n);
                    pred_ys_cov(:,:,n+1) = marginal_ys_cov(:,:,n) + dt(n)*Q;
                end
            end
        end

        function [marginal_ys_mean, marginal_ys_cov, llh_xs] = computeKalmanSmoother_Free(obs,params,prior)
            % Run the kalman filter for the free model predicting the sequence of marginal distributions
            % p(x_t | y_{1:t}, w_{1:t}) as well as an estimate of the likelihood p(x_{1:T} | w_{1:T})
            % [out]
            %
            %
            %  llh_xs = log(p(x_{1:T} | w_{1:T}));
            N = size(obs,1);

            [Fmarginal_ys_mean, Fmarginal_ys_cov, llh_xs, pred_ys_mean, pred_ys_cov] = ...
                    PairInteractionMCMCModel.computeKalmanFilter_Free(obs,params,prior);
            marginal_ys_mean=zeros(4,N); % smoothed marginal in time p(y(t) |  x_(1:T), w_(1:T))
            marginal_ys_cov=zeros(4,4,N);
            marginal_ys_mean(:,N) = Fmarginal_ys_mean(:,N); %For both filtering and smoothing last element is the same
            marginal_ys_cov(:,N) = Fmarginal_ys_cov(:,N);
            for n=N-1:-1:1
                J = Fmarginal_ys_cov(:,:,n)/pred_ys_cov(:,:,n+1);
                marginal_ys_mean(:,n) = Fmarginal_ys_mean(:,n) + J*(marginal_ys_mean(:,n+1)-pred_ys_mean(:,n+1)); 
                marginal_ys_cov(:,:,N) = Fmarginal_ys_cov(:,:,N) +J*(marginal_ys_cov(:,:,n+1) - pred_ys_cov(:,:,n+1))*J';
            end
        end

        function [marginal_vs_mean, marginal_vs_cov, llh_xs, pred_vs_mean, pred_vs_cov] = computeUnscentedKalmanFilter_Bound(obs,params,prior,plotOn)
            % Run the kalman filter for the bound model predicting the sequence of marginal distributions
            % p(x_t | y_{1:t}, w_{1:t}) as well as an estimate of the likelihood p(x_{1:T} | w_{1:T})
            % [out]
            %
            %
            %  llh_xs = log(p(x_{1:T} | w_{1:T}));
            N = size(obs,1);
            ts = obs(:,1);
            dt = diff(ts);
            xs = obs(:,2:5);
            ws = obs(:,6:9);
            % marginal in time p(v(t) |  x_(1:t), w_(1:t)) this is what we are computing.  
            % Gaussian in V-space
            marginal_vs_mean = zeros(4,N); 
            marginal_vs_cov = zeros(4,4,N);
            % predictive density  p(v(t) |  x_(1:t-1), w_(1:t-1)).
            % Gaussian in V-space
            pred_vs_mean = zeros(4,N); 
            pred_vs_cov = zeros(4,4,N);
            llh_xs = 0; %log-likelihood of the data log(p(x_{1:T} | w_{1:T})).  Built up temporally.
            [pred_vs_mean(:,1),pred_vs_cov(:,:,1)] = PairInteractionMCMCModel.pos_prior_distribution_bound(params,prior);
            UT=UnscentedTransform();
            for n=1:N
                %Approximate p(x_t | y_t)=N(x_t| est_y_mean, S) as normal with UT 
                [~,p] = chol(pred_vs_cov(:,:,n));
                assert(p==0);
                [est_xs_mean, est_xs_cov, x_sigmaPts] = UT.estimate(@PairInteractionMCMCModel.transVinv_row, pred_vs_mean(:,n), pred_vs_cov(:,:,n));
                S = est_xs_cov+diag(ws(n,:).^2);
                %Compute the cross-covariance using the x_sigmaPts and sigma pts derived from the exact values 
                % for pred_vs in V-space.
                cross_cov = UT.estimateCrossModelCov(x_sigmaPts, est_xs_mean,  pred_vs_mean(:,n), pred_vs_cov(:,:,n));
                K = cross_cov/S;
                residual = xs(n,:)' - est_xs_mean;
                % Use Bayes's rule for gaussians to get p(y_t | x_{1:t}, w_{1:t}) the marginal filtered probability
                marginal_vs_mean(:,n) = pred_vs_mean(:,n) + K*residual;
                marginal_vs_cov(:,:,n) =  pred_vs_cov(:,:,n) - K*S*K';              
                llh_xs = llh_xs + PairInteractionMCMCModel.multiGaussianLLH(PairInteractionMCMCModel.transVinv_row(pred_vs_mean(:,n)) - xs(n,:)', S);
                %Compute the predicted values for the next iter using this iter's values skip on last iter
                if n<N
                    E = exp(-params.gamma*dt(n));
                    Q = diag([2*dt(n)*params.DAB, 2*dt(n)*params.DAB, params.sigma_rho^2*(1-E^2), 2*dt(n)*params.Dphi]);
                    pred_vs_mean(:,n+1) = marginal_vs_mean(:,n).*[1;1;E;1] + [0;0;params.rho*(1-E);0]; % A_t*mu_t+B_t*u_t - update handles OU process for r
                    pred_vs_cov(:,:,n+1) = marginal_vs_cov(:,:,n) + Q;
                end
            end
        end
        
        function debugRDCaptureSurvivalProb(obs,p,params)
            RD = RDCapture(params.DA+params.DB, params.rho_bind, params.lambda_bind);

            times = obs(:,1);
            T = size(obs,1); %Number of times
            dt = diff(times);
            ys = p(:,2:5);
            vs = PairInteractionMCMCModel.transV(ys); %Transform to V-space
            rs = vs(:,3);

            survivalProb = RD.computeSurvivalProb(rs(1:end-1),dt);
            survivalProb_C = RD.computeSurvivalProb_C(rs(1:end-1),dt);
            mu =  RD.compute_mu(rs(1:end-1),dt);
            mu_C =  RD.compute_mu_C(rs(1:end-1),dt);
            nu =  RD.compute_nu(dt);
            nu_C =  RD.compute_nu_C(dt);
            figure();
            subplot(3,2,1);
            plot(1:T-1,survivalProb,'k-','DisplayName','Matlab');
            hold('on');
            plot(1:T-1,survivalProb_C,'r-','DisplayName','C');
            title('Survival Probability: Q(r,t)');
            legend('location','best');

            subplot(3,2,2);
            plot(1:T-1,abs(survivalProb-survivalProb_C),'bo','MarkerFaceColor',[0,0,1]);
            hold('on');
            set(gca(),'YScale','log');
            title('Q-Delta');
            
            subplot(3,2,3);
            plot(1:T-1,mu,'k-','DisplayName','Matlab');
            hold('on');
            plot(1:T-1,mu_C,'r-','DisplayName','C');
            title('Capture Target Conditional Prob: mu(r,t)');
            legend('location','best');

            subplot(3,2,4);
            plot(1:T-1,abs(mu-mu_C),'bo','MarkerFaceColor',[0,0,1]);
            hold('on');
            set(gca(),'YScale','log');
            title('mu-Delta');

            subplot(3,2,5);
            plot(1:T-1,nu,'k-','DisplayName','Matlab');
            hold('on');
            plot(1:T-1,nu_C,'r-','DisplayName','C');
            title('Capture-to-Capture Target Conditional Prob: nu(t)');
            legend('location','best');

            subplot(3,2,6);
            plot(1:T-1,abs(nu-nu_C),'bo','MarkerFaceColor',[0,0,1]);
            hold('on');
            set(gca(),'YScale','log');
            title('nu-Delta');
        end


        function [llhSum, llh] = computeSingleParticleLLH_Distance(obs,p,params,prior, plotOn)
            %For now ignore the prior.  We'll come back to it
            RD = RDCapture(params.DA+params.DB, params.rho_bind, params.lambda_bind);
            OU = OUProcess(params.rho, params.sigma_rho, params.gamma, params.rho);
            T = size(obs,1); %Number of times
            times = obs(:,1);
            dt = diff(times);
            llh = zeros(T,9);
            ys = p(:,2:5);
            vs = PairInteractionMCMCModel.transV(ys); %Transform to V-space
            rs = vs(:,3);
            xs = obs(:,2:5);
            states = p(:,1);
            
            [~,llh(1,1:5)] = PairInteractionMCMCModel.llhG0(p(1,1),p(1,2:5),params,prior);

            %XVariance
            obs_ds = ys - xs;
            obs_vars = obs(:,6:9).^2;
            llh(:,6:9) = PairInteractionMCMCModel.gaussianLLH(obs_ds,obs_vars); %Observation variances

            %Free Model Ys
            free_ds = diff(ys); %Free distances 
            free_vars = repmat(dt,1,4).*repmat(2*[params.DA, params.DA, params.DB, params.DB],T-1,1); %Free variances 2*DA*dt or 2*DB*dt
            llh_free = PairInteractionMCMCModel.gaussianLLH(free_ds,free_vars);

            %Bound Model Vs = (c_x,c_y,r,phi)
            %Bound center
            llh_bound = zeros(T-1,4);
            bound_C_ds = diff(vs(:,1:2));
            bound_C_vars = repmat(dt,1,2).*repmat(2*params.DAB,T-1,2); %2*DAB*dt
            llh_bound(:,1:2) = PairInteractionMCMCModel.gaussianLLH(bound_C_ds, bound_C_vars);
            
            %Bound distance r
            rLLH = OU.sampleLLH(times,rs);           
            llh_bound(:,3) = rLLH(2:end) - log(rs(2:end)); %OU and log determinant of V^-1 [log(1/r)] 
            
            %Bound angle phi
            kappa = 1./(2*params.Dphi*dt);
            llh_bound(:,4) = PairInteractionMCMCModel.vonMisesLLH(diff(vs(:,4)), kappa); % VonMises

            %Zs
            %Koff
            expCDF = exp(-params.koff*dt);
            llh_BB = log(expCDF);
            llh_BF = log(1-expCDF);
            
            %Kon
            survivalProb = RD.computeSurvivalProb(rs(1:end-1),dt);
%             PairInteractionMCMCModel.debugRDCaptureSurvivalProb(obs,p,params)

            llh_FF = log(survivalProb);
            llh_FB = log(1-survivalProb);

            %Prior for state Z
            % if states(1) %Bound
            %    llh(1,1) = log(prior.state(2));
            % else %Free
            %    llh(1,1) = log(prior.state(1));                
            % 
%             FtoB_R_LLH = PairInteractionMCMCModel.gaussianLLH(0, params.sigma_rho^2);

            for t=2:T
                if states(t) % Bound at t
                    if states(t-1) % Bound at t-1
                        llh(t,1)=llh_BB(t-1); %B->B
                        llh(t,2:5) = llh_bound(t-1,:);
                    else
                        llh(t,1)=llh_FB(t-1); %F->B
                        %Initial conditions are a little different for F->B transition since we have no
                        %prior for r or phi given z_t=bound and z_t-1 = free
%                         llh(t,2:3) = PairInteractionMCMCModel.gaussianLLH(vs(t,1:2)-vs(t-1,1:2),2*dt(t-1)*(params.DA+params.DB)*[1,1]);
                        llh(t,2:3) = llh_bound(t-1,1:2);
                        FtoB_R_LLH = PairInteractionMCMCModel.gaussianLLH(abs(rs(t)-params.rho), params.sigma_rho^2);
                        llh(t,4) = FtoB_R_LLH - log(rs(t)); %Equlibrium for distance r
                        center = vs(t-1,1:2);
                        cholC = eye(2)*sqrt(2*dt(t-1)*(params.DA+params.DB));
                        [muPhi, kappaPhi] = PairInteractionMCMCModel.estimateVonMisesKappa(center, cholC,8,sqrt(2));
                        llh(t,5) = PairInteractionMCMCModel.vonMisesLLH(vs(t,4)-muPhi,kappaPhi);
                    end
                else % Free at time t
                    if states(t-1)
                        llh(t,1)=llh_BF(t-1); %B->F
                        llh(t,2:5) = llh_free(t-1,:);
                    else
                        llh(t,1)=llh_FF(t-1); % F-F
                        llh(t,2:5) = llh_free(t-1,:);
                    end
                end
%                 llhz = PairInteractionMCMCModel.propogatorLLH(ys(t,:),dt(t-1),states(t-1),ys(t-1,:),params);
%                 llh1 = llhz(states(t)+1);
%                 llh2 = sum(llh(t,2:5));
% 
%                 assert(abs(llh1-llh2)>eps);
            end

            llhSum = sum(llh(:));

            if nargin==5 && plotOn
                figure();
                subplot(3,2,1);
                plot(times,vs(:,3),'r-','DisplayName','Dist: r');
                hold();
                yl=ylim();
                for t=2:T
                    if states(t)
                        xs=times([t-1,t,t,t-1])+dt(t-1)/2;
                        ys=[0,0, yl(2), yl(2)];
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                end
                title('Distance');

                legend('location','best');
                xlabel('time (s)');
                subplot(3,2,2);
                plot(times,cos(vs(:,4)),'m-','DisplayName','cos(Phi)');
                legend('location','best');
                xlabel('time (s)');
                hold();
                yl=ylim();
                for t=2:T
                    if states(t)
                        xs=times([t-1,t,t,t-1])+dt(t-1)/2;
                        ys=yl([1,1,2,2]);
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                end
                title('Relative angle Phi');

                subplot(3,2,3);
                plot(times,llh(:,1),'k-','DisplayName','LLH:Z');
                hold();
                plot(times,sum(llh(:,2:5),2),'r-','DisplayName','LLH:Y');
                plot(times,sum(llh(:,6:9),2),'b-','DisplayName','LLH:X');
                legend('location','best');
                xlabel('time (s)');
                yl=ylim();
                for t=2:T
                    if states(t)
                        xs=times([t-1,t,t,t-1])+dt(t-1)/2;
                        ys=yl([1,1,2,2]);
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                end
                title('Log-likelihood');

                subplot(3,2,4);
                plot(times(2:end),sum(llh_free,2),'r-','DisplayName','LLH:Y free');
                hold();
                plot(times(2:end),sum(llh_bound,2),'b-','DisplayName','LLH:Y bound');
                plot(times(2:end),sum(llh_bound(:,1:2),2),'g-','DisplayName','LLH:Y_C bound');
                plot(times(2:end),llh_bound(:,3),'m-','DisplayName','LLH:Y_r bound');
                plot(times(2:end),llh_bound(:,4),'-','Color',[1,.3,0],'DisplayName','LLH:Y_phi bound');
                ylim([-15,10]);
                legend('location','best');
                xlabel('time (s)');
                title('LLH Y: Free vs Bound');

                subplot(3,2,5);
                plot(times(2:end),llh_FF,'r-','DisplayName','LLH:Z_F->Z_F');
                hold();
                plot(times(2:end),llh_FB,'r--','DisplayName','LLH:Z_F->Z_B');
                plot(times(2:end),llh_BB,'b-','DisplayName','LLH:Z_B->Z_B');
                plot(times(2:end),llh_BF,'b--','DisplayName','LLH:Z_B->Z_F');
                legend('location','best');
                xlabel('time (s)');
                title('LLH Binding and Unbinding');
                ylim([-25,1]);
                yl=ylim();
                for t=2:T
                    if states(t)
                        xs=times([t-1,t,t,t-1])+dt(t-1)/2;
                        ys=yl([1,1,2,2]);
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                end
            end
        end
        
        function Obs = convertPosteriorSamplesToObservation(samples)
            % [in]
            %  samples (Nsteps,17,Nsamples) array produced from posterior sampling
            % [out]
            %  Obs - cell array of observation vectors [t ax ay bx by SEax SEay SEbx SEby trueState]
            Nsamples = size(samples,3);
            Obs = cell(Nsamples,1);
            for n=1:Nsamples
                Obs{n} = samples(:,[1,10:17,2]);                
            end
        end        
        

        function p_rs = pairR_PDF_rect(rs,w,h)
            %  p(|A-B|=r)  where A,B~Uniform(Rect(l,b,w,h))
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
            %  p(|A-B|=r)  where A,B~Uniform(Circle(ceter,a))
            % [in]
            %   r - vector Nx1 of distances to evaluate pdf at
            %   a - radius of circle
            % [out]
            %  p_r - vector Nx1 of pdf vals for each r value
            theta = 2*acos(r/(2*a));
            p_r = (2*r).*(theta - sin(theta))/(pi*a^2);
        end
       
        function p_r = pairR_PDF_2Dgauss(r,Sigma)
            %  p(|A-B|=r)  where A,B~Gauss(ceter,Sigma)
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
            q4=q2*q2;
            p_r = (1/(q*omega))*(1+q2)*r.*exp(-(1+q2)^2/(4*q2*omega)*r.^2).*besseli(0,(1-q4)/(4*q2*omega)*r.^2);
        end

        function llhs=multiGaussianLLH(ds,Sigma)
            % Compute LLH for a sequence of k dimensional samples from a single
            % [in]
            % ds: size:[k,N] each row is a distance vector in R^k
            % SigmaL size:[k,k]
            % [out]
            % llh: vector [1,N] Log-likelihood for each indiviudal ds
            % llh_all: [optional] size:[k,N] Log-likelihoods along each dimension
            llhs=-.5*(size(ds,1)*PairInteractionMCMCModel.LOG2PI + log(det(Sigma)) + sum(ds.*(Sigma\ds)));
        end

        function llhs=gaussianLLH(ds,vars)
            % Compute LLH for a vector or matrix of independet 1D gaussians with distances: ds and sigma: sigmas
            % [in]
            %   ds - size:[M,N] distances
            %   vars - size:[M,N] variances
            % [out]
            %   llhs: size:[M,N] Log-likelihoods for each independent 1D gaussian
            llhs = -.5*(PairInteractionMCMCModel.LOG2PI + log(vars) + ds.^2./vars);
        end

        function llhs=vonMisesLLH(ds,kappa)
            % log-liklihood for vonMises distribution for diffusion in angle
            % [in]
            %   ds - size:[N] vector independent of angular displacements to draw from the same vonMises
            %           distribution
            %   kappa - angular distribution parameter Use kappa=1/(2*Dphi*dt) for the vonMises spread
            %   parameter for a diffusing angle
            % [out]
            %   llhs = size:[N] log-likelihood of vonMises distribution 
            %             kappa = 1./(2*D); %VonMises variance
            llhs = kappa(:).*cos(ds(:)) - log(2*pi) - PairInteractionMCMCModel.logBesselI0(kappa(:)); % VonMises
        end

        function z = logBesselI0(kappa)
            z = zeros(size(kappa));
            for n=1:numel(kappa)
                if kappa(n) >= 3.75
                    y = 3.75 / kappa(n);
                    v1 = (0.01328592 + (0.00225319 + (-0.00157565 + (0.00916281 + (-0.02057706 + ...
                         (0.02635537 + (0.00392377 * y - 0.01647633) * y) * y) * y) * y) * y) * y) * y;
                    z(n) = kappa(n) - 0.5 * log(kappa(n)) + log(v1 + 0.39894228);
                elseif kappa(n) == 0.0
                    z(n) = 0.0;
                else
                    y = (kappa(n)/3.75)^2;
                    v2 = (3.5156229 + (3.0899424 + (1.2067492 + (0.2659732 + (0.0360768 + ...
                          0.0045813 * y) * y) * y) * y) * y) * y;
                    z(n) = log(1.0 + v2);
                end
            end
        end

        function [angles,llhs] = sampleVonMises(mu, kappa,N)
            
            % if kappa is small, treat as uniform distribution
            if kappa < 1e-6
                angles = 2*pi*rand(1,N);
                llhs = repmat(-PairInteractionMCMCModel.LOG2PI,1,N);
            else
                % other cases
                a = 1 + sqrt(1+4*kappa.^2);
                b = (a - sqrt(2*a))/(2*kappa);
                r = (1 + b^2)/(2*b);

                angles = zeros(1,N);
                for j = 1:N
                  while true
                      z = cos(pi*rand());
                      f = (1+r*z)/(r+z);
                      c = kappa*(r-f);
                      u = rand();
                      if u < c*(2-c) || 0 <= log(c)-log(u)+1-c
                         break
                      end
                  end
                  angles(j) = mod(angle(exp(1i*(mu +  sign(rand() - 0.5) * acos(f)))),2*pi);
                end
                llhs = PairInteractionMCMCModel.vonMisesLLH(mu-angles,kappa);
            end
        end

        function [mu, kappa] = estimateVonMisesKappa(center, cholC, N, alpha)
            % Estimate the kappa parameter for the vonMises distribution from a sample of
            % angles
            % [in]
            %  center - size:[2,1] center positions of Gaussian 
            %  cholC - size:[2,2] cholesky decomposition of covariance matrix
            %  N - [optional] number of samples >4 [default:8]
            %  alpha - [optional] radius out from center to use [default: sqrt(2)]
            % [out]
            %  mu - center point for the vonMises distribution (the angle of m) for convenience
            %  kappa - estimate of the spread of kappa
            if nargin<4
                alpha = sqrt(2);
            end
            if nargin<3
                N=8;
            elseif N<4
                error('PairInteractionMCMCModel:estimateVonMisesKappa','Require N>=4 points for reasonable accuracy');
            end
            samp_angles = linspace(0,2*pi*(N-1)/N,N);
            samp_points = alpha*[cos(samp_angles);sin(samp_angles)];
            samp = repmat(center(:),1,N) + cholC*samp_points;
            observedAngles = atan2(samp(2,:),samp(1,:));
            kappa = PairInteractionMCMCModel.estimateVonMisesKappaSample(observedAngles);
            mu = atan2(center(2),center(1));
        end

        function kappa = estimateVonMisesKappaSample(observedAngles)
            % Estimate the kappa parameter for the vonMises distribution from a sample of
            % angles
            % [in]
            %  observedAngles - vector of angles that were observed
            % [out]
            %  kappa - estimate of the spread of kappa
            N=numel(observedAngles);
            r = abs(sum(exp(1j*observedAngles)))/N;
            if r < 0.53
                kappa = 2*r + r^3 +(5/6)*r^5;
            elseif r < 0.85
                kappa = -.4 + 1.39*r + 0.43/(1-r);
            else
                kappa = 1/(r^3-4*r^2+3*r);
            end
            if N<15
                kappa = max(kappa-2/(N*kappa),0);    
            else
                kappa = (N-1)^3*kappa/(N^3+N);
            end
        end

        function vs = transV(ys)
            % Implement the Y->V transform
            % [in]
            %  ys - [Nx4xK] cols [ax,ay,bx,by]
            % [out]
            %  vs - [Nx4xK] cols [cx,cy,r,theta]
            assert(size(ys,2)==4);
            vs = [ .5*(ys(:,1:2,:) + ys(:,3:4,:)), ...
                   sqrt(sum((ys(:,1:2,:) - ys(:,3:4,:)).^2,2)), ...
                   mod(atan2(ys(:,2,:) - ys(:,4,:),ys(:,1,:) - ys(:,3,:))+2*pi,2*pi)];
        end

        function rs = transVr(ys)
            % Implement the Y->Vr transform
            % [in]
            %  ys - [Nx4xK] cols [ax,ay,bx,by]
            % [out]
            %  rs - [NxK] 
            assert(size(ys,2)==4);
            rs = sqrt(sum((ys(:,1:2,:) - ys(:,3:4,:)).^2,2));
        end

        function phis = transVphi(ys)
            % Implement the Y->Vphi transform
            % [in]
            %  ys - [Nx4xK] cols [ax,ay,bx,by]
            % [out]
            %  phis - [NxK] 
            assert(size(ys,2)==4);
            phis = mod(atan2(ys(:,2,:) - ys(:,4,:),ys(:,1,:) - ys(:,3,:)),2*pi);
        end

        function ys = transVinv(vs)
            % Probably flip this later
            % Implement the V->Y transform
            % [in]
            %  vs - [Nx4xK] cols [cx,cy,r,theta]
            % [out]
            %  ys - [Nx4xK] cols [ax,ay,bx,by]
            if numel(vs)==4
                ys=zeros(size(vs));
                cosphi = cos(vs(4));
                sinphi = sin(vs(4));
                d = .5*vs(3);
                ys(1) = vs(1)+d*cosphi;
                ys(2) = vs(2)+d*sinphi;
                ys(3) = vs(1)-d*cosphi;
                ys(4) = vs(2)-d*sinphi;
            else
                assert(size(vs,2)==4);
                sinphi=sin(vs(:,4,:));
                cosphi=cos(vs(:,4,:));
                ys = [vs(:,1,:) + .5*vs(:,3,:).*cosphi,...
                      vs(:,2,:) + .5*vs(:,3,:).*sinphi,...
                      vs(:,1,:) - .5*vs(:,3,:).*cosphi,...
                      vs(:,2,:) - .5*vs(:,3,:).*sinphi];
            end
        end


        function ys = transVinv_row(vs)
            % Implement the V->Y transform where vs and ys store components in rows
            % [in]
            %  vs - [4xN] rows [cx,cy,r,theta]
            % [out]
            %  ys - [4xN] rows [ax,ay,bx,by]
            assert(size(vs,1)==4);
            sinphi=sin(vs(4,:));
            cosphi=cos(vs(4,:));
            ys = [vs(1,:) + .5*vs(3,:).*cosphi;...
                  vs(2,:) + .5*vs(3,:).*sinphi;...
                  vs(1,:) - .5*vs(3,:).*cosphi;...
                  vs(2,:) - .5*vs(3,:).*sinphi];
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

        function [new_mean, new_cov] = gaussianBayes(yt_mean, yt_var, xt_mean, xt_var)
            xt_inv = inv(diag(xt_var));
            yt_inv = inv(yt_var);
            new_cov_inv = yt_inv + xt_inv;
            new_cov = inv(new_cov_inv);
            new_mean = new_cov_inv\(diag(xt_var)\xt_mean(:)+yt_var\yt_mean);
        end

        function plotSE(mu,sigma,color,width,label)
            alphas=linspace(0,2*pi,100);
            ps1=[mu(1)+sigma(1)*cos(alphas); mu(2)+sigma(2)*sin(alphas)];
            ps3=[mu(1)+3*sigma(1)*cos(alphas); mu(2)+3*sigma(2)*sin(alphas)];
            plot(ps1(1,:),ps1(2,:),'-','Color',color,'LineWidth',width,'DisplayName',sprintf('%s 1-sigma',label));
            plot(ps3(1,:),ps3(2,:),'-','Color',color,'LineWidth',width/3,'DisplayName',sprintf('%s 3-sigma',label));
        end
    end

    methods (Access=protected)
        function id=addData(obj, D)
            obj.Data = [obj.Data; obj.DefaultData];
            names = fieldnames(obj.DefaultData);
            for k = 1:numel(names)
                name=names{k};
                if isfield(D,name)
                    obj.Data(end).(name) = D.(name);
                end
            end
            id = numel(obj.Data);
        end

        function id = addExperiment(obj, E)
            obj.Experiment = [obj.Experiment; obj.DefaultExperiment];
            names = fieldnames(obj.DefaultExperiment);
            for k = 1:numel(names)
                name=names{k};
                if isfield(E,name)
                    obj.Experiment(end).(name) = E.(name);
                end
            end
            id = numel(obj.Experiment);
        end

        function id = addResults(obj, R)
            obj.Results = [obj.Results; obj.DefaultResults];
            names = fieldnames(obj.DefaultResults);
            for k = 1:numel(names)
                name=names{k};
                if isfield(R,name)
                    obj.Results(end).(name) = R.(name);
                end
            end
            id = numel(obj.Results);
        end

        function estimatePositionPrior(obj, exp_id, posPriorShape)
            % Compute and set the prior for given experiment using all data already availible.
            % [in]
            %   exp_id: index of single experiment. 
            %   posPriorShape: one of {'Gaussian','Circular','Rectangular'}
            % 
            if nargin<3
                posPriorShape='Gaussian';
            end
            
            % Find initial positions
            filter = [obj.Data(:).experimentId] == exp_id;
            allPos = cellmatfun(@(O) max(O(:,2:5)+O(:,6:9),O(:,2:5)-O(:,6:9)), {obj.Data(filter).obs}'); %allPos has cols [ax ay bx by]
            allPos = [allPos(:,1:2);allPos(:,3:4)]; % now cols are [xpos ypos]
            
            mu = mean(allPos);
            displacement = allPos-repmat(mu,size(allPos,1),1);
            % Initialize prior
            prior = obj.DefaultPrior;
            prior.posShape = posPriorShape;
            switch lower(posPriorShape)
                case {'gaussian','gauss'}
                    prior.posPoint = mu; % center
                    prior.posExtent = cov(allPos); % 2x2 covariance matrix
                case {'circular','circ'}
                    prior.posPoint = mu; % center
                    prior.posExtent = max(sqrt(sum(displacement^2,2))); % radius
                case {'rectangular','rect'}
                    prior.posPoint = min(allPos); % [left, bottom]
                    prior.posExtent = max(allPos)-prior.posPoint; %[w, h]
            end
            obj.Experiment(exp_id).prior = prior;
        end
    end
    
    methods
        function ndata = get.Ndata(obj)
            ndata=numel(obj.Data);
        end
    end
end

