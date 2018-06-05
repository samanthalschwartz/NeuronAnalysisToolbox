classdef PairInteractionMarkovModel < handle
    % A Hidden Markov model based dimer detection class.

    % This simplifies the Bayesian network to allow representing it as an auto-regresive HMM
    % on which efficient exact inference is possible with forwards-backwards and other HMM algorithms
    % These are typically linear time.

    % Outputs:
    % obj.statesMAP: Size:[n,1] The single most likely a-posteriori state sequence
    % obj.statesMPM: Size:[n,1] The maximum of the posterior marginals for each time.
    % obj.samplePosterior(nSamples): Size:[n,nSamples] A sample set of sequences from the posterior distributions=
    
    properties (Constant=true)
        MODELS={'Distance','Diffusive'};
        DEFAULT_RHO = 0.050; %microns - bond distance
        DEFAULT_SIGMA_REG = 0.005; %microns - 

        DefaultParams=struct('Distance', struct(...
                                        'D1',0.050,... %(um^2/s) monomer diffusion constant 
                                        'D2',0.050,... %(um^2/s) dimer diffusion constant 
                                        'rho',0.050,... %(um) dimer bond distance between fluorophores
                                        'kon',1e-1,... %(1/s) (on) F->B rate (Fudge)
                                        'koff',1e-1,...%(1/s) (off) B->F rate.  (True kinetic rate)
                                        'sigmaReg',0.010... %(um) The standard deviation of the registration error
                                        ),...
                             'Diffusive', struct(...
                                        'DA',0.050,... %(um^2/s) particle A diffusion constant 
                                        'DB',0.050,... %(um^2/s) particle B diffusion constant 
                                        'DAB',0.050,... %(um^2/s) dimer diffusion constant 
                                        'rho',0.050,... %(um) dimer bond distance between fluorophores
                                        'kon',1e-1,... %(1/s) (on) F->B rate (Fudge)
                                        'koff',1e-1,...%(1/s) (off) B->F rate.  (True kinetic rate)
                                        'sigmaReg',0.010... %(um) The standard deviation of the registration error
                                        ));

    end
    properties
        %User supplied data
        data; %cell array of pair-format matricies.  1 matrix for each tracked pair.
              % pair matrix column format is [t x1 y1 x2 y2 sx1 sy1 sx2 sy2].
        statesTrue; % A cell-array of known true states assignments, if availible. (i.e., from simulation)
                    % leave empty if data is from a real experimental source.
        Ndata; %number of data items
        Nsteps; % size:[Ndata,1]  Number of steps for each data point

        %User supplied params
        modelName;
        modelParams;
        deltaT; %Time step in s
        prior; %Prior on states [F B];  Must sum to 1
        
        %Computed params
        evidence; % size:[n,2], cols:[F_i, B_i] - Evidence for free at time i.  This is what Bayes rules gives us.
        transitions; % transition probabilities
        fwdFilt; % size:[n,2], cols:[alpha_i_F, alpha_i_B]  
        bkwFilt; % 
        marginal; %
        marginal2Slice; %
        statesMAP; %
        statesMPM; %Maximum of posterior marginals        
    end

    methods
        function obj = PairInteractionMarkovModel(deltaT, varargin)
            obj.deltaT = deltaT;
            if nargin>0
                obj.loadData(varargin{:})
            end
        end
        
        function loadData(obj,newData)
            % data is a cellarray of matrices
            % column format is [t x1 y1 x2 y2 sx1 sy1 sx2 sy2 trueState].
            % true state is optional so each matrix should have 9 or 10 cols.
            % all distance units are in microns, times are in seconds.
            % trueState==0 -- Free
            % trueState==1 -- Bound
            newData = makecell(newData);
            if ~all(cellfun(@(D) size(D,2)>=9,newData))
                error('PairInteractionModel:loadData','Incorrect data matricies format');
            end
            obj.data = [obj.data; newData];
            if size(newData{1},2)>=10
                obj.statesTrue = [obj.statesTrue; cellmap(@(D) D(:,10), newData)];
            end
            obj.Ndata = numel(obj.data);
            obj.Nsteps = cellfun(@(O) size(O,1), obj.data);
        end

        function initializeModelParams(obj, model, newParams)
            if nargin<3
                newParams = obj.DefaultParams.(model);
            else
                switch model
                    case 'Distance'
                        newParams = obj.makeModelParams_Distance(newParams);
                    case 'Diffusive'
                        newParams = obj.makeModelParams_Diffusive(newParams);
                    otherwise
                        error('PairInteractionMarkovModel:initializeModelParams','Unknown model %s', model);
                end
            end
            obj.modelName = model;
            obj.modelParams.(obj.modelName) = newParams;
        end

        function calculateModel(obj, name, params, prior)
            if nargin==1 && isempty(obj.modelName)
                error('PairInteractionMarkovModel:calculateModel','No model name set');
            end
            if isempty(obj.prior)
                obj.prior = [.5 .5];
            end
            if nargin>=2 
                obj.modelName = name;
            end
            if nargin>=3
                obj.params = params;
            end
            if nargin>=4
                obj.prior = prior;
            end
            switch obj.modelName
                case {'Distance','Diffusive'}
                    obj.calculateModelHMM();
                otherwise
                    error('PairInteractionMarkovModel:calculateModel','Unknown model %s', model);
            end
        end
        
        function plotDistance(obj,idx)
            figure('Color',[1,1,1],'Position',[0,0, 600, 400]);
            T = obj.data{idx};
            N = size(T,1);
            P = zeros(N,18);
            P(:,[1,2,3,6,7,10,11,12,13]) = T(:,1:9);
            certainty=0.95;
            D = PairAnalysis.makePairDists(P,certainty);
            plot(D(:,1),D(:,2),'-r','LineWidth',2,'DisplayName','Dist');
            hold on;
            plot([0, T(end,1)],[0.05 0.05], 'k:','DisplayName','Binding Radius');
            xs = [D(:,1)', fliplr(D(:,1)')];
            ys = [D(:,3)', fliplr(D(:,4)')];
            xlim([T(1,1), T(end,1)]);
%             fill(xs,ys,'r','EdgeColor','none','FaceAlpha',0.4,'DisplayName',sprintf('Dist (%.1f%% conf.)',certainty*100));
            xlabel('Time (s)');
            ylabel('Dist (um)');
            legend('Location','best');
        end
        
        function plotHMMMarginal(obj, idx)
            T = obj.data{idx};
            N = size(T,1);
            E = obj.evidence{idx};
            figure('Color',[1,1,1],'Position',[0,0, 600, 400]);
            fwdFilter = obj.fwdFilt{idx};
            bkwFilter = obj.bkwFilt{idx};
            marg = obj.marginal{idx};
            pB1 = E(:,2);
            pBf = fwdFilter(:,2);
            pBb = bkwFilter(:,2);
            pBm = marg(:,2);
            plot(1:N, pB1,'k-','DisplayName', 'Immediate evidence for B.  $P(z_t=B | x_t)$');
            hold;
            plot(1:N, pBf,'r-','DisplayName', 'Forwards filtered evidence for B.  $P(z_t=B | x_{1:t})$');
            plot(1:N, pBb,'b-','DisplayName', 'Backwards filtered evidence for B. $P(x_{t+1:t}| z_t=B)$');
            plot(1:N, pBm,'g-','DisplayName', 'Marginal (smoothed) evidence for B.  $P(z_t=B | x_{1:T})$');
            ylim([0,1]);
            L=legend('location','best');
            L.Interpreter='latex';
%             title('Fwd-Backwd Marginal Evidence --- Distance Model (Low-Nam, et. al.)');
            xlabel('Observation index');
            ylabel('Probability of Bound');
        end

        function plotHMMEvidence(obj, idx)
            T = obj.data{idx};
            N = size(T,1);
            figure('Color',[1,1,1],'Position',[0,0, 600, 350]);
            E = obj.evidence{idx};
            plot(1:N,log(E(:,1)),'r-','DisplayName','Free evidence $P\left(z_t=F \mid x_t, x_{t-1}\right)$');
            hold();
            plot(1:N,log(E(:,2)),'b-','DisplayName','Bound evidence $P\left(z_t=B \mid x_t\right)$');
            xlim([1,N]);
            yl = ylim();
            if ~isempty(obj.statesTrue)
                stateSeq = obj.statesTrue{idx};
                state = stateSeq(1);
                lastIdx = 1;
                for i=1:size(T,1)
                    if state ~= stateSeq(i)
                        if state
                            xs=[lastIdx,i+1, i+1, lastIdx];
                            ys=[yl(1),yl(1), yl(2), yl(2)];
                            h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                            h.Annotation.LegendInformation.IconDisplayStyle='off';
                        end
                        lastIdx=i;
                        state = stateSeq(i);
                    end
                end
                ylim(yl);
            end
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize = 12;
            title('Log-Evidence (Normalized)');
            xlabel('Observation index');
            ylabel('Log-Probability');
        end
        
        function plotHMMForwardsEvidence(obj, idx)
            T = obj.data{idx};
            N = size(T,1);
            figure('Color',[1,1,1],'Position',[0,0, 600, 350]);
            E = obj.evidence{idx};
            
            plot(1:N, E(:,2),'k-','DisplayName', 'Immediate evidence for B.  $P(z_t=B \mid \mathbf{x}_t)$');
            hold;
            plot(1:N, obj.fwdFilt{idx}(:,2),'-','Color',[0.8,0,0.7],'LineWidth',2,'DisplayName', 'Forwards filtered evidence for B.  $P(z_t=B \mid \mathbf{x}_{1:t})$');
            xlim([1,N]);

            yl = ylim();
            stateSeq = obj.statesTrue{idx};
            state = stateSeq(1);
            lastIdx = 1;
            for i=1:size(T,1)
                if state ~= stateSeq(i)
                    if state
                        xs=[lastIdx,i+1, i+1, lastIdx];
                        ys=[yl(1),yl(1), yl(2), yl(2)];
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                    lastIdx=i;
                    state = stateSeq(i);
                end
            end
            ylim(yl);
            
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize = 12;
            xlabel('Observation index');
            ylabel('Probability of Bound');
        end

        function plotHMMBackwardsEvidence(obj, idx)
            T = obj.data{idx};
            N = size(T,1);
            figure('Color',[1,1,1],'Position',[0,0, 600, 350]);
            E = obj.evidence{idx};
            
            plot(1:N, E(:,2),'k-','DisplayName', 'Immediate evidence for B.  $P(z_t=B \mid \mathbf{x}_t)$');
            hold;
            plot(1:N, obj.bkwFilt{idx}(:,2),'-','Color',[1, 0.5, 0],'LineWidth',2,'DisplayName', 'Backwards filtered evidence for B.  $P(\mathbf{x}_{t+1:T}\mid z_t=B)$');
            xlim([1,N]);

            yl = ylim();
            stateSeq = obj.statesTrue{idx};
            state = stateSeq(1);
            lastIdx = 1;
            for i=1:size(T,1)
                if state ~= stateSeq(i)
                    if state
                        xs=[lastIdx,i+1, i+1, lastIdx];
                        ys=[yl(1),yl(1), yl(2), yl(2)];
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                    lastIdx=i;
                    state = stateSeq(i);
                end
            end
            ylim(yl);
            
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize = 12;
            xlabel('Observation index');
            ylabel('Probability of Bound');
        end

        
        function plotHMMMarginalEvidence(obj, idx)
            T = obj.data{idx};
            N = size(T,1);
            figure('Color',[1,1,1],'Position',[0,0, 600, 350]);
            E = obj.evidence{idx};
            
            plot(1:N, E(:,2),'k-','DisplayName', 'Immediate evidence for B.  $P(z_t=B \mid \mathbf{x}_t)$');
            hold;
            plot(1:N, obj.fwdFilt{idx}(:,2),'-','Color',[0.8,0,0.7],'LineWidth',2,'DisplayName', 'Forwards filtered evidence for B.  $P(z_t=B \mid \mathbf{x}_{1:t})$');
            plot(1:N, obj.bkwFilt{idx}(:,2),'-','Color',[1, 0.5, 0],'LineWidth',2,'DisplayName', 'Backwards filtered evidence for B.  $P(\mathbf{x}_{t+1:T}\mid z_t=B)$');
            plot(1:N, obj.marginal{idx}(:,2),'-','Color',[0, 1, 0],'LineWidth',2,'DisplayName', 'Marginal (smoothed) evidence for B.  $P(z_t=B \mid \mathbf{x}_{1:T})$');
            xlim([1,N]);

            if ~isempty(obj.statesTrue)
                yl = [0,1];
                stateSeq = obj.statesTrue{idx};
                state = stateSeq(1);
                lastIdx = 1;
                for i=1:size(T,1)
                    if state ~= stateSeq(i)
                        if state
                            xs=[lastIdx,i+1, i+1, lastIdx];
                            ys=[yl(1),yl(1), yl(2), yl(2)];
                            h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                            h.Annotation.LegendInformation.IconDisplayStyle='off';
                        end
                        lastIdx=i;
                        state = stateSeq(i);
                    end
                end
                ylim(yl);
            end
            
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize = 12;
            xlabel('Observation index');
            ylabel('Probability of Bound');
        end
        
        function plotHMMInferredStates(obj, idx)
            Traj = obj.data{idx};
            N = size(Traj,1);
            figure('Color',[1,1,1],'Position',[0,0, 600, 350]);
            E = obj.evidence{idx};
            d = sqrt((Traj(:,2)-Traj(:,4)).^2 + (Traj(:,3)-Traj(:,5)).^2); %observed distance
            subplot(2,1,1);
            plot(1:N,d,'-','Color',[1,0,0],'LineWidth',1,'DisplayName', 'Observed distance');
            hold();
            plot(1:N, obj.marginal{idx}(:,2),'-','Color',[0, 1, 0],'LineWidth',2,'DisplayName', 'Marginal (smoothed) evidence for B.  $P(z_t=B \mid \mathbf{x}_{1:T})$');
            plot([1,N],[.5 .5],'k:','LineWidth',0.75,'DisplayName','Marginal Likelihood 50\% threshold');
            yl = ylim();
            xlim([1,N]);
            if ~isempty(obj.statesTrue)
                stateSeq = obj.statesTrue{idx};
                state = stateSeq(1);
                lastIdx = 1;
                for i=1:size(Traj,1)
                    if state ~= stateSeq(i)
                        if state
                            xs=[lastIdx,i+1, i+1, lastIdx];
                            ys=[yl(1),yl(1), yl(2), yl(2)];
                            h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                            h.Annotation.LegendInformation.IconDisplayStyle='off';
                        end
                        lastIdx=i;
                        state = stateSeq(i);
                    end
                end
                ylim(yl);
            end
            
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize=12;
            xlabel('Observation index');
            ylabel('Distance (um)');
            
            subplot(2,1,2);
            T = obj.transitions;
            llhMPE = PairInteractionMarkovModel.logLikelihood_HMM(obj.statesMAP{idx},obj.prior,E,T);
            llhMPM = PairInteractionMarkovModel.logLikelihood_HMM(obj.statesMPM{idx},obj.prior,E,T);
            stairs(1:N, obj.statesMAP{idx}+0.04, 'Color',[1,0.0,0.0],'Linewidth',3,'DisplayName',...
                    sprintf('[LLH = %.6g] MAP (Max \\emph{a Posteriori}) (Viterbi). $\\mathrm{arg}\\max_{z_{1:T}} P\\left( z_{1:T} \\mid \\mathbf{x}_{1:t} \\right)$',llhMPE));
            hold();
            xlim([1,N]);
            stairs(1:N, obj.statesMPM{idx}+0.08, 'Color',[0,1,0.5],'Linewidth',3,'DisplayName',...
                sprintf('[LLH = %.6g] MPM (Max posterior marginal).',llhMPM));
            if ~isempty(obj.statesTrue)
                trueStateSeq = obj.statesTrue{idx};
                llhTru = PairInteractionMarkovModel.logLikelihood_HMM(trueStateSeq,obj.prior,E,T);
                stairs(1:N, trueStateSeq, 'Color','k','Linewidth',1,'DisplayName',...
                    sprintf('[LLH = %.6g] True state',llhTru));
            end
            
            set(gca(),'YTick',[0,1]);
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize = 12;
            ylim([0,2]);
            xlabel('Observation index');
            ylabel('State (shifted to see)');            
        end
        
        
        function plotHMMPosteriorSample(obj, idx)
            Traj = obj.data{idx};
            N = size(Traj,1);
            figure('Position',[0,0, 600, 350]);
            E = obj.evidence{idx};
            d = sqrt((Traj(:,2)-Traj(:,4)).^2 + (Traj(:,3)-Traj(:,5)).^2); %observed distance
            subplot(2,1,1);
            plot(1:N,d,'k-','DisplayName', 'Observed distance');
            hold();
            plot(1:N, obj.marginal{idx}(:,2),'-','Color',[0, 1, 0],'LineWidth',2,'DisplayName', 'Marginal (smoothed) evidence for B.  $P(z_t=B \mid \mathbf{x}_{1:T})$');
            yl = ylim();
            stateSeq = obj.statesTrue{idx};
            state = stateSeq(1);
            lastIdx = 1;
            for i=1:size(Traj,1)
                if state ~= stateSeq(i)
                    if state
                        xs=[lastIdx,i+1, i+1, lastIdx];
                        ys=[yl(1),yl(1), yl(2), yl(2)];
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                    lastIdx=i;
                    state = stateSeq(i);
                end
            end
            ylim(yl);
            
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize=12;
            xlabel('Observation index');
            ylabel('Distance (um)');
            
            subplot(2,1,2);
            trueStateSeq = obj.statesTrue{idx};
            T = obj.transitions;
            llhTru = PairInteractionMarkovModel.logLikelihood_HMM(trueStateSeq,obj.prior,E,T);
            llhMPE = PairInteractionMarkovModel.logLikelihood_HMM(obj.statesMAP{idx},obj.prior,E,T);
            llhMPM = PairInteractionMarkovModel.logLikelihood_HMM(obj.statesMPM{idx},obj.prior,E,T);
            stairs(1:N, obj.statesMAP{idx}+0.04, 'Color',[1,0.0,0.0],'Linewidth',3,'DisplayName',...
                    sprintf('[LLH = %.6g] MAP (Max \\emph{a Posteriori}) (Viterbi). $\\mathrm{arg}\\max_{z_{1:T}} P\\left( z_{1:T} \\mid \\mathbf{x}_{1:t} \\right)$',llhMPE));
            hold();
            stairs(1:N, obj.statesMPM{idx}+0.08, 'Color',[0,1,0.5],'Linewidth',3,'DisplayName',...
                sprintf('[LLH = %.6g] MPM (Max posterior marginal).',llhMPM));
            stairs(1:N, trueStateSeq, 'Color','k','Linewidth',1,'DisplayName',...
                sprintf('[LLH = %.6g] True state',llhTru));
            
            
            set(gca(),'YTick',[0,1]);
            L=legend('location','best');
            L.Interpreter='latex';
            L.FontSize = 12;
            ylim([0,2]);
            xlabel('Observation index');
            ylabel('State (shifted to see)');            
        end
        
        function plotHMMObservation(obj, idx)
            T = obj.data{idx};
            N = size(T,1);
            d = sqrt((T(:,2)-T(:,4)).^2 + (T(:,3)-T(:,5)).^2); %observed distance
            E = obj.evidence{idx};
            figure('Position',[0,0, 600, 1000]);
            subplot(4,1,1);
            plot(1:N,d,'k-','DisplayName', 'Observed distance');
            hold();
            L=legend('location','best');
            L.Interpreter='latex';
            title('Observed Distance');
            xlabel('Observation index');
            ylabel('Distance (um)');
            
            subplot(4,1,2);
            plot(1:N,log(E(:,1)),'r-','DisplayName','Free evidence $P(z_t=F | x_t, x_{t-1})$');
            hold();
            plot(1:N,log(E(:,2)),'b-','DisplayName','Bound evidence $P(z_t=B | x_t)$');
            L=legend('location','best');
            L.Interpreter='latex';
            title('Log-Evidence --- Distance Model (Low-Nam, et. al.)');
            xlabel('Observation index');
            ylabel('Log-Probability');
            
            fwdFilter = obj.fwdFilt{idx};
            bkwFilter = obj.bkwFilt{idx};
            marg = obj.marginal{idx};
            subplot(4,1,3);
            pB1 = E(:,2);
            pBf = fwdFilter(:,2);
            pBb = bkwFilter(:,2);
            pBm = marg(:,2);
            plot(1:N, pB1,'k-','DisplayName', 'Immediate evidence for B.  $P(z_t=B | x_t)$');
            hold;
            plot(1:N, pBf,'r-','DisplayName', 'Forwards filtered evidence for B.  $P(z_t=B | x_{1:t})$');
            plot(1:N, pBb,'b-','DisplayName', 'Backwards filtered evidence for B. $P(x_{t+1:t}| z_t=B)$');
            plot(1:N, pBm,'g-','DisplayName', 'Marginal (smoothed) evidence for B.  $P(z_t=B | x_{1:T})$');
            L=legend('location','best');
            L.Interpreter='latex';
            title('Fwd-Backwd Marginal Evidence --- Distance Model (Low-Nam, et. al.)');
            xlabel('Observation index');
            ylabel('Probability of Bound');
            
            subplot(4,1,4);
            T = obj.transitions;
            llhMPE = PairInteractionMarkovModel.logLikelihood_HMM(obj.statesMAP{idx},obj.prior,E,T);
            llhMPM = PairInteractionMarkovModel.logLikelihood_HMM(obj.statesMPM{idx},obj.prior,E,T);
            stairs(1:N, obj.statesMAP{idx}+0.02, 'Color',[1,0.5,0.5],'Linewidth',2,'DisplayName',...
                    sprintf('[LLH=%.6g] MAP (Max a-posteriori) (viterbi). $\\max_{\\vec{z}} P(\\vec{z}|\\vec{x}_{1:t})$',llhMPE));
            hold();
            stairs(1:N, obj.statesMPM{idx}+0.04, 'Color',[0,1,0.5],'Linewidth',1,'DisplayName',...
                sprintf('[LLH=%.6g] MPM (Max posterior marginal).',llhMPM));
            if ~isempty(obj.statesTrue)
                trueStateSeq = obj.statesTrue{idx};
                llhTru = PairInteractionMarkovModel.logLikelihood_HMM(trueStateSeq,obj.prior,E,T);
                stairs(1:N, trueStateSeq, 'Color','k','Linewidth',1,'DisplayName',...
                    sprintf('[LLH=%.6g] True state',llhTru));
            end
            L=legend('location','best');
            L.Interpreter='latex';
            ylim([0,2]);
            title('Prediction --- Distance Model (Low-Nam, et. al.)');
            xlabel('Observation index');
            ylabel('State (shifted to see)');
            
        end
        
        function [S, S_llh, W] = particleSMC(obj, obs_idx, Nsamples)
            
        end
        
        function [S, S_llh, P, P_llh] = sampleMCMC(obj, obs_idx, Nsamples)
            %Initially just vary the states, eventually also vary input parameters
            Nburnin = ceil(Nsamples*0.2);
            Ntotal = Nsamples +  Nburnin;
            initS = obj.statesMAP{obs_idx};
            K = numel(initS);
            initS = false(K,1);
            E = obj.evidence{obs_idx};
            T = obj.transitions;
            initS_llh = obj.logLikelihood_HMM(initS, obj.prior, E, T);            
            
            S = zeros(Ntotal,K);
            S_llh = zeros(Ntotal,1);
            P = zeros(Ntotal, K);
            P_llh = zeros(Ntotal, 1);
            S(1,:) = initS;
            P(1,:) = initS;
            S_llh(1) = initS_llh;
            P_llh(1) = initS_llh;
            Naccepted = 0;
            Nrejected = 0;
            for n = 2:Ntotal
                pos = randi(K,1,1);
                P(n,:) = S(n-1,:);
                P(n,pos) = ~P(n,pos);
                P_llh(n) = obj.logLikelihood_HMM(P(n,:), obj.prior, E, T);
                if P_llh(n) >= S_llh(n-1) || log(rand()) < P_llh(n) - S_llh(n-1)
                    %accept
                    Naccepted = Naccepted + 1;
                    S(n,:) = P(n,:);
                    S_llh(n) = P_llh(n);
                else
                    %reject
                    Nrejected = Nrejected + 1;
                    S(n,:) = S(n-1,:);
                    S_llh(n) = S_llh(n-1);
                end
            end
            fprintf('Accepted %i/%i samples. [Rate: %.3g]\n',Naccepted, Ntotal, Naccepted/Ntotal);
            S = S(Nburnin:end,:);
            S_llh = S_llh(Nburnin:end);
            P = P(Nburnin:end,:);
            P_llh = P_llh(Nburnin:end);
        end
        
        
        function plotMCMCSamples(obj, S, S_llh, P, P_llh)
            f=figure('Position',[10,10,900,700]);
            N = size(S,1);
            plot(1:N, S_llh, 'o', 'MarkerSize',2,'Color', [0,0,1], 'DisplayName','Accepted');
            hold('on');
            plot(1:N, P_llh, 'o', 'MarkerSize',2, 'Color',[0.5,0,0], 'DisplayName','Proposed');
            xlabel('Sequence');
            ylabel('Log-likelihood');
            legend('Location','Best');
        end
        
        function [En1, Enjk, Enj] = calculateExpectedCounts(obj)
            En1 = zeros(1,2);
            Enjk = zeros(2,2);
            Enj = zeros(1,2);
            N=numel(obj.data);
            for k=1:N
                marg = obj.marginal{k};
                marg2 = obj.marginal2Slice{k};
                En1 = En1 + marg(1,:);
                Enj = Enj + sum(marg);
                Enjk = Enjk + sum(marg2,3);
            end
        end

        function [T,kon,koff] = estimateTransitions_knownState(obj)
            T = zeros(2,2);
            for n=1:numel(obj.data)
                states = obj.statesTrue{n}+1;
                for t=2:numel(states)
                    T(states(t-1),states(t)) = T(states(t-1),states(t)) + 1;
                end
            end
            T = T ./ repmat(sum(T,2),1,2);
            kon = log(1-T(1,2))/(-obj.deltaT);
            koff = log(1-T(2,1))/(-obj.deltaT);
        end
        
        function [newT, rates] = maximizeTransitionProbMPM(obj, model)
            Ps = obj.modelParams.(model);
            ND = numel(obj.data);
            function llh = ensembleLLH(ks)
                llh=0;
                T = PairInteractionMarkovModel.transitionRate_Markov(obj.deltaT, ks(1), ks(2));
                for k=1:ND
                    E = obj.evidence{k};
                    [~,~,marg] = PairInteractionMarkovModel.computeMarginalFB(obj.prior, E, T);
                    mpm = PairInteractionMarkovModel.computeMPMstates(marg);
                    llh = llh - PairInteractionMarkovModel.logLikelihood_HMM(mpm, obj.prior, E, T);
                end
            end
            function stop = outFun(x, vals, state)
                stop=false;
                fprintf('[%i] kon:%.9g koff:%.9g\n',vals.iteration,x(1),x(2));
            end
            problem.objective = @ensembleLLH;
            problem.x0 = [Ps.kon, Ps.koff];
            problem.lb = [0, 0];
            problem.solver = 'fmincon';
            problem.options = optimoptions('fmincon');
            problem.options.Algorithm = 'interior-point';
            problem.options.Hessian = 'bfgs';
            problem.options.OutputFcn= @outFun;
            [rates, maxLLH, flag, out, lambda, grad, hess] = fmincon(problem);
            newT = obj.transitionRate_Markov(obj.deltaT, rates(1), rates(2));
            obj.modelParams.(model).kon = rates(1);
            obj.modelParams.(model).koff = rates(2);
        end

        
        function [newT, newPrior] = estimateTransitions(obj)
            N=numel(obj.data);
            [En1, Enjk, Enj] = obj.calculateExpectedCounts();
            newT = Enjk ./ repmat(sum(Enjk,2), 1,2);
            newPrior = En1./N;
            ND = numel(obj.data);
            obj.transitions = newT;
            obj.prior = newPrior;
            for k=1:ND           
                [obj.fwdFilt{k}, obj.bkwFilt{k}, obj.marginal{k}] = ...
                    obj.computeMarginalFB(obj.prior, obj.evidence{k}, obj.transitions);
                obj.marginal2Slice{k} = obj.computeMarginal2Slice(obj.prior, obj.evidence{k}, obj.transitions);
                obj.statesMAP{k} = obj.computeViterbi(obj.prior, obj.evidence{k}, obj.transitions);
                obj.statesMPM{k} = obj.computeMPMstates(obj.marginal{k}); % maximum of posterior marginals.  The marginally most probable.
            end
        end
    end

    methods (Access=protected)
        function calculateModelHMM(obj)
            ND = numel(obj.data);
            obj.evidence = cell(ND,1);
            obj.transitions = cell(ND,1);
            obj.fwdFilt = cell(ND,1);
            obj.bkwFilt = cell(ND,1);
            obj.marginal = cell(ND,1);
            obj.marginal2Slice = cell(ND,1);
            obj.statesMAP = cell(ND,1);
            obj.statesMPM = cell(ND,1);
            P = obj.modelParams.(obj.modelName);
            obj.deltaT = min(diff(obj.data{1}(:,1)));
            obj.transitions = obj.transitionRate_Markov(obj.deltaT,P.kon, P.koff);
            for k=1:ND
                switch obj.modelName
                    case 'Distance'
                        obj.evidence{k} = obj.evidence_Distance(obj.data{k},P);
                    case 'Diffusive'
                        obj.evidence{k} = obj.evidence_Diffusive(obj.data{k},P);
                    otherwise
                        error('PairInteractionMarkovModel:initializeModelParams','Unknown model %s', obj.modelName);
                end
                [obj.fwdFilt{k}, obj.bkwFilt{k}, obj.marginal{k}] = ...
                    obj.computeMarginalFB(obj.prior, obj.evidence{k}, obj.transitions);
                obj.marginal2Slice{k} = obj.computeMarginal2Slice(obj.prior, obj.evidence{k}, obj.transitions);
                obj.statesMAP{k} = obj.computeViterbi(obj.prior, obj.evidence{k}, obj.transitions);
                obj.statesMPM{k} = obj.computeMPMstates(obj.marginal{k}); % maximum of posterior marginals.  The marginally most probable.
            end
        end
    end
    
    methods (Static=true, Access=public)

        function params = makeModelParams_Distance(params_vec)
            params.D = params_vec(1);
            params.kon = params_vec(2);
            params.koff = params_vec(3);
            params.rho = PairInteractionMarkovModel.DEFAULT_RHO;
            params.sigmaReg = PairInteractionMarkovModel.DEFAULT_SIGMA_REG;            
        end

        function params = makeModelParams_Diffusive(params_vec)
            params.DA = params_vec(1);
            params.DB = params_vec(2);
            params.DAB = params_vec(3);
            params.kon = params_vec(4);
            params.koff = params_vec(5);
            params.rho = PairInteractionMarkovModel.DEFAULT_RHO;
            params.sigmaReg = PairInteractionMarkovModel.DEFAULT_SIGMA_REG;            
        end
        
        function T = transitionRate_Markov(deltaT,kon, koff)
            T=zeros(2,2);
            T(1,2,:) = 1-exp(-deltaT*kon);
            T(2,1,:) = 1-exp(-deltaT*koff);
            T(1,1,:) = 1-T(1,2);
            T(2,2,:) = 1-T(2,1);
        end
        
        function [E,varB,varF] = evidence_Distance(data, params)
            % The E(i,j) = log-evidence for state j at time i.
            % free=state0, bound=state1
            N = size(data,1);
            E = zeros(N,2);
            d2 = (data(:,2)-data(:,4)).^2 + (data(:,3)-data(:,5)).^2; %observed distance squared
            d = sqrt(d2);
            dT = diff(data(:,1));
            varA = max(data(:,6:7),[],2).^2; % Effective variance in A as max of X&Y varainces
            varB = max(data(:,8:9),[],2).^2; % Effective variance in B as max of X&Y varainces
            
            varB =  0.5*(varA+varB) + params.sigmaReg^2 + params.rho^2;
            varF =  0.5*(varA(2:end) + varB(2:end) + varA(1:end-1) + varB(1:end-1))...
                    + params.sigmaReg^2 + 4*dT*sqrt(params.D1^2+params.D2^2);
            E(:,2) = log(d) - log(2*pi*varB) + -d2./(2*varB); %Bound evidence
            nsteps=360;
            theta = linspace(0,2*pi,nsteps);
            %Free evidence
            E(2:end,1) = log(d(2:end)) - log(2*pi*varF) + ...
                log(1/nsteps*sum(exp(repmat(-1./(2*varF),1,nsteps) .* ...
                (repmat(d2(2:end)+d2(1:end-1),1,nsteps)-2*d(2:end).*d(1:end-1)*sin(theta))),2));
            E(1,1) = log(2*d(1))-2*log(max(d)); 
            for k=1:N
                E(k,:) =  exp(E(k,:) - PairInteractionMarkovModel.logSum(E(k,:)));
            end
        end

        function E = evidence_Diffusive(data, params)
            % The E(i,j) = log-evidence for state j at time i.
            % free=state0, bound=state1
            N = size(data,1);
            E = zeros(N,2);
            exposureT = min(diff(data(:,1)));
            icp = InteractionChangePoint(params.DA, params.DB, params.DAB, exposureT, params.rho);
            P = zeros(N, 18);
            P(:,[1,2,3,6,7,10,11,14,15]) = data(:,1:9);
            icp.initializePair(P);
            d = sqrt(sum( (P(:,2:3)-P(:,6:7)).^2,2));
            varA = max(data(:,6:7),[],2).^2; % Effective variance in A as max of X&Y varainces
            varB = max(data(:,8:9),[],2).^2;
            varB =  0.5*(varA+varB) + params.sigmaReg^2 + params.rho^2;
            E(1,2) = log(d(1)) - log(2*pi*varB(1)) + -d(1)^2./(2*varB(1));
            E(1,1) = log1p(-exp(E(1,2)));
            
            E(2:end,1) = diag(icp.freeLLH,-1);
            E(2:end,2) = diag(icp.boundLLH,-1);                   
            for k=1:N
                E(k,:) =  exp(E(k,:) - PairInteractionMarkovModel.logSum(E(k,:)));
            end
        end

        function [alpha, beta, gamma, Z] = computeMarginalFB(prior, E, T)
            % Use the forwards algorithm to return the filtered belief state
            % Evidence and alpha/beta/gamma vectors are in log-space.
            N = size(E,1);
            K = size(E,2);
            alpha = zeros(N,K);
            beta = zeros(N,K);
            Z = zeros(N,1);
            alpha(1,:) = E(1,:) .* prior(:)';
            Z(1) = sum(alpha(1,:));
            alpha(1,:) = alpha(1,:) ./ Z(1);
            beta(N,:) = [0.5,0.5];
            for n=2:N
                alpha(n,:) = E(n,:) .* alpha(n-1,:) * T;
                beta(N-n+1,:) = E(N-n+2,:) .* beta(N-n+2,:) * T;
                Z(n) = sum(alpha(n,:));
                alpha(n,:) = alpha(n,:) ./ Z(n); 
                beta(N-n+1,:) = beta(N-n+1,:) ./ sum(beta(N-n+1,:)); 
            end
            gamma = alpha .* beta;
            gamma = gamma ./ repmat(sum(gamma,2),1,2);
        end
        
        function eta = computeMarginal2Slice(prior, E, T)
            N = size(E,1);
            K = size(E,2);
            [alpha, beta] = PairInteractionMarkovModel.computeMarginalFB(prior, E, T);
            eta = zeros(K,K,N-1);
            for n=1:N-1
                eta(:,:,n) = T .* repmat(alpha(n,:)',1,2) .* repmat(beta(n+1,:) .* E(n+1,:),2,1);
                eta(:,:,n) = eta(:,:,n) ./ repmat(sum(eta(:,:,n),2),1,2);
            end
        end
        
        function [Zs, xi, alpha] = samplePosterior(nSamples, E, T, fwdFilter)
            % [out]
            %  Zs = vector of samples row=t, col=sample#, 0=free, 1=bound
            %  xi = 
            K = size(E,2);
            N = size(E,1);
            alpha = fwdFilter;
            xi = zeros(K,K,N-1); % p(z_t | z_{t+1}, x_{1:t}).  row = z_t, col = z_t+1, slice = t (1..N-1)
            Zs = false(N,nSamples); % sampled state sequence 0=free 1=bound, row=t col=sample#
            Zs(N,:) = rand(1,nSamples) >= alpha(N,1);
            for n=N-1:-1:1
                xi(:,:,n) = repmat(E(n+1,:).*alpha(n+1,:),2,1) .* T .* repmat(alpha(n,:)',1,2);
                Q = sum(xi(:,:,n),2);
                xi(:,:,n) = xi(:,:,n) - repmat(Q,2,1);
                W = squeeze(xi(1,:,n)); %P(z_t=F | z_{t+1}, x_{1:t}).
                Zs(n,:) = rand(1,nSamples) >= W(Zs(n+1,:)+1);
            end
        end
        
        function [mpe, delta, A] = computeViterbi(prior, E, T)
            K = size(E,2);
            N = size(E,1);
            A = zeros(N,K);
            mpe = zeros(N,1);
            delta = zeros(N,K);
            delta(1,:) = prior .* E(1,:);
            delta(1,:) = delta(1,:) ./ sum(delta(1,:));
            for n=2:N
                [delta(n,:),A(n,:)] = max(repmat(delta(n-1,:)',1,2) .* T .* repmat(E(n,:),2,1));                
                delta(n,:) = delta(n,:) ./ sum(delta(n,:));
            end
            [~,mpe(N)] = max(delta(n,:));
            for n=N-1:-1:1
                mpe(n) = A(n+1,mpe(n+1));
            end
            mpe = mpe-1; % 0=free, 1=bound
        end

        function llh = logLikelihood_HMM(stateSeq, prior, E, T)
            N = numel(stateSeq);
            stateSeq = stateSeq+1;
            llh = log(prior(stateSeq(1))) + log(E(1,stateSeq(1)));
            for n=2:N
                llh = llh + log(E(n,stateSeq(n))) + log(T(stateSeq(n-1), stateSeq(n)));
            end
        end
        
       function logS = logSum(v)
            a=max(v,[],2);
            b=min(v,[],2);
            logS = a + log1p(exp(b-a));
        end

        function mpm = computeMPMstates(marginal)
            mpm = marginal(:,1) < 0.5;
        end
      
    end
    
end

