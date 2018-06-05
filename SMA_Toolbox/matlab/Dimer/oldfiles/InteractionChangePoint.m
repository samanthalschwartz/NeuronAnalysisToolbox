%
% InteractionChangePoint - Detect inetaractions between pairs of particles, using Bayesian changepoint analysis.
% This relies on the data format established in the PairAnalysis class for storing single particle trajectory
% pairs.
%
% Mark J. Olah (mjo@cs.unm.edu)
% 07/2015
%

classdef InteractionChangePoint < IfaceMixin
    
    properties
        D_A;
        D_B;
        D_AB;        
        exposureT;
        rho; %This is the gaussian sigma (SE) for the dimer distance
        
        beta; % The penalty for adding a changepoint

        nObs; %The number of simultaneous observations of both molecules in the pair.  
              %This is the number of rows in obj.pair
       

        pair; % Matrix of the pair tajectory information in the format defined by PairAnalysis
        dist; % Matrix of the pair distance calculations as defined by PairAnalysis
       
        
             
        freeLLH; % matrix size:[nObs,nObs] element [i,j] where i<j gives the likelihood of 
                       % observations i:j occurring when the particles are free
        boundLLH; % matrix size:[nObs,nObs] element [i,j] where i<j gives the likelihood of 
                       % observations i:j occurring when the particles are bound
    end
    
    properties (SetAccess=protected)
        initialized=false;
    end
    
    methods
        function obj = InteractionChangePoint(varargin)
            % Initialization of the InteractionChangePoint for a trajectory pair.   This
            % works on a single pair at a time.
            %
            % [in] pair - A trajectory pair given in the PairAnalysis matrix format.
            obj = obj@IfaceMixin(@InteractionChangePoint_Iface);
            if nargin>0
                obj.initializeParameters(varargin{:})
            end
            
        end
        
        function initializeParameters(obj, D_A, D_B, D_AB, exposureT, rho)
            %Same arguments as constructor.  This allows the class to be reinitialized with a new pair
            obj.D_A =  double(D_A);
            obj.D_B = double(D_B);
            obj.D_AB = double(D_AB);
            obj.exposureT = double(exposureT);
            obj.rho = double(rho);
            if obj.initialized
                obj.closeIface()
            end
            obj.openIface(obj.D_A, obj.D_B, obj.D_AB, obj.exposureT, obj.rho); %Create C++ object & open interface
            %Reset pair computational data
            obj.pair = [];
            obj.nObs = [];
            obj.dist = [];
            obj.freeLLH = [];
            obj.boundLLH = [];
            obj.initialized = true;            
        end
        
        function initializePair(obj, pair)
            %Given a pair matix compute relvent distances and liklihoods for later analysis
            PairAnalysis.checkPairFormat(pair);
            obj.pair = pair;
            obj.nObs = size(pair,1);
            obj.beta = -log(obj.nObs);
            obj.dist = PairAnalysis.makePairDists(obj.pair);
            obj.freeLLH = obj.computeFreePairLikelihood(obj.pair);
            obj.boundLLH = obj.computeBoundPairLikelihood(obj.pair);
        end
    
        function stateSeq = interactionEventTimesToStateSequence(obj, bindingEvents, unbindingEvents)
            % Given a vector of binding event times and unbinding event times, generate a discrete sequence
            % of states where stateSeq(i)=0 if free at the end of time sequence i
            %                 stateSeq(i)=1 if bound at the end of state sequence i.
            % 
            stateSeq = zeros(obj.nObs,1);
            nBindings = numel(bindingEvents);
            events = [bindingEvents(:); unbindingEvents(:)];
            [sevents, sevents_idx] = sort(events);
            tmin = obj.pair(1,1);
            tmax = obj.pair(end,1)+obj.exposureT;
            next_event_begin_idx = 1;
            for i=1:numel(events)
                if sevents(i)<tmin || sevents(i)>tmax
                    error('InteractionChangePoint:interactionEventsToStateSequence','Event time %f is not valid',events(i));
                end
                event_time_idx = find(sevents(i)<obj.pair(:,1), 1,'first')-1;
                stateSeq(next_event_begin_idx:event_time_idx-1) = sevents_idx(i)>nBindings;                
                next_event_begin_idx = event_time_idx;
            end
        end        

        function [state0, cps, LLH] = segmentProductMLE(obj, beta)
            [state0, cps, LLH]=obj.call('segmentProductMLE',obj.freeLLH, obj.boundLLH, beta);
            cps = cps+1;
        end
        
        function [LLH0, LLH1, cps] = compute1CPLLH(obj, seq_beg, seq_end)
            if nargin<3
                seq_end = obj.nObs;
            end
            if nargin<2
                seq_beg = 1;
            end
            cps = seq_beg+1:seq_end-1;
            seq_beg = int64(seq_beg)-1; % 0-based indexing;
            seq_end = int64(seq_end)-1; % 0-based indexing
            LLH0 = obj.call('compute1CPLLH',obj.pair, seq_beg, seq_end, false);
            LLH1 = obj.call('compute1CPLLH',obj.pair, seq_beg, seq_end, true);
        end

        function [LLH0, LLH1] = compute2CPLLH(obj, seq_beg, seq_end)
            if nargin<3
                seq_end = obj.nObs;
            end
            if nargin<2
                seq_beg = 1;
            end
            seq_beg = int64(seq_beg)-1; % 0-based indexing;
            seq_end = int64(seq_end)-1; % 0-based indexing
            LLH0 = obj.call('compute2CPLLH',obj.pair, seq_beg, seq_end, false);
            LLH1 = obj.call('compute2CPLLH',obj.pair, seq_beg, seq_end, true);
        end

        function f = plot1CPLLH(obj)
            [LLH0, LLH1, cps] = obj.compute1CPLLH();
            f=figure();
            subplot(2,1,1);
            plot(cps,LLH0,'r-','DisplayName','F->B')
            hold on;
            plot(cps,LLH1,'b-','DisplayName','B->F')
            plot([1,obj.nObs],repmat(obj.freeLLH(end,1),1,2),'g-','DisplayName','Free');
            plot([1,obj.nObs],repmat(obj.boundLLH(end,1),1,2),'m-','DisplayName','Bound');
            yl=ylim();
            [mLLH0, mLLH0x] = max(LLH0);
            [mLLH1, mLLH1x] = max(LLH1);
            plot([mLLH0x+1, mLLH0x+1], yl, 'r:','DisplayName', 'MLE F->B');
            plot([mLLH1x+1, mLLH1x+1], yl, 'b:','DisplayName', 'MLE B->F');
            hold off;
            xlabel('Change Point');
            ylabel('LLH');
            legend('location','best');

            subplot(2,1,2);            
            D = PairAnalysis.makePairDists(obj.pair);
            xs=1:obj.nObs-1;
            hold on;
            plot(xs,D(xs,2),'r-');
            state0 = mLLH1>mLLH0
            if state0
                cps = mLLH1x+1;
            else
                cps = mLLH0x+1;
            end
            stateSeq = obj.cpsToStateSeq(cps, state0, 1, obj.nObs);
            yl = ylim();
            stairs(1:obj.nObs,stateSeq*yl(2)*0.95)
            hold off;
        end

        function f = plot2CPLLH(obj)
            [LLH0, LLH1] = obj.compute2CPLLH();
            f=figure();
            subplot(2,1,1);
            hold on;
            [mLLH0, mLLH0x] = max(LLH0(:));
            [mLLH1, mLLH1x] = max(LLH1(:));
            if mLLH0>mLLH1
                state0 = false;
                [i,j] = ind2sub(size(LLH0),mLLH0x);
                cps = [1+j, i];
            else
                state0 = true;
                [i,j] = ind2sub(size(LLH1),mLLH1x);
                cps = [1+j, i];
            end
            plot([1,obj.nObs],repmat(obj.freeLLH(end,1),1,2),'g-','DisplayName','Free');
            plot([1,obj.nObs],repmat(obj.boundLLH(end,1),1,2),'m-','DisplayName','Bound');
            yl=ylim();
            for cp1=1:obj.nObs-2
                possiblecps=(cp1+1:obj.nObs-2);
                if ~state0 && cp1==cps(1)
                    plot(possiblecps,LLH0(possiblecps,cp1),'r-','LineWidth',3, 'DisplayName','F->B');
                    plot([cps(2), cps(2)], yl, 'r:', 'DisplayName', 'MLE F->B');
                else
%                     plot(possiblecps,LLH0(possiblecps,cp1),'r-','LineWidth',0.5,'DisplayName','F->B');
                end
                if state0 && cp1==cps(1)
                    plot(possiblecps,LLH1(possiblecps,cp1),'b-','LineWidth',3, 'DisplayName','B->F');
                    plot([cps(2), cps(2)], yl, 'b:', 'DisplayName', 'MLE B->F');
                else
%                     plot(possiblecps,LLH1(possiblecps,cp1),'b-','DisplayName','B->F');
                end
            end
            hold off;
            xlabel('Change Point');
            ylabel('LLH');
%             legend('location','best');

            subplot(2,1,2);            
            D = PairAnalysis.makePairDists(obj.pair);
            xs=1:obj.nObs-1;
            hold on;
            plot(xs,D(xs,2),'r-');
            stateSeq = obj.cpsToStateSeq(cps, state0, 1, obj.nObs);
            yl = ylim();
            stairs(1:obj.nObs,stateSeq*yl(2)*0.95)
            hold off;
        end

        
        function f=plotSegmentProductMLE(obj, beta)
            if nargin==1
                beta=obj.beta;
            end
            [state0,cps,LLH] = obj.segmentProductMLE(beta);
%             LLHfull=obj.computeSequenceLikelihood(state0,cps);
            xs = (1:obj.nObs)';
            f=figure();
            plot(xs(1:end-1), diag(obj.freeLLH,-1),'b-','DisplayName','FreeModel-1');
            hold on;
            plot(xs(1:end-2), diag(obj.freeLLH,-2),'b-','DisplayName','FreeModel-2');
            plot(xs(1:end-3), diag(obj.freeLLH,-3),'b-','DisplayName','FreeModel-3');
            plot(xs(1:end-1), diag(obj.boundLLH,-1),'r-','DisplayName','BoundModel-1');
            plot(xs(1:end-2), diag(obj.boundLLH,-2),'r-','DisplayName','BoundModel-2');
            plot(xs(1:end-3), diag(obj.boundLLH,-3),'r-','DisplayName','BoundModel-3');
            stateSeq = obj.cpsToStateSeq( cps, state0);
%             D = PairAnalysis.makePairDists(obj.pair);

%             plot(xs,D(xs,2),'r-');
            yl = ylim();
            hold on;
            stairs(xs,(stateSeq*(yl(2)-yl(1))*0.95+yl(1)*0.975), 'g-');
            hold off;
                        
            xlabel('Jump number');
            ylabel('Model LLH');
            legend('Location','best');
        end
        
%         function LLH = sequenceLikelihood(obj, stateSeq)
%             % stateSeq : vector length:nObs giving state at end of each frame aquisition.  
%             %               stateSeq(i)=0 if free
%             %               stateSeq(i)=1 if bound
%             state = stateSeq(1);
%             cps = find(diff(stateSeq));
%             cps = [1; cps; obj.nObs]; % all change points
%             LLH = 0;
%             for i=1:numel(cps)-1
%                 if state %bound
%                     LLH = LLH + obj.boundLLH(cps(i+1),cps(i)); 
%                 else %free
%                     LLH = LLH + obj.freeLLH(cps(i+1),cps(i)); 
%                 end
%                 state = ~state;
%             end
%         end
%         
%         function [cp, LLH] = maxLikelihood1CPseq(obj, init_state, seq_beg, seq_end)
%             % seq_beg<cp<=seq_end
%             if nargin<4
%                 seq_end = obj.nObs;
%             end
%             if nargin<3
%                 seq_beg = 1;
%             end
%             cps = seq_beg+1:seq_end-1;
%             if init_state %initial state == bound
%                 LLH_seq = obj.boundLLH(cps, seq_beg) + flip(obj.freeLLH(seq_end, cps))'; %state 1 -> state 0
%             else
%                 LLH_seq = obj.freeLLH(cps, seq_beg) + flip(obj.boundLLH(seq_end, cps))'; %state 0 -> state 1
%             end
%             f=figure();
%             plot(cps, LLH_seq,'r-');
%             [LLH, cp] = max(LLH_seq);
%             cp = cp + seq_beg;
%         end
%         
%         function [cps, LLH] = maxLikelihood2CPseq(obj, init_state, seq_beg, seq_end)
%             % 1<cp<=nObs
%             if nargin<4
%                 seq_end = obj.nObs;
%             end
%             if nargin<3
%                 seq_beg = 1;
%             end
%             if init_state
%                 A=obj.boundLLH;
%                 B=obj.freeLLH;
%             else
%                 A=obj.freeLLH;
%                 B=obj.boundLLH;
%             end
%             N = seq_end-seq_beg-2; %number of options for first change point
%             LLHs=zeros(N,1);
%             cps=zeros(N,2);
%             cp1s = seq_beg+1:seq_end-2;
%             cps(:,1)=cp1s;
%             for n = 1:N
%                 cp2s = cp1s(n)+1:seq_end-1;
%                 LLH_seq = A(cp1s(n),seq_beg) + B(cp2s,cp1s(n)) + flip(A(seq_end,cp2s))';                
%                 [LLHs(n), cp2idx] = max(LLH_seq);
%                 cps(n,2) = cp2s(cp2idx);
%             end
%             f=figure();
%             subplot(1,2,1);
%             plot(cp1s, LLHs,'r-');
%             xlabel('cp1');
%             ylabel('LLH');
%             subplot(1,2,2);
%             plot(cps(:,1),cps(:,2),'b-');
%             xlabel('cp1');
%             ylabel('cp2');
%             
%             [LLH, cp1idx] = max(LLHs);
%             cps = cps(cp1idx,:);
%         end
        
        function stateSeq = cpsToStateSeq(obj, cps, init_state, seq_beg, seq_end)
            if nargin<4
                seq_beg = 1;
            end
            if nargin<5
                seq_end = obj.nObs;
            end
            N = seq_end-seq_beg+1;
            cps = [seq_beg;cps(:);seq_end];
            stateSeq = false(N,1);
            for n=1:numel(cps)-1
                stateSeq(cps(n):cps(n+1)) = xor(init_state, mod(n+1,2));
            end
        end
        
        
        function f = plot1CPMLE(obj, seq_state0, seq_beg, seq_end)
            if nargin<4
                seq_end = obj.nObs;
            end
            if nargin<3
                seq_beg = 1;
            end
            cps = seq_beg+1:seq_end-1;
            f=figure();
            subplot(1,2,1);
            hold on;
            subplot(1,2,2);
            hold on;
            for cp = cps
                [LLH, comps] = obj.computeSequenceLikelihoodDebug(seq_state0, cp, seq_beg, seq_end);
                comp_LLH_x = -0.5*(log(2*pi)+log(comps(5,:,1))+(comps(3,:,1)-comps(4,:,1)).^2./comps(5,:,1));
                comp_LLH_y = -0.5*(log(2*pi)+log(comps(5,:,2))+(comps(3,:,2)-comps(4,:,2)).^2./comps(5,:,2));
                xs=1:size(comps,2);
                n=sprintf('cp=%i (x)',cp);
                c=[1,0,0]*((numel(cps)-cp+1)/numel(cps));
                subplot(1,2,1);
                plot(xs, comp_LLH_x,'-','Color',c,'DisplayName',n);
                n=sprintf('cp=%i (y)',cp);
                c=[0,0,1]*((numel(cps)-cp+1)/numel(cps));
                subplot(1,2,2);
                plot(xs, comp_LLH_y,'-','Color',c,'DisplayName',n);
                sLLH=sum([comp_LLH_x,comp_LLH_y]);
            end
            subplot(1,2,1);
            hold off;
            xlabel('Distribution');
            ylabel('Likelihood Contribution');
%             legend('location','best');
            subplot(1,2,2);
            hold off;
            xlabel('Distribution');
            ylabel('Likelihood Contribution');
%             legend('location','best');
        end
        
        function f = plot2CPMLE(obj, init_state, seq_beg, seq_end)
            if nargin<4
                seq_end = obj.nObs;
            end
            if nargin<3
                seq_beg = 1;
            end
            [cps, ~] = obj.maxLikelihood2CPseq(init_state, seq_beg, seq_end);
            stateSeq = obj.cpsToStateSeq( cps, init_state, seq_beg, seq_end);
            D = PairAnalysis.makePairDists(obj.pair);
            f=figure();
            xs=seq_beg:seq_end;
            plot(xs,D(xs,2),'r-');
            hold on;
            yl = ylim();
            stairs(xs,stateSeq*yl(2)*0.95)
            hold off;
        end
        
        function plotTransitionProdLLH(obj)
            f=figure();
            ax=axes();
            hold on;
            for e = seq_ends
                xs=2:e;
                plot(xs,obj.freeLLH(xs,1)+obj.freeLLH(e,xs)','r-');
                llh=obj.freeLLH(e,1);
                plot([1,e],[llh, llh],'k:');
            end
            hold off;
            xlabel('Free->Free Changepoint');
            ylabel('LLH');
        end

%         function [stateSeq, LLH] = maximumLikelihoodSequence(obj, nCPs)
%             
%         end
        
        function [cpPos, LLH] = maximumLikelihoodChangepointPosition(obj, seqBegin, seqEnd, initialState)
            % find the single changepoint that maximizes the LLH of the sequence from seqBegin to seqEnd
            % when changing from initalState to ~initialState at pos cpPos.  state(cpPos)==~initialState
            N = seqEnd-seqBegin+1;
            assert(N>1);
            if initialState
                initMat = obj.boundLLH;
                finalMat = obj.freeLLH;
            else
                initMat = obj.freeLLH;
                finalMat = obj.boundLLH;
            end
            LLHs = initMat(seqBegin+1:seqEnd ,seqBegin) + finalMat(seqEnd, seqBegin+1:seqEnd)';
            [LLH, cpPos] = max(LLHs);            
        end

        function [DmleA,DmleB] = estimateDMLE(obj)
            assert(~isempty(obj.freeLLH)); %check precompute is done
            DmleA = zeros(obj.nObs, obj.nObs);
            DmleB = zeros(obj.nObs, obj.nObs);
            for i=1:obj.nObs
                for j=i+1:obj.nObs
                    dest_a = DEstimator(obj.pair(i:j,[2,3]), obj.pair(i:j,1), obj.pair(i:j,[10,11]).^2, obj.exposureT);
                    dest_b = DEstimator(obj.pair(i:j,[6,7]), obj.pair(i:j,1), obj.pair(i:j,[14,15]).^2, obj.exposureT);
                    DmleA(j,i) = dest_a.MLE();
                    DmleB(j,i) = dest_b.MLE();
                end
            end
        end

        function LLH = computeFreePairLikelihood(obj, pair)
            PairAnalysis.checkPairFormat(pair);
            LLH = obj.call('freeSegmentLLH', pair);
        end
        
        function LLH = computeBoundPairLikelihood(obj, pair)
            PairAnalysis.checkPairFormat(pair);
            LLH = obj.call('boundSegmentLLH', pair);
        end

        function LLH = computeSequenceLikelihood(obj, beta, seq_state0, seq_cps, seq_begin, seq_end)
            if nargin<5
                seq_begin=1;
            end
            if nargin<6
                seq_end=obj.nObs;
            end
            seq_state0 = logical(seq_state0);
            seq_cps = int32(seq_cps)-1; %convert to 0-based indexing
            seq_begin = int32(seq_begin)-1; % 0-based indexing;
            seq_end = int32(seq_end)-1; %0-bsed indexing
            LLH = obj.call('sequenceLLH', obj.pair, seq_begin, seq_end, seq_state0, seq_cps);
            LLH = LLH + beta*length(seq_cps);
        end

        function [LLH,LLH_components] = computeSequenceLikelihoodDebug(obj,  seq_state0, seq_cps, seq_begin, seq_end)
            if nargin<4
                seq_begin=1;
            end
            if nargin<5
                seq_end=obj.nObs;
            end
            seq_state0 = logical(seq_state0);
            seq_cps = int32(seq_cps)-1; %convert to 0-based indexing
            seq_begin = int32(seq_begin)-1; % 0-based indexing;
            seq_end = int32(seq_end)-1; %0-bsed indexing
            [LLH,LLH_components] = obj.call('sequenceLLH_debug', obj.pair, seq_begin, seq_end, seq_state0, seq_cps);
            
        end
        
        
        
        
        function [mle_state0, mle_cps] = exhaustiveMLESequence(obj, seq_begin, seq_end)
             if nargin<2
                seq_begin=1;
            end
            if nargin<3
                seq_end=obj.nObs;
            end
            possible_cps=seq_begin+1:seq_end-1;
            maxLLH=-inf;
            mle_cps=[];
            mle_state0=[];
            max_cps=2;
            for ncps = 1:min(max_cps,numel(possible_cps))
                all_cps = nchoosek(possible_cps,ncps);
                fprintf('Considering %i sequences of %i cps\n',size(all_cps,1),ncps);
                for n = 1:size(all_cps,1)
                    LLH0 = obj.computeSequenceLikelihood(0, all_cps(n,:), seq_begin, seq_end);
                    if LLH0>maxLLH
                        maxLLH=LLH0;
                        mle_cps = all_cps(n,:);
                        mle_state0 = 0;
                    end
                    LLH1 = obj.computeSequenceLikelihood(1, all_cps(n,:), seq_begin, seq_end);
                    if LLH1>maxLLH
                        maxLLH=LLH1;
                        mle_cps = all_cps(n,:);
                        mle_state0 = 1;
                    end                    
                end
            end
        end
        
        
        function f=plotPairMatrixDistances(obj)
            f=figure();
            certainty=0.95;
            D = PairAnalysis.makePairDists(obj.pair, certainty);
            plot(obj.pair(:,end),D(:,2),'-r','DisplayName','Dist');
            hold on;
            xs = [obj.pair(:,end)', fliplr(obj.pair(:,end)')];
            ys = [D(:,3)', fliplr(D(:,4)')];
            fill(xs,ys,'r','EdgeColor','none','FaceAlpha',0.4,'DisplayName',sprintf('Dist (%.1f%% conf.)',certainty*100));
            hold off;
            xlabel('Time (frames)');
            ylabel('Dist (um)');
            legend('Location','best');            
        end
        
        function f=plotLengthNSegmentComparison(obj, segLength)
            %Plot relative likelihood for each segments of segLength being from
            % free or bound model.
            %  segLength : 1<= segLength <= nObs-2.  [default =1]
            if nargin==1
                segLength=1;
            end
            assert(all(1<=segLength) && all(segLength<=obj.nObs-2));
            f=figure();
            ax=axes();
            ax.Box='on';
            hold on;
            for n=1:numel(segLength)
                gamma = 0.15+0.85*n/numel(segLength);
                xs = (1:obj.nObs-segLength(n))';
                plot(xs, diag(obj.freeLLH,-segLength(n)), '-','Color',[0,0,gamma],'DisplayName',sprintf('Free Len:%i',segLength(n)));
                plot(xs, diag(obj.boundLLH,-segLength(n)),'-','Color',[gamma,0,0],'DisplayName',sprintf('Bound Len:%i',segLength(n)));
            end
            hold off;
            xlabel('Jump number');
            ylabel('Model LLH');
            legend('Location','best');
        end

        
%         function [changePoints, accLogBayesFactors, rejLogBayesFactors] = estimateGreedyCPs(obj, seq_begin, seq_end, logBayesFactorThreshold)           
%             logBF = obj.greedyLogBayesFactor(data);
%             if logBF>=logBayesFactorThreshold
%                 cp = obj.estimateChangePointLocation(data);
%                 data1=data(1:cp-1);
%                 data2=data(cp:end);
%                 [cP1, accLogBF1, rejLogBF1]=IntensityCPA.estimateChangePoints(data1,logBayesFactorThreshold);
%                 [cP2, accLogBF2, rejLogBF2]=IntensityCPA.estimateChangePoints(data2,logBayesFactorThreshold);
%                 accLogBayesFactors=[accLogBF1 logBF accLogBF2];
%                 rejLogBayesFactors=[rejLogBF1 rejLogBF2];
%                 changePoints=[cP1 cp cP2+cp-1];
%             else
%                 changePoints=[];
%                 rejLogBayesFactors=logBF;
%                 accLogBayesFactors=[];
%             end
%         end

        function f = plotGreedyChangePointsEstimate(obj,logBayesFactorThreshold, seq_beg, seq_end)
            if nargin<4
                seq_end = obj.nObs;
            end
            if nargin<3
                seq_beg = 1;
            end
            [LLH, state0, changePoints, accLogBayesFactors, rejLogBayesFactors] = obj.greedyChangePointsEstimate(logBayesFactorThreshold);
            stateSeq = obj.cpsToStateSeq( changePoints, state0, seq_beg, seq_end);
            D = PairAnalysis.makePairDists(obj.pair);
            f=figure();
            xs=seq_beg:seq_end;
            plot(xs,D(xs,2),'r-');
            hold on;
            yl = ylim();
            stairs(xs,stateSeq*yl(2)*0.95)
            hold off;
        end
            
        function [LLH, state0, changePoints, accLogBayesFactors, rejLogBayesFactors] = greedyChangePointsEstimate(obj,logBayesFactorThreshold)
            [logBF, ~, mle_state0, ~] = obj.greedyLogBayesFactor();
            if logBF>=logBayesFactorThreshold
                %At lest 1cp.  This starts the recursion                
                state0 = mle_state0;
                [changePoints, accLogBayesFactors, rejLogBayesFactors] = obj.greedyChangePointsEstimateRecurse(logBayesFactorThreshold, state0, 1, obj.nObs);
            else
                %No changepoints predicted choose the better of free or bound
                changePoints = [];
                accLogBayesFactors = [];
                rejLogBayesFactors = logBF;
                state0 = obj.boundLLH(end,1) > obj.freeLLH(end,1);
            end
            LLH = obj.computeSequenceLikelihood(state0, changePoints);
        end
        
        function [changePoints, accLogBayesFactors, rejLogBayesFactors] = greedyChangePointsEstimateRecurse(obj, logBayesFactorThreshold, seq_state0, seq_beg, seq_end)
            [logBF, mle_cps, ~] = greedyLogBayesFactorFixedState(obj, seq_state0, seq_beg, seq_end);
            if logBF>=logBayesFactorThreshold
                %At lest 1cp. Continue recursion
                if numel(mle_cps)==1
                    %one cp predicted
                    [cP1, accLBF1, rejLBF1] = obj.greedyChangePointsEstimateRecurse(logBayesFactorThreshold, seq_state0, seq_beg, mle_cps);
                    %now if cP1>0 we should maybe re-test our CP???
                    next_state =  xor(seq_state0, mod(numel(cP1),2)==0);
                    [cP2, accLBF2, rejLBF2] = obj.greedyChangePointsEstimateRecurse(logBayesFactorThreshold, next_state, mle_cps, seq_end);
                    accLogBayesFactors = [accLBF1, logBF, accLBF2];
                    rejLogBayesFactors = [rejLBF1, rejLBF2];
                    changePoints = [cP1, mle_cps, cP2];
                else
                    %two cps predicted
                    [cP1, accLBF1, rejLBF1] = obj.greedyChangePointsEstimateRecurse(logBayesFactorThreshold, seq_state0, seq_beg, mle_cps(1));
                    %now if cP1>0 we should maybe re-test our CP???
                    next_state =  xor(seq_state0, mod(numel(cP1),2)==0);
                    [cP2, accLBF2, rejLBF2] = obj.greedyChangePointsEstimateRecurse(logBayesFactorThreshold, next_state, mle_cps(1), mle_cps(2));
                    
                    next_state =  xor(next_state, mod(numel(cP2),2)==0);
                    [cP3, accLBF3, rejLBF3] = obj.greedyChangePointsEstimateRecurse(logBayesFactorThreshold, next_state, mle_cps(2), seq_end);
                    
                    accLogBayesFactors = [accLBF1, logBF, accLBF2, logBF, accLBF3];
                    rejLogBayesFactors = [rejLBF1, rejLBF2, rejLBF3];
                    changePoints = [cP1, mle_cps(1), cP2, mle_cps(2), cP3];
                end
            else
                %No changepoints predicted choose the better of free or bound
                changePoints = [];
                accLogBayesFactors = [];
                rejLogBayesFactors = logBF;
            end
            
        end
        
        function [logBF, mle_cps, mle_state0, mle_llh] = greedyLogBayesFactor(obj, seq_beg, seq_end)
            if nargin<3
                seq_end = obj.nObs;
            end
            if nargin<2
                seq_beg = 1;
            end
            prior = [0.5 0.5];
            log_prior = log(prior);
            log_pH1 = logSum(log_prior+[obj.freeLLH(seq_end,seq_beg),  obj.boundLLH(seq_end,seq_beg)]);
            
            [llh_1cp_0, llh_1cp_1, cps1] = obj.compute1CPLLH(seq_beg, seq_end);
            [llh_2cp_0, llh_2cp_1] = obj.compute2CPLLH(seq_beg, seq_end);
            
            log_pH2_1 = -log(numel(cps1)) -log(2) + logSum([llh_1cp_0; llh_1cp_1]);
            log_pH2_2 = -log(0.5*numel(cps1)*(numel(cps1)-1)) -log(2) + logSum([llh_2cp_0(:); llh_2cp_1(:)]);
            
            log_pH2 = -log(2) + logSum([log_pH2_1, log_pH2_2]);
            
            logBF = log_pH2 - log_pH1;
            [mle1_0, mle1_0_idx] = max(llh_1cp_0);
            [mle1_1, mle1_1_idx] = max(llh_1cp_1);
            [mle2_0, mle2_0_idx] = max(llh_2cp_0(:));
            [mle2_1, mle2_1_idx] = max(llh_2cp_1(:));
            mle_llh = max([mle1_0, mle1_1, mle2_0, mle2_1]);
            switch mle_llh
                case mle1_0
                    mle_state0 = 0;
                    mle_cps = cps1(mle1_0_idx);
                case mle1_1
                    mle_state0 = 1;
                    mle_cps = cps1(mle1_1_idx);
                case mle2_0
                    mle_state0 = 0;
                    [i,j] = ind2sub(size(llh_2cp_0),mle2_0_idx);
                    mle_cps = [seq_beg+j, seq_beg+i];
                case mle2_1
                    mle_state0 = 1;
                    [i,j] = ind2sub(size(llh_2cp_1),mle2_1_idx);
                    mle_cps = [seq_beg+j, seq_beg+i];
            end
        end
        
         
        function [logBF, mle_cps, mle_llh] = greedyLogBayesFactorFixedState(obj, state0, seq_beg, seq_end)
            if nargin<4
                seq_end = obj.nObs;
            end
            if nargin<3
                seq_beg = 1;
            end
            if seq_beg+1==seq_end;
                logBF=-inf;
                mle_cps=[];
                if state0
                    mle_llh = obj.boundLLH(seq_end, seq_beg);
                else
                    mle_llh = obj.freeLLH(seq_end, seq_beg);
                end
                return
            end
            if state0
                log_pH1 = obj.boundLLH(seq_end,seq_beg);
                [~, llh_1cp, cps1] = obj.compute1CPLLH(seq_beg, seq_end);
                [~, llh_2cp] = obj.compute2CPLLH(seq_beg, seq_end);
            else
                log_pH1 = obj.freeLLH(seq_end,seq_beg);
                [llh_1cp, ~, cps1] = obj.compute1CPLLH(seq_beg, seq_end);
                [llh_2cp, ~] = obj.compute2CPLLH(seq_beg, seq_end);
            end
            
            log_pH2_1 = -log(numel(cps1))+logSum(llh_1cp);
            log_pH2_2 = -log(0.5*numel(cps1)*(numel(cps1)-1))+logSum(llh_2cp(:));            
            log_pH2 = -log(2) + logSum([log_pH2_1, log_pH2_2]);
            logBF = log_pH2 - log_pH1;
            
            [mle1, mle1_idx] = max(llh_1cp);
            [mle2, mle2_idx] = max(llh_2cp(:));
            mle_llh = max([mle1, mle2]);
            if mle1>mle2
                mle_cps = cps1(mle1_idx);
            else
                [i,j] = ind2sub(size(llh_2cp),mle2_idx);
                mle_cps = [seq_beg+j, seq_beg+i];
            end
        end
        
        function [state0, cps, llh] = uphillSearch(obj, seq_beg, seq_end)
            if nargin<3
                seq_end = obj.nObs;
            end
            if nargin<2
                seq_beg = 1;
            end
            state0 = obj.freeLLH(seq_end,seq_beg) < obj.boundLLH(seq_end,seq_beg);
            cps=[];
            if state0
                llh=obj.boundLLH(seq_end,seq_beg);
            else
                llh=obj.freeLLH(seq_end,seq_beg);
            end
            niters=3;
            for K=1:niters
                %try removing any neighboring cps
                if numel(cps)>2
                    removal_llhs=zeros(numel(cps)-1);
                    for n=1:numel(cps)-1
                        test_cps = setdiff(cps, cps(n:n+1));
                        removal_llhs(n) = obj.computeSequenceLikelihood(state0, test_cps, seq_beg, seq_end);                        
                    end
                    [mx_removal_llh, mx_removal_llh_idx] = max(removal_llhs);
                    if mx_removal_llh > llh % removal will increase LLH
                        removal_cps=cps(mx_removal_llh_idx:mx_removal_llh_idx+1);
                        fprintf('[Iter:%i] Removing cps[LLH: %g->%g]: %s\n',K, llh,  mx_removal_llh, num2str(removal_cps));
                        cps = setdiff(cps, removal_cps);
                    end                        
                end
                %try propsing a single new change point that may remove the next point
                Npts = seq_end-seq_beg-1; %number of possible change points
                next_llh = zeros(Npts,1);
                next_cps = cell(Npts,1);
                for n=seq_beg+1:seq_end-1
                    idx = find(n==cps, 1, 'first');
                    if isempty(idx) %n is not already a cp
                        idx = find(cps>n, 1, 'first');
                        if isempty(idx) %there is no larger cp so add it
                            next_cps{n-seq_beg} = [cps, n];
                        else % there is a larger cp so remove that one and propose this new one
                            next_cps{n-seq_beg} = cps;
                            next_cps{n-seq_beg}(idx) = n;                                
                        end
                    else % n is a cp.  propose deleting it?
                        next_cps{n-seq_beg} = cps;
                        next_cps{n-seq_beg}(idx) = [];                        
                    end
                    next_llh(n-seq_beg) = obj.computeSequenceLikelihood(state0, next_cps{n-seq_beg}, seq_beg, seq_end);
                end
                [mx_next_llh, mx_next_llh_idx] = max(next_llh);
                if mx_next_llh > llh %propose new cps
                    fprintf('[Iter:%i] Next CP Insert Success [LLH: %g->%g]: %i\n', K, llh, mx_next_llh, mx_next_llh_idx+seq_beg);
                    cps = next_cps{mx_next_llh_idx};
                    llh = mx_next_llh;
                end
            end
        end

        %% Sanity checking code
        function LLH = computeFreePairLikelihoodDEstimator(obj, pair)
            % compute the free LLH matrix for all segments i:j. Using the DEstimator package
            % This is slower and only for testing/verification of computeFreePairLikelihood() which is
            % much faster.
            PairAnalysis.checkPairFormat(pair);
            N = size(pair,1);
            LLH = zeros(N,N);
            for i=1:N
                for j=i+1:N
                    dest_a = DEstimator(pair(i:j,[2,3]), pair(i:j,1), pair(i:j,[10,11]).^2, obj.exposureT);
                    dest_b = DEstimator(pair(i:j,[6,7]), pair(i:j,1), pair(i:j,[14,15]).^2, obj.exposureT);
                    LLH(j,i) = dest_a.LLH(obj.D_A) + dest_b.LLH(obj.D_B);                    
                end
            end
        end
        
        function LLH = compute1JumpBoundLLH(obj)
            % This should be exactly equivalent to diag(obj.boundLLH,-1) as produced by C++
            N = size(obj.pair,1); 
            LLH = zeros(1,N-1);
            dt = diff(obj.pair(:,1));
            [vMax, ~] = obj.computeVariance(obj.D_AB, dt, obj.pair(:,10), obj.exposureT, 1e-8);
            [vMay, ~] = obj.computeVariance(obj.D_AB, dt, obj.pair(:,11), obj.exposureT, 1e-8);
            [vMbx, ~] = obj.computeVariance(obj.D_AB, dt, obj.pair(:,14), obj.exposureT, 1e-8);
            [vMby, vD] = obj.computeVariance(obj.D_AB, dt, obj.pair(:,15), obj.exposureT, 1e-8);
            
            for n=1:N-1
                seg = [n,n+1];
                %xs
                obs_a = obj.pair(seg,2);
                obs_b = obj.pair(seg,6);
                var_a = vMax(seg);
                var_b = vMbx(seg);
                beta = obj.rho^2 + var_b;
                gamma = beta + var_a;
                kappa = (obs_a.*beta + obs_b.*var_a) ./ gamma;
                zeta = beta.*var_a ./ gamma;
                mu = kappa(1);
                eta = zeta(1) + vD(n); %eta1
                alpha = zeta(2) + eta; %alpha2                
                LLH(n) = LLH(n) - 0.5*( 2*log(2*pi) + sum(log(gamma)) +  sum((obs_a-obs_b).^2./gamma) + ...
                                          log(2*pi) + log(alpha) + (mu - kappa(2))^2/alpha);
                                   
                %ys
                obs_a = obj.pair(seg,3);
                obs_b = obj.pair(seg,7);
                var_a = vMay(seg);
                var_b = vMby(seg);
                beta = obj.rho^2 + var_b;
                gamma = beta + var_a;
                kappa = (obs_a.*beta + obs_b.*var_a) ./ gamma;
                zeta = beta.*var_a ./ gamma;
                mu = kappa(1);
                eta = zeta(1) + vD(n); %eta1
                alpha = zeta(2) + eta; %alpha2                
                LLH(n) = LLH(n) -0.5*( 2*log(2*pi) + sum(log(gamma)) +  sum((obs_a-obs_b).^2./gamma) + ...
                                       log(2*pi) + log(alpha) + (mu - kappa(2))^2/alpha);

            end
        end
        
        
    end
    
    methods (Static=true)
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
    end
end
