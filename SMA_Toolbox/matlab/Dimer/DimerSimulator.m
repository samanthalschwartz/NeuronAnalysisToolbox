classdef DimerSimulator < handle
    % The DimerSimulator generates undelying trajectories at a high temporal resolution.  These
    % trajectories are stored in simTrajectories property.  Then it can be used to create PairAnalysis
    % objects by simulating the observational error.
    properties (Constant=true)
        LOCALIZATION_ERRORS = {'OU',... % OU process sampled errors
                               'Fixed',... % Fixed a particular sigma
                               'Normal',... %normally distributed (in absolute value).
                               'Gamma'}; %gamma distributed
        DIMER_ORIENTATION = {'RandomRigid',... %random on every step
                             'Normal',... %Normally distributed with mean 0
                             'FixedRigid',... %preserve initial orientation
                             'OU'}; %drifts with the Ornstein-Uhlenbeck process. 
        INTERACTION_TYPE = {...
                        'Bound', ... % Always bound
                        'Flyby', ... %Non-interacting nearby routes at random time (mimics expexted "free" pair)
                        'FixedFlyby', ... %Non-interacting crossing routes at time T/2 exactly.
                        'Fixed', ... % A single fixed duration
                        'Exponential', ... %Random single binding period at random position in tracjetory
                        'ReactionSim'}; % Brownian reaction simulation selected interacting pairs
    end
    
    properties
        simDeltaT = 1e-3; % The simulation time step should be at least as fast as the fastet frame rate to simulate
        simNSteps; %Number of simulation steps

        % Changeable Parameters controlling the dimer model
        DA=0.050; % micron^2/s Diffusion constant of molecule A
        DB=0.050; % micron^2/s Diffusion constant of molecule B
        DAB=0.035; % micron^2/s Diffusion constant of molecule AB 
        dimerOrientation = 'OU';
        dimerRigidLength = 0.050; %[microns] Physical distance between flurophore centers
        dimerLengthStddev = 0.025; %[microns] [Used for: 'Normal' and 'OU' Orientation] Standard deviation of seperation distance
        dimerRelaxationRate = 1e1 ; %[1/s] Relaxation rate for OU model
        dimerOrientationDriftRate = 10*2*pi; %rad/s standard deviation of drift rate for angular diffusion

        % Changeable Parameters controlling the localization error model
        locError = 'OU'; %Obersrvational error from localization process.  We assume symmetric X&Y for now.
        locErrorMean = 0.015; % microns
        locErrorStddev = 0.005; %microns [used for: 'Normal']
        locErrorGammaShape = 6; %shape parameter for gamma function [used for: 'Gamms']
        overlayErrorStddev = 0.005; % micron.  Observational error from overlay (applied to particle B)
        
        interactionType = 'Fixed';
        interactionFixedTime = 0.8; % s [used f
        interactionKoffRate = 1e-1; %s^-1 [used for: interactionType='Exponential' and 'ReactionSim']
        %Reaction simulator control comparison
        reactionSimBindRadius = 0.050;
        reactionSimBindProb = 0.30;
        reactionSimUnbindRadius = 0.100;        
    end

    properties (Dependent=true)
        Ntrajectories; % Number of simulated trajectories sotored
        simMaxTime; %Maximum simulated time
    end

    methods
        function obj = DimerSimulator(nSteps, deltaT)
            obj.simNSteps = nSteps;
            obj.simDeltaT = deltaT;
            %Stimate good values for common parameters
            obj.interactionFixedTime = obj.simMaxTime/4;
            obj.interactionKoffRate = 4/obj.simMaxTime;

        end

        function obs = observeTrajectories(obj, Ts, frameT)
            % [in]
            %   Ts: Cell array of simulated trajectory (noisless) column fmt: [t ax ay bx by trueState]
            %   frameT: Frame time in seconds.
            % [out]
            %  obs: Cell array of observation matricies with column fmt: [t obs_ax obs_ay obs_bx obs_by Sax Say Sbx Sby trueState true_ax true_ay true_bx true_by]
            %       where trueState is 0=Free 1=Bound at observation time t as
            %             listed in corresponding observation matrix row
            %       true positions mean that obs(:,10:14) is the representation of a "particle" for the DBN
            nFrames = floor(obj.simMaxTime/frameT);            
            obs = cell(numel(Ts),1);
            trueState = cell(numel(Ts),1);
            for k=1:numel(Ts)
                T=Ts{k};
                N=size(T,1);
                obs{k} = zeros(nFrames,14);
                trueState{k} = false(nFrames,1);
                curFrameIdx = 1; % The index of the sampled frame we are working on
                curFrameStartT = 0; %The real time when the current frame exposure starts
                curFrameEndT = frameT; %The real time when the current frame exposure will end
                curFrameStartLoc = 1; %The first localization of the current frame we will eventuall include in the average
                for n=1:N
                    if n==N || T(n,1)>curFrameEndT
                        %Finished the current frame, record it
                        %Choosing to set the state based on the start of the observation frame.
                        %this is consisetent with the idea the the start of a frame is the frameTime
                        % and ideal exposure times are 0 even when frame times are finite.
                        obs{k}(curFrameIdx,1) = curFrameStartT;
                        obs{k}(curFrameIdx,2:5) = mean(T(curFrameStartLoc:n-1,2:5));
                        trueState{k}(curFrameIdx) = T(curFrameStartLoc,6);
                        if curFrameIdx == nFrames
                            break;
                        else
                            curFrameIdx = curFrameIdx+1;
                            curFrameStartT = curFrameEndT;
                            curFrameEndT = curFrameStartT + frameT;
                            curFrameStartLoc = n;
                        end
                    end
                end
                switch obj.locError
                    case 'OU'
                        OU=OUProcess(obj.locErrorMean,obj.locErrorStddev,10*frameT);
                        ou_samp = OU.sampleGillespie(obs{k}(:,1), 2);
                        errors = ou_samp(:,[1,1,2,2]); 
                    case 'Fixed'
                        errors = ones(nFrames,4).*obj.locErrorMean;
                    case 'Normal'
                        errors = abs(randn(nFrames,4).*obj.locErrorStddev + obj.locErrorMean);
                    case 'Gamma'
                        errors = gamrnd(obj.locErrorGammaShape, obj.locErrorMean/obj.locErrorGammaShape, nFrames, 4);
                    otherwise
                        error('DimerSimulator:observeTrajectories','Unknown localization Error Method: %s', obj.locError);
                end
                obs{k}(:,11:14) = obs{k}(:,2:5); % Save true positions
                obs{k}(:,6:9) = errors;
                obs{k}(:,2:5) = obs{k}(:,2:5) + randn(nFrames,4).*errors;
                obs{k}(:,10) = trueState{k};
            end            
        end

        function Ts = simulateTrajectories(obj, N)
            %Append simulated trajectories to those already generated.
            switch obj.interactionType
                case 'Bound'
                    Ts = obj.simulateTrajectories_Bound(N);
                case 'Flyby'
                    Ts = obj.simulateTrajectories_Flyby(N);
                case 'Fixed'
                    Ts = obj.simulateTrajectories_Fixed(N);
                otherwise
                    error('DimerSimulator:simulateTrajectories','Unknown interaction simulation type: %s', obj.interactionType);
            end
        end        
    end

    methods (Access=protected)
        function Ts = simulateTrajectories_Bound(obj, N)
            % Molecules are always bound
            trueState = ones(obj.simNSteps, 1); % Mol
            times = (0:obj.simNSteps-1)'*obj.simDeltaT;
            Ts = cell(N,1);
            for n=1:N
                [Pa,Pb] = obj.simulateBoundDimer(obj.simNSteps);                
                Ts{n} = [times, Pa, Pb, trueState];
            end
        end

        function Ts = simulateTrajectories_Flyby(obj, N)
            % Molecules are never bound but they are exactly at the origin at the midpoint                    
            % Simulate diffusive motion backwards and forwards from here
            trueState = zeros(obj.simNSteps, 1); % Mol
            flybyStep = floor(obj.simNSteps/2); % The step that both have the same position (0,0)
            times = (0:obj.simNSteps-1)'*obj.simDeltaT;
            Ts = cell(N,1);
            for n=1:N
                Pa_pre = obj.diffusionSim(flybyStep,obj.DA,obj.simDeltaT);
                Pb_pre = obj.diffusionSim(flybyStep,obj.DB,obj.simDeltaT);
                Pa_post = obj.diffusionSim(obj.simNSteps-flybyStep+1,obj.DA,obj.simDeltaT);
                Pb_post = obj.diffusionSim(obj.simNSteps-flybyStep+1,obj.DB,obj.simDeltaT);
                Pa=[flipud(Pa_pre); Pa_post(2:end,:)];
                Pb=[flipud(Pb_pre); Pb_post(2:end,:)];
                Ts{n} = [times, Pa, Pb, trueState];
            end
        end

        function Ts = simulateTrajectories_Fixed(obj, N)
            % Molecules always bind at the exact same time and for the same duration controlled by
            % interactionFixedTime
            if obj.interactionFixedTime > obj.simMaxTime
                error('DimerSimulator:simulateTrajectories_Fixed', 'Fixed interaction time %g(s) too large compared to max time %g(s)',obj.interactionFixedTime, obj.simMaxTime);
            end
            NBoundSteps = ceil(obj.interactionFixedTime / obj.simDeltaT);
            NPreBoundSteps = floor((obj.simNSteps - NBoundSteps) /2);
            NPostBoundSteps = obj.simNSteps - NPreBoundSteps - NBoundSteps;
            trueState = zeros(obj.simNSteps, 1); % Mol
            trueState(NPreBoundSteps+1:NPreBoundSteps+NBoundSteps,:)=1;            
            times = (0:obj.simNSteps-1)'*obj.simDeltaT;
            Ts = cell(N,1);
            for n=1:N                
                [Pa_bound,Pb_bound] = obj.simulateBoundDimer(NBoundSteps);
                Pa_pre = obj.diffusionSim(NPreBoundSteps+1,obj.DA,obj.simDeltaT,Pa_bound(1,:));
                Pb_pre = obj.diffusionSim(NPreBoundSteps+1,obj.DB,obj.simDeltaT,Pb_bound(1,:));
                Pa_post = obj.diffusionSim(NPostBoundSteps+1,obj.DA,obj.simDeltaT,Pa_bound(end,:));
                Pb_post = obj.diffusionSim(NPostBoundSteps+1,obj.DB,obj.simDeltaT,Pb_bound(end,:));
                Pa=[flipud(Pa_pre); Pa_bound(2:end-1,:); Pa_post];
                Pb=[flipud(Pb_pre); Pb_bound(2:end-1,:); Pb_post];
                Ts{n} = [times, Pa, Pb, trueState];
            end
        end

        function [Pa, Pb] = simulateBoundDimer(obj, Nsteps)
            Pd = obj.diffusionSim(Nsteps,obj.DAB,obj.simDeltaT,[0,0]);
            switch obj.dimerOrientation
                case 'Normal'
                    Pa=Pd;
                    Pb=Pa + randn(Nsteps,2).*obj.dimerLengthStddev;
                case 'RandomRigid'
                    phi = rand(Nsteps,1)*2*pi;
                    v = [cos(phi), sin(phi)]*(obj.dimerRigidLength/2);
                    Pa = Pd-v;
                    Pb = Pd+v;
                case 'FixedRigid'
                    phi = rand(1,1)*2*pi;
                    v = repmat([cos(phi), sin(phi)]*(obj.dimerRigidLength/2),Nsteps,1);
                    Pa = Pd-v;
                    Pb = Pd+v;
                case 'OU'
                    phi = mod(cumsum([rand(1,1)*2*pi; randn(Nsteps-1,1).*(2*obj.dimerOrientationDriftRate*obj.simDeltaT)]),2*pi);
                    gamma = obj.dimerRelaxationRate;
                    sigma = obj.dimerLengthStddev;
                    dt = obj.simDeltaT;
                    mu = obj.dimerRigidLength;
                    x0 = mu;
                    ou = OUProcess(mu,sigma,gamma,x0);
                    ts = ((1:Nsteps)-1)*dt;
                    r = ou.sampleGillespie(ts,1);
                    v = [cos(phi).*(r/2), sin(phi).*(r/2)];
                    Pa = Pd-v;
                    Pb = Pd+v;
                otherwise
                    error('DimerSimulator:simulateBoundDimer','Unknown Dimer Orientation Model: %s',obj.dimerOrientation);
            end
        end
    end
    
    methods (Static=true)
        function P = diffusionSim(N,D,dT,start_pos)
            % (in) D - Diffusion constant 
            % (in) times - Vector of times size Nx1
            % (in) start_pos - Start position [X,Y]: (optional: default [0,0])
            % (out) P - postions size Nx2
            if nargin<4
                start_pos = [0,0];
            else
                start_pos = start_pos(:)';
            end
            P = cumsum([start_pos; randn(N-1,2).* sqrt(2*D*dT)]);
        end
        
        function plotDMLE(Ts,exposureT)
            % When Ts is the direct trajectory simulation, exposureT=0.
            Ts = makecell(Ts);
            if nargin<2
                exposureT = min(diff(Ts{1}(:,1)));
            end
            figure()
            xlabel('D (um^2/s)');
            ylabel('LLH');
            Ds = logspace(-5,1,1e3);
            if size(Ts,2)<=6
                TAs = cellmap(@(T) [T(:,[1,2,3]),zeros(size(T,1),2)], Ts);
                TBs = cellmap(@(T) [T(:,[1,4,5]),zeros(size(T,1),2)], Ts);
            elseif size(Ts,3)>=9
                TAs = cellmap(@(T) T(:,[1,2,3,6,7]), Ts);
                TBs = cellmap(@(T) T(:,[1,4,5,8,9]), Ts);
            end
            TDs = cellmap(@(T) [T(:,1), mean(T(:,[2,4]),2), mean(T(:,[3,5]),2)], Ts);
            LLH_A = DEstimator.computeEnsembleLLH(Ds, TAs, exposureT);
            LLH_B = DEstimator.computeEnsembleLLH(Ds, TBs, exposureT);
%             LLH_D = DEstimator.computeEnsembleLLH(Ds, TDs, exposureT);
            [mle_A, mle_A_llh] = DEstimator.computeEnsembleMLE(TAs, exposureT);
            [mle_B, mle_B_llh] = DEstimator.computeEnsembleMLE(TBs, exposureT);
%             [mle_D, mle_D_llh] = DEstimator.computeEnsembleMLE(TDs, exposureT);
            hold('on');
            plot(Ds,LLH_A,'r-','DisplayName','LLH(D): Particle A');
            plot(Ds,LLH_B,'b-','DisplayName','LLH(D): Particle B');
%             plot(Ds,LLH_D,'g-','DisplayName','LLH(D): Dimer (mean pos)');
            plot(mle_A, mle_A_llh,'rs','DisplayName',sprintf('MLE(D): Particle A = %.6g um^2/s', mle_A));
            plot(mle_B, mle_B_llh,'bs','DisplayName',sprintf('MLE(D): Particle B = %.6g um^2/s', mle_B));
%             plot(mle_D, mle_D_llh,'gs','DisplayName',sprintf('MLE(D): Dimer (mean pos) = %.6g um^2/s', mle_D));
            legend('location','best');            
            set(gca(),'XScale','log','Yscale','log');
        end
        

       function plotDistances(Traj, Obs)
            F=figure();
            F.Color=[1 1 1];
            certainty=0.95;
            obsPair = zeros(size(Obs,1), 17);
            obsPair(:,1:3) = Obs(:,1:3);
            obsPair(:,6:7) = Obs(:,4:5);
            obsPair(:,10:11) = Obs(:,6:7);
            obsPair(:,14:15) = Obs(:,8:9);


            D = PairAnalysis.makePairDists(obsPair, certainty);
            plot(D(:,1),D(:,2),'-r','Marker','.','DisplayName','Dist');
            hold on;
            true_dists = sqrt(sum( (Traj(:,[2,3]) - Traj(:,[4,5])).^2,2));
            nSteps = size(Traj,1);
            deltaT = Traj(2,1)-Traj(1,1);
            plot( (0:nSteps-1)'*deltaT,true_dists,'.','Color',[0.5,0.2,0],'MarkerSize',2.5,'LineWidth',0.5,'DisplayName','True Dist.');
            xs = [D(:,1)', fliplr(D(:,1)')];
            ys = [D(:,3)', fliplr(D(:,4)')];
            fill(xs,ys,'r','EdgeColor','none','FaceAlpha',0.4,'DisplayName',sprintf('Dist (%.1f%% conf.)',certainty*100));
            yl = ylim();
            stateSeq = logical(Traj(:,end));
            state = stateSeq(1);
            lastIdx = 1;
            for i=1:size(Traj,1)
                if state ~= stateSeq(i)
                    if state
                        xs=[Traj(lastIdx,1),Traj(i,1)+deltaT, Traj(i,1)+deltaT, Traj(lastIdx,1)];
                        ys=[0,0, yl(2), yl(2)];
                        h=fill(xs,ys, 'c','EdgeColor','None', 'FaceAlpha',0.25);
                        h.Annotation.LegendInformation.IconDisplayStyle='off';
                    end
                    lastIdx=i;
                    state = stateSeq(i);
                end
            end
            ylim(yl);
            xlabel('Time (s)');
            ylabel('Dist (um)');
            legend('Location','best');
            
        end

       function browseObservedTrajectories(Ts, Obs)
            F=figure('Position',[30,30,800,700]);
            Ts = makecell(Ts);
            Obs = makecell(Obs);
            N = numel(Ts);
            assert(numel(Obs) == numel(Ts));
            if N>1
                uicontrol('Style','slider','Min',1, 'Max',N,'Value',1,'SliderStep',[1/(N-1),1/(N-1)],...
                          'Units','normalized','Position',[0.05, 0, 0.8, 0.035], 'Callback',@slider_CB);
            end
            whitebg(F);
            F.Color = [0 0 0];
            xmin = min(cellfun(@(T) min(min(T(:,[2,4]))), [Ts; Obs]));
            xmax = max(cellfun(@(T) max(max(T(:,[2,4]))), [Ts; Obs]));
            ymin = min(cellfun(@(T) min(min(T(:,[3,5]))), [Ts; Obs]));
            ymax = max(cellfun(@(T) max(max(T(:,[3,5]))), [Ts; Obs]));
            tmin = min(cellfun(@(T) T(1,1), Ts));
            tmax = max(cellfun(@(T) T(end,1), Ts));
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            zlim([tmin,tmax]);
            grid('on'); grid('minor');
            view(-45,30);
            asp(1) = (xmax-xmin)/(ymax-ymin);
            asp(2) = 1;
            asp(3) = mean([asp(1),asp(2)]);
            aH=[];
            bH=[];
            oaH=[];
            obH=[];
            dH=[];
            dBs={};
            drawBonds=true; %Change this to control bond drawing
            plotT(Ts{1}, Obs{1});
            title(sprintf('Observed Trajectory %i',1));
            pbaspect(asp);
            xlabel('X (um)');
            ylabel('Y (um)');
            zlabel('T (s)');
            ax=gca();
            ax.TickDir='out';
            ax.Box='on';
            ax.BoxStyle='Full';
            function slider_CB(hObj,~)
                k=round(hObj.Value);
                plotT(Ts{k}, Obs{k});
                title(sprintf('Trajectory %i',k));
            end
            function plotT(T,O)
                if ishandle(aH); delete(aH); end
                if ishandle(bH); delete(bH); end
                if ishandle(oaH); delete(oaH); end
                if ishandle(obH); delete(obH); end
                if ishandle(dH); delete(dH); end
                times = T(:,[1,1]);
                otimes = O(:,[1,1]);
                aH=surface(T(:,[2,2]),T(:,[3,3]),times,'EdgeColor',[0.5,0,0]);
                hold('on');
                bH=surface(T(:,[4,4]),T(:,[5,5]),times,'EdgeColor',[0,0.3,0.8]);
                oaH=surface(O(:,[2,2]),O(:,[3,3]),otimes,'Marker','o','LineWidth',1,'MarkerEdgeColor',[.8,.2,0],...
                                                         'MarkerFaceColor',[1,0,0],'EdgeColor',[1,0,0]);
                obH=surface(O(:,[4,4]),O(:,[5,5]),otimes,'Marker','o','LineWidth',1,'MarkerEdgeColor',[0,.2,.8],...
                                                         'MarkerFaceColor',[0,0,1],'EdgeColor',[0,0,1]);
                bndTs=find(T(:,6));
                if drawBonds
                    for i=1:numel(dBs)
                        if ishandle(dBs{i}); delete(dBs{i}); end
                    end
                    dBs=cell(numel(bndTs),1);
                    for i=1:numel(bndTs)
                        b = T(bndTs(i),:);                        
                        dBs{i}=surface(b([2 2;4 4]),b([3 3;5 5]),b([1 1;1 1]),'LineWidth',1,'EdgeColor',[0,1,0],'EdgeAlpha',0.55); 
                    end
                end
            end            
        end        

        function browseTrajectories(Ts)
            F=figure('Position',[30,30,800,700]);
            Ts = makecell(Ts);
            N = numel(Ts);
            if N>1
                uicontrol('Style','slider','Min',1, 'Max',N,'Value',1,'SliderStep',[1/(N-1),1/(N-1)],...
                    'Units','normalized','Position',[0.05, 0, 0.8, 0.035], 'Callback',@slider_CB);
            end
            whitebg(F);
            F.Color=[0 0 0];
            xmin = min(cellfun(@(T) min(min(T(:,[2,4]))), Ts));
            xmax = max(cellfun(@(T) max(max(T(:,[2,4]))), Ts));
            ymin = min(cellfun(@(T) min(min(T(:,[3,5]))), Ts));
            ymax = max(cellfun(@(T) max(max(T(:,[3,5]))), Ts));
            tmin = min(cellfun(@(T) min(T(:,1)), Ts));
            tmax = max(cellfun(@(T) max(T(:,1)), Ts));
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            zlim([tmin,tmax]);
            grid('on');
            grid('minor');
            box('on');
            view(-45,30);
            asp(1) = (xmax-xmin)/(ymax-ymin);
            asp(2) = 1;
            asp(3) = mean([asp(1),asp(2)]);
            aH=[];
            bH=[];
            dH=[];
            dBs={};
            drawBonds=true; %Change this to control bond drawing
            plotT(Ts{1});
            title(sprintf('Trajectory %i',1));
            pbaspect(asp);
            xlabel('X (um)');
            ylabel('Y (um)');
            zlabel('T (s)');
            ax=gca();
            ax.TickDir='out';
            ax.BoxStyle='Full';
            function slider_CB(hObj,~)
                k=round(hObj.Value);
                plotT(Ts{k});
                title(sprintf('Trajectory %i',k));
            end
            function plotT(T)
                if ishandle(aH); delete(aH); end
                if ishandle(bH); delete(bH); end
                if ishandle(dH); delete(dH); end
                times = T(:,[1,1]);
                aH=surface(T(:,[2,2]),T(:,[3,3]),times,'EdgeColor',[1,0,0]);
                hold('on');
                bH=surface(T(:,[4,4]),T(:,[5,5]),times,'EdgeColor',[0,0,1]);
                dX = mean(T(:,[2,4]),2);
                dY = mean(T(:,[3,5]),2);                
%                 dH=surface([dX,dX], [dY,dY], times,'EdgeColor',[0,1,0],'EdgeAlpha','Flat','AlphaData',T(:,[6,6]));
                bndTs=find(T(:,6));
                if drawBonds
                    for i=1:numel(dBs)
                        if ishandle(dBs{i}); delete(dBs{i}); end
                    end
                    dBs=cell(numel(bndTs),1);
                    for i=1:numel(bndTs)
                        b = T(bndTs(i),:);                        
                        dBs{i}=surface(b([2 2;4 4]),b([3 3;5 5]),b([1 1;1 1]),'EdgeColor',[0,1,0],'EdgeAlpha',0.45); 
                    end
                end
            end            
        end        
    end
    
    methods
        % Dependent property accessors
        function val = get.Ntrajectories(obj)
            val =  numel(obj.simTrajectories);
        end
        function val = get.simMaxTime(obj)
            val =  obj.simDeltaT * obj.simNSteps;
        end
    end       
end

 
