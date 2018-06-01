classdef GeometricDist_MLEandMCMC < handle
% This class fits a vector of positive intergers, representing ontime duration
% in frames, to a Geometric Distribution model with N=1,2 or 3 components.
% It uses either Maximum Likelihood or Markov Chain Monte Carlo approaches
% for parameter estimation. 
%
% Inputs to function calls are
% '1comp': Fits all data to the same 1 component model (single P value representing offrate for all data)
% '2comp': Fits all data to the same 2 component model (single A1 value representing fraction and P1,P2 value representing offrate for all data)
% '3comp': Fits all data to the same 3 component model (single A1,A2,value representing fraction P1,P2,P3 value representing offrate for all data)
% (): No input defaults to 1 component independent model
% constrained probability boolean vector: [1],[1 0],[0 1 1] where length of vector indicates number of fit components
% Also use GeometricDist_MLEandMCMC.MODELS to see a list of all the possible fitting models.
%
% Example starter code:
%--------
% gd = GeometricDist_MLEandMCMC();
% gd.setOffset(5);
% gd.setframerate(10);
% gd.K = 5;
% gd.simulate(.8,[2.581,0.625],10000);
% model = '2-comp';
% [tmle, errtmle, fmle] = gd.fit_MLE(model);
%-----
%
% Some more information:
% gd = GeometricDist_MLEandMCMC(); %Initial the class object
% gd.Offset = 5; %Set the offset value (in frames). % all input tracks are longer than this value. ie: Offset = 5 means minimun input track length is 6
% gd.setframerate(10); %Input the framerate in frames/second. This is an optional parameter. If this is not set, the default output is in units of frames.
% %--- when inputing your own data-----
% gd.D = {Data1,Data2}; %Cell array of input data. Each data set is a
% vector of positive integers representing tracklengths in frames.
% %--- to simulate a dataset ----
% gd.K = 2; %number of datasets to simulate. 
% gd.simulate([.3],[2.5,0.6],60000); Input alpha, p1,p2 values. If framerate is set, then p1,p2 values assumed to be in units of 1/sec. Otherwise assumed to be in units of 1/frames.
% %-------
% gd.setConditionNames({'TestCondition'}); %This allows you to give each
% condition a name for displaying.
% % Analysis methods:
% [tmle, errtmle, fmle] = gd.fit_MLE(model); % This uses maximum likelihood to fit the parameters.
% [t, err_t,f] = gd.fit_MCMC(model); % This uses MCMC to fit the parameters.
% % To get more information about the MCMC use the following plotting
% functions:
% [samples, samplesLLH, proposed_samples, proposed_samplesLLH]  = gd.MCMC(model);
% gd.plotMCMC(model,samples) % returns the distribution of values for each parameter
% gd.plotMCMCchain(model,samples, proposed_samples) % returns the full chain including rejected points for each parameter
% gd.plotProfileLLH(model)

% [ModelProbabilities,selectedModel,LogEvidence_Array,bestTheta,allsamples] = evidence_Test(obj,models)
% % -- this method calculates the evidence for each model given the data
% and compares to determine the model with the maximum evidence.
%
% ----- S.L. Schwartz, M.J. Olah, K.A. Lidke. University of New Mexico 2015

    properties (Constant=true, Hidden=true)
        COLORS = {[1,0,0],[0,1,0],[0,0,1],[1,0.5,0],[0.5,0,1],[0,0.5,1],...
            [0, 0, 0], [0.5, 0, 0], [0.3, 0.7, 0.2]};
    end
    properties (Constant=true)
        MODELS = {'1-comp','1-comp ind',...
            '2-comp','2-comp ind','2-comp constrain-P1','2-comp constrain-P2','2-comp constrain-P1P2',...
            '3-comp','3-comp ind','3-comp constrain-P1','3-comp constrain-P2','3-comp constrain-P3',...
            '3-comp constrain-P1P2','3-comp constrain-P1P3','3-comp constrain-P2P3','3-comp constrain-P1P2P3'};
    end
    properties
        K; % number of conditions. -- need to change and make this private. just need to update in simulate function to take in variable
        MCMCSampleNum = 10000;
        MCMCBurnFrac = 0.25;
        MCMCSigma = 0.05;
        plotDiagnostics = false; % set this to display info about optimizations (both for MLE and MCMC approach)
        Verbose = false; % set this to display progress and MLE/SSE values on command line
    end  
    properties (GetAccess = 'public', SetAccess = 'private')
        % These properties must be set using the 'setData' method
        Dinput; % this is the input data: cell array of vectors of track lengths. integer values >=1. D{} = n x 1 vector, where n is number of datapoints
        D; % this is the data used for analysis with Offset applied: cell array of vectors of track lengths. integer values >=1. D{} = n x 1 vector, where n is number of datapoints
        Ntot; % Total number of samples
        N; % vector size:k Number of samples per-condition
        Stot; % Sum of the total data
        S; %vector size:k Sum of the data for each condition
%         K; % number of conditions. -- need to update simulate function
%         before making this property private
        % This property must be set using the 'setConditionNames' method
        Offset = 0; % all tracks are longer than this value (value in frames). ie: Offset = 5 means minimun track length is 6
        ConditionNames; %cell array of condition names associated with obj.D, used for plotting etc numel(ConditionNames) = numel(D);
        framerate = 0; % frame rate (frames/second): This is used to convert from frames to generate offrate output in sec^-1. If this is set then it assumes simulate() input is given in sec^-1 and uses sec^-1 for plotting 
        rateunits = 'Frames';
        rateunits_html = 'Frames'
    end
    methods
        function obj = GeometricDist_MLEandMCMC(varargin)
            if nargin>0
                obj.setData(varargin{:});
            end
        end
        function setData(obj,Dinput)
            Dinput = makecell(Dinput);
            obj.K = numel(Dinput);
            Dinput = cellmap(@(d) d(:),Dinput); %Make each D value a column vector
            obj.Dinput = Dinput;
            allD = vertcat(Dinput{:}); %all the data in one column
            assert(all(allD == fix(allD)),'make sure all data is an integer'); %make sure all data is an integer
            assert(all(allD>0),'make sure all data is non-zero'); % make sure all data is non-zero
            assert(all(isfinite(allD)),'make sure all data is finite'); %make sure all data is finite
            D = Dinput; % now D is a copy of Dinput
            % don't include any data that is shorter than Offset
            if ~isempty(obj.Offset)
                for ii = 1:numel(D)
                    D{ii}(D{ii}<=obj.Offset) = [];
                end
            end
            %At this point data is OK
            obj.D = D;
            obj.N = cellfun(@numel,obj.D)';
            obj.Ntot = sum(obj.N);
            obj.S = cellfun(@sum,obj.D)';
            obj.Stot = sum(obj.S);
        end 
        function clearData(obj)
            obj.D = [];
        end
        function setOffset(obj,offset)
            obj.Offset = offset;
            if ~isempty(obj.Dinput)                
                obj.setData(obj.Dinput);
            end
        end
        function setConditionNames(obj,ConditionNames)
           % this associates data plots with ConditionNames
           % numel(ConditionNames) == numel(obj.D)
           % note - any other changes to plotting output naming should be done
           % inside getModelInfo where output.modelParamDesc is set
           if isempty(obj.D)
               display('Need to add some data using obj.setData first');
           elseif isempty(ConditionNames)
               obj.ConditionNames = cellmap(@(k) sprintf('Condition #%i',k),1:obj.K);
           else
            assert(numel(ConditionNames) == numel(obj.D),'Size of ConditionNames needs to match size of obj.D');
            obj.ConditionNames = ConditionNames;
           end
        end
        function setframerate(obj,fr)
           % this sets the frame rate and the appropriate frame rate units
           % assumes frame rate 'fr' is given in frames/second
           % -- if you set the framerate but want to convert back to
           % reporting GMM parameters in frames set fr = 0
           assert(isnumeric(fr),'frame rate input must be a number');
           assert(numel(fr) == 1,'there can only be one frame rate input given per object');
           assert(fr>=0,'frame rate input must be >  0');
           obj.framerate = fr;
           if fr==0  
               obj.rateunits = 'Frames';
               obj.rateunits_html = 'Frames';
           else
               obj.rateunits = 'Offrate (sec^{-1})';
               obj.rateunits_html = '<html>Offrate (sec<sup>-1</sup>)</html>'; 
           end
        end
        %% Functions for Model Specifics
        function output = getModelInfo(obj,modelName,fixedP_vals)
            % Now generate output based on the model-----------------------------------
            %--- Outputs for all model examples -------
            % output.modelName = '1-comp';
            % output.Nalphas = [];
            % output.Ncomp = [];
            % ouput.Nparam = [];
            % output.modelParamDesc = [];
            %
            % output.problem.x0 = []; Must be a column vector for fmincon
            % call
            % output.problem.lb = [];
            % output.problem.ub = [];
            % output.problem.Aineq = [];
            % output.problem.Bineq = [];
            %
            % output.MCMC.thetainit = [];
            % output.MCMC.fixedProb = [];
            %----------------------------------------
            if isempty(obj.ConditionNames)
                paramNames = cell(1,obj.K);
                for nn = 1:obj.K
                    paramNames{nn} = ' ';
                end
            else
                paramNames = obj.ConditionNames;
            end
            output.modelName = modelName;
            if nargin<3
                output.fixedP_vals = [];
            else
                output.fixedP_vals = fixedP_vals;
            end
            switch modelName
                case '1-comp'                   %(1)
                    output.Nalphas = 0;
                    output.Np1s = 1;
                    output.Np2s = 0;
                    output.Np3s = 0;
                    output.Ncomp = 1;
                    output.Nparam = 1;
                    output.modelParamDesc = ...
                        cellmap(@(k) [sprintf('P-cond%i',k) '-' paramNames{k}],1:obj.K);
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '1-comp ind'               %(2)
                    output.Nalphas = 0;
                    output.Np1s = obj.K;
                    output.Np2s = 0;
                    output.Np3s = 0;
                    output.Ncomp = 1;
                    output.Nparam = obj.K;
                    output.modelParamDesc = ...
                        cellmap(@(k) [sprintf('P-cond%i',k) '-' paramNames{k}],1:obj.K);
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '2-comp'                   %(3)
                    output.Nalphas = 1;
                    output.Np1s = 1;
                    output.Np2s = 1;
                    output.Np3s = 0;
                    output.Ncomp = 2;
                    output.Nparam = 3;
                    output.modelParamDesc = ...
                        {'Alpha', 'P1','P2'};
                    p1= 2*obj.Ntot / (obj.Ntot+obj.Stot);
                    p2= (obj.Ntot/2) / (obj.Ntot/2+obj.Stot);
                    theta_init=[0.5;p1(:);p2(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    output.problem.Aineq = [0,-1,1];
                    output.problem.Bineq = -10^-3;
                case '2-comp ind'               %(4)
                    output.Nalphas = obj.K;
                    output.Np1s = obj.K;
                    output.Np2s = obj.K;
                    output.Np3s = 0;
                    output.Ncomp = 2;
                    output.Nparam = 3*obj.K;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P2-cond%i',k) '-' paramNames{k}],1:obj.K)];
                    p1s= (2*obj.N) ./ (obj.N+obj.S);
                    p2s= (obj.N/2) ./ (obj.N/2+obj.S);
                    alphas = obj.S./obj.Stot;
                    theta_init=[alphas(:);p1s(:);p2s(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    A=zeros(obj.K, output.Nparam);
                    b=repmat(-1e-3,obj.K,1);
                    for k=1:obj.K
                        A(k,obj.K+k) = -1;
                        A(k,2*obj.K+k) = 1;
                    end
                    output.problem.Aineq = A;
                    output.problem.Bineq = b;
                case '2-comp constrain-P1'      %(5)
                    output.Nalphas = obj.K;
                    output.Np1s = 1;
                    output.Np2s = obj.K;
                    output.Np3s = 0;
                    output.Ncomp = 2;
                    output.Nparam = 2*obj.K+1;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P1',...
                        cellmap(@(k) [sprintf('P2-cond%i',k) '-' paramNames{k}],1:obj.K)];
                    p1 = 2*obj.Ntot / (obj.Ntot+obj.Stot);
                    p2s= (obj.N/2) ./ (obj.N/2+obj.S);
                    alphas = obj.S./obj.Stot;
                    theta_init=[alphas(:);p1(:);p2s(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    A=zeros(obj.K, output.Nparam);
                    b=zeros(obj.K,1);
                    for k=1:obj.K
                        A(k,obj.K+1) = -1;
                        A(k,obj.K+1+k) = 1;
                    end
                    output.problem.Aineq = A;
                    output.problem.Bineq = b;
                case '2-comp constrain-P2'      %(6)
                    output.Nalphas = obj.K;
                    output.Np1s = obj.K;
                    output.Np2s = 1;
                    output.Np3s = 0;
                    output.Ncomp = 2;
                    output.Nparam = 2*obj.K+1;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P2'];
                    p1s= (2*obj.N) ./ (obj.N+obj.S);
                    p2= (obj.Ntot/2) / (obj.Ntot/2+obj.Stot);%same value repeated obj.K times
                    alphas = obj.S./obj.Stot;
                    theta_init=[alphas(:);p1s(:);p2(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    A=zeros(obj.K, output.Nparam);
                    b=zeros(obj.K,1);
                    for k=1:obj.K
                        A(k,obj.K+k) = -1;
                        A(k,2*obj.K+1) = 1;
                    end
                    output.problem.Aineq = A;
                    output.problem.Bineq = b;
                case '2-comp constrain-P1P2'    %(7)
                    output.Nalphas = obj.K;
                    output.Np1s = 1;
                    output.Np2s = 1;
                    output.Np3s = 0;
                    output.Ncomp = 2;
                    output.Nparam = obj.K+2;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha-cond%i',k) '-' paramNames{k}],1:obj.K), 'P1','P2'];
                    p1= 2*obj.Ntot / (obj.Ntot+obj.Stot);
                    p2= (obj.Ntot/2) / (obj.Ntot/2+obj.Stot);
                    alphas = obj.S./obj.Stot;
                    theta_init=[alphas(:);p1(:);p2(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    output.problem.Aineq = [zeros(1, output.Nparam-2), -1, 1];
                    output.problem.Bineq = -1e-5;
                case '3-comp'                   %(8)
                    output.Nalphas = 2;
                    output.Np1s = 1;
                    output.Np2s = 1;
                    output.Np3s = 1;
                    output.Ncomp = 3;
                    output.Nparam = 5;
                    output.modelParamDesc = ...
                        {'Alpha1', 'Alpha2', 'P1', 'P2', 'P3'};
                    p1= 3*obj.Ntot / (2*obj.Ntot+obj.Stot);
                    p2= 2*obj.Ntot / (obj.Ntot+obj.Stot);
                    p3= (obj.Ntot/3) / (obj.Ntot/3+obj.Stot);
                    assert(p1>p2 && p2>p3,'Ps must be in decreasing order(largest p1>p2>p3');
                    theta_init=[0.3; 0.3;p1(:);p2(:);p3(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    %Linear inequality constraints force:
                    %  (1)  alpha1+alpha2<1
                    %  (2)  -theta1+theta2<0
                    %  (3)  -theta2+theta3<0
                    output.problem.Aineq = [1 1 0 0 0; 0 0 -1 1 0; 0 0 0 -1 1];
                    output.problem.Bineq = [1-eps;-eps;-eps];
                case '3-comp ind'               %(9)
                    output.Nalphas = obj.K*2;
                    output.Np1s = obj.K;
                    output.Np2s = obj.K;
                    output.Np3s = obj.K;
                    output.Ncomp = 3;
                    output.Nparam = 5*obj.K;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P3-cond%i',k)  '-' paramNames{k}],1:obj.K)];
                    p1= 3*obj.N ./ (2*obj.N+obj.S);
                    p2= 2*obj.N ./ (obj.N+obj.S);
                    p3= (obj.N/3) ./ (obj.N/3+obj.S);
                    alpha = 0.3*ones(2*obj.K,1);
                    assert(all(p1>p2) && all(p2>p3),'Ps must be in decreasing order(largest p1>p2>p3');
                    theta_init=[alpha(:); p1(:);p2(:);p3(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    %Linear inequality constraints force:
                    %  (K)  alpha1+alpha2<1
                    %  (K)  -theta1+theta2<0
                    %  (K)  -theta2+theta3<0
                    A=zeros(3*obj.K, output.Nparam);
                    b=zeros(3*obj.K,1);
                    for k=1:obj.K %implement the alpha1+alpha2<1 constraints
                        A(k,k) = 1;
                        A(k,k+obj.K) = 1;
                        b(k)=1-eps;
                    end
                    for k=1:obj.K %implement the p1>p2 constraints
                        A(obj.K+k, 2*obj.K+k) = -1; %p1
                        A(obj.K+k, 3*obj.K+k) = 1;  %p2
                        b(obj.K+k)=-eps;
                    end
                    for k=1:obj.K %implement the p2>p3 constraints
                        A(2*obj.K+k, 3*obj.K+k) = -1; %p2
                        A(2*obj.K+k, 4*obj.K+k) = 1;  %p3
                        b(2*obj.K+k)=-eps;
                    end
                    output.problem.Aineq = A;
                    output.problem.Bineq = b;
                case '3-comp constrain-P1'      %(10)
                    output.Nalphas = obj.K*2;
                    output.Np1s = 1;
                    output.Np2s = obj.K;
                    output.Np3s = obj.K;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + obj.K*2 + 1;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P1',...
                        cellmap(@(k) [sprintf('P2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P3-cond%i',k) '-' paramNames{k}],1:obj.K)];
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '3-comp constrain-P2'      %(11)
                    output.Nalphas = obj.K*2;
                    output.Np1s = obj.K;
                    output.Np2s = 1;
                    output.Np3s = obj.K;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + obj.K*2 + 1;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P2',...
                        cellmap(@(k) [sprintf('P3-cond%i',k) '-' paramNames{k}],1:obj.K)];
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '3-comp constrain-P3'      %(12)
                    output.Nalphas = obj.K*2;
                    output.Np1s = obj.K;
                    output.Np2s = obj.K;
                    output.Np3s = 1;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + obj.K*2 + 1;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P3'];
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '3-comp constrain-P1P2'    %(13)
                    output.Nalphas = obj.K*2;
                    output.Np1s = 1;
                    output.Np2s = 1;
                    output.Np3s = obj.K;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + obj.K + 2;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P1',...
                        'P2',...
                        cellmap(@(k) [sprintf('P3-cond%i',k) '-' paramNames{k}],1:obj.K)];
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '3-comp constrain-P1P3'    %(14)
                    output.Nalphas = obj.K*2;
                    output.Np1s = 1;
                    output.Np2s = obj.K;
                    output.Np3s = 1;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + obj.K + 2;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P1',...
                        cellmap(@(k) [sprintf('P2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P3'];
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '3-comp constrain-P2P3'    %(15)
                    output.Nalphas = obj.K*2;
                    output.Np1s = obj.K;
                    output.Np2s = 1;
                    output.Np3s = 1;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + obj.K + 2;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('P1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P2',...
                        'P3'];
                    output.problem.x0 = [];
                    output.problem.lb = [];
                    output.problem.ub = [];
                    output.problem.Aineq = [];
                    output.problem.Bineq = [];
                case '3-comp constrain-P1P2P3'  %(16)
                    output.Nalphas = obj.K*2;
                    output.Np1s = 1;
                    output.Np2s = 1;
                    output.Np3s = 1;
                    output.Ncomp = 3;
                    output.Nparam = obj.K*2 + 3;
                    output.modelParamDesc = ...
                        [ cellmap(@(k) [sprintf('Alpha1-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Alpha2-cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'P1','P2','P3'];
                    p1= 3*obj.Ntot / (2*obj.Ntot+obj.Stot);
                    p2= 2*obj.Ntot / (obj.Ntot+obj.Stot);
                    p3= (obj.Ntot/3) / (obj.Ntot/3+obj.Stot);
                    alpha = 0.3*ones(2*obj.K,1);
                    assert(p1>p2 && p2>p3,'Ps must be in decreasing order(largest p1>p2>p3');
                    theta_init=[alpha(:);p1(:);p2(:);p3(:)];
                    output.problem.x0 = theta_init;
                    output.problem.lb = zeros(output.Nparam,1)+eps;
                    output.problem.ub = ones(output.Nparam,1)-eps;
                    %Linear inequality constraints force:
                    %  (K)  alpha1+alpha2<1
                    %  (1)  -theta1+theta2<0
                    %  (1)  -theta2+theta3<0
                    A=zeros(obj.K+2, output.Nparam);
                    b=zeros(obj.K+2,1);
                    for k=1:obj.K %implement the alpha1+alpha2<1 constraints
                        A(k,k) = 1;
                        A(k,k+obj.K) = 1;
                        b(k)=1-eps;
                    end
                    %implement the p1>p2 constraint
                    A(end-1,end-2) = -1;
                    A(end-1,end-1) = 1;
                    b(end-1) = -eps;
                    %implement the p2>p3 constraint
                    A(end,end-1) = -1;
                    A(end,end) = 1;
                    b(end) = -eps;
                    output.problem.Aineq = A;
                    output.problem.Bineq = b;
            end
        end
        function [alpha1, alpha2, alpha3, p1, p2, p3] = getValsFromTheta(obj,theta,model)
            % ----- Possible Models -------
            %     '1-comp'
            %     '1-comp ind'
            %     '2-comp'
            %     '2-comp ind'
            %     '2-comp constrain-P1'
            %     '2-comp constrain-P2'
            %     '2-comp constrain-P1P2'
            %     '3-comp'
            %     '3-comp ind'
            %     '3-comp constrain-P1'
            %     '3-comp constrain-P2'
            %     '3-comp constrain-P3'
            %     '3-comp constrain-P1P2'
            %     '3-comp constrain-P1P3'
            %     '3-comp constrain-P2P3'
            %     '3-comp constrain-P1P2P3'
            % -----------------------------
            switch model
                case '1-comp'                   %(1)
                    p1= theta;
                    p2 = [];
                    alpha1 = [];
                    alpha2 = [];
                    alpha3 = [];
                    p3 = [];
                case '1-comp ind'               %(2)
                    p1= theta;
                    p2 = [];
                    alpha1 = [];
                    alpha2 = [];
                    alpha3 = [];
                    p3 = [];
                case '2-comp'                   %(3)
                    alpha1 = theta(:,1);
                    alpha2 = [];
                    p1= theta(:,2);
                    p2 = theta(:,3);
                    alpha3 = [];
                    p3 = [];
                case '2-comp ind'               %(4)
                    alpha1= theta(:,1:obj.K);
                    alpha2 = [];
                    p1 = theta(:,obj.K+1:obj.K*2);
                    p2 = theta(:,obj.K*2+1:end);
                    alpha3 = [];
                    p3 = [];
                case '2-comp constrain-P1'      %(5)
                    alpha1= theta(:,1:obj.K);
                    alpha2 = [];
                    p1 = theta(:,obj.K+1);
                    p2 = theta(:,obj.K+2:end);
                    alpha3 = [];
                    p3 = [];
                case '2-comp constrain-P2'      %(6)
                    alpha1= theta(:,1:obj.K);
                    alpha2 = [];
                    p1 = theta(:,obj.K+1:obj.K*2);
                    p2 = theta(:,obj.K*2+1);
                    alpha3 = [];
                    p3 = [];
                case '2-comp constrain-P1P2'    %(7)
                    alpha1 = theta(:,1:obj.K);
                    alpha2 = [];
                    p1 = theta(:,obj.K+1);
                    p2 = theta(:,obj.K+2);
                    alpha3 = [];
                    p3 = [];
                case '3-comp'                   %(8)
                    alpha1 = theta(:,1);
                    alpha2 = theta(:,2);
                    alpha3 = 1-alpha1-alpha2;
                    p1 = theta(:,3);
                    p2 = theta(:,4);
                    p3 = theta(:,5);
                case '3-comp ind'               %(9)
                    alpha1 = theta(:,1:obj.K);
                    alpha2 = theta(:,obj.K+1:obj.K*2);
                    alpha3 = 1-alpha1-alpha2;
                    p1 = theta(:,obj.K*2+1:obj.K*3);
                    p2 = theta(:,obj.K*3+1:obj.K*4);
                    p3 = theta(:,obj.K*4+1:obj.K*5);
                case '3-comp constrain-P1'      %(10)
                case '3-comp constrain-P2'      %(11)
                case '3-comp constrain-P3'      %(12)
                case '3-comp constrain-P1P2'    %(13)
                case '3-comp constrain-P1P3'    %(14)
                case '3-comp constrain-P2P3'    %(15)
                case '3-comp constrain-P1P2P3'  %(16)
                    alpha1 = theta(:,1:obj.K);
                    alpha2 = theta(:,obj.K+1:obj.K*2);
                    alpha3 = 1-alpha1-alpha2;
                    p1 = theta(:,obj.K*2+1);
                    p2 = theta(:,obj.K*2+2);
                    p3 = theta(:,obj.K*2+3);
            end
        end
        function [out_alpha1,out_alpha2, out_alpha3, out_p1,out_p2,out_p3] = thetasqareTheta(obj,in_alpha1,in_alpha2,in_alpha3,in_p1,in_p2,in_p3)
            datasize = obj.K;
            [out_alpha1] = makesize(datasize,in_alpha1);
            [out_alpha2] = makesize(datasize,in_alpha2);
            [out_alpha3] = makesize(datasize,in_alpha3);
            [out_p1] = makesize(datasize,in_p1);
            [out_p2] = makesize(datasize,in_p2);
            [out_p3] = makesize(datasize,in_p3);
            function [output] = makesize(datsize,input)
                output = input;
                if ~isempty(output)
                    if ~(size(output,2)==datsize)
                        output = repmat(input,1,datsize);
                    end
                end
            end
        end
        %% Extra Functions to get model specific info
        function Nparams = modelNParam(obj,modelName)
            output = getModelInfo(obj,modelName);
            Nparams = output.Nparam;
        end
        function Nalphas = modelNalphas(obj,modelName)
            output = getModelInfo(obj,modelName);
            Nalphas = output.Nalphas;
        end
        function Ncomp = modelNcomps(obj,modelName)
            output = getModelInfo(obj,modelName);
            Ncomp = output.Ncomp;
        end
        function modelParamDesc = modelParamDesc(obj,modelName)
            output = getModelInfo(obj,modelName);
            modelParamDesc = output.modelParamDesc;
        end
        %% Actual Optimization Calls: all P params reported in frames      
        function [theta_mle,mle_llh,std_err,modelInfo] = MLE(obj, constrainedP_bool,fixedP_vals)
            if nargin<3
                fixedP_vals = [];
            end
            if ischar(constrainedP_bool)
                modelName = constrainedP_bool;
            else
            modelName = obj.getModelName(constrainedP_bool);
            end
            modelInfo = obj.getModelInfo(modelName,fixedP_vals);
            modelInfo.problem.objective = @(theta) -obj.fmincon_LLH(modelName,theta);
            switch modelInfo.Ncomp
                case 1
                    % to correct for Offset - create new object and
                    % subtract Offset from all track values
%                     gd_1comp = GeometricDist_MLEandMCMC();
%                     gd_1comp.K = obj.K;
%                     gd_1comp.setData(cellmap(@(x) x-obj.Offset, obj.D));
                    switch modelName
                        case '1-comp'
                            theta_mle = 1 / (obj.Stot/obj.Ntot - obj.Offset);     
                        case '1-comp ind'
                            theta_mle = 1 ./ (obj.S./obj.N - obj.Offset)';      
                    end
                    [mle_llh, mle_llh_conds] = obj.LLH(modelName, theta_mle);   
                    std_err = 1./sqrt(obj.Ntot./(theta_mle.^2) + (obj.Stot-obj.Ntot*obj.Offset)./(1-theta_mle).^2);
                case {2,3}
                    % Offset corrected in LLH call for 2 and 3 component
                    % models
                    [theta_mle, mle_llh, std_err] = obj.fmincon_MLE(modelInfo);
                    theta_mle = theta_mle(:)';
                    mle_llh = mle_llh(:)';
                    std_err = std_err(:)';
            end
        end
        function [mle_llh] = MLE_profile(obj, constrainedP_bool,fixedP_vals,eqVec)
            if nargin<3
                fixedP_vals = [];
            end
            if ischar(constrainedP_bool)
                modelName = constrainedP_bool;
            else
            modelName = obj.getModelName(constrainedP_bool);
            end
            modelInfo = obj.getModelInfo(modelName,fixedP_vals);
            modelInfo.problem.objective = @(theta) -obj.fmincon_LLH(modelName,theta);
            Aeq = zeros(size(eqVec)); Aeq(eqVec>0) = 1; Beq = eqVec(eqVec>0);
            modelInfo.problem.Aeq = Aeq;
            modelInfo.problem.Beq = Beq;
            switch modelInfo.Ncomp
                case 1
                    % to correct for Offset - create new object and
                    % subtract Offset from all track values
%                     gd_1comp = GeometricDist_MLEandMCMC();
%                     gd_1comp.K = obj.K;
%                     gd_1comp.setData(cellmap(@(x) x-obj.Offset, obj.D));
                    switch modelName
                        case '1-comp'
                            theta_mle = 1 / (obj.Stot/obj.Ntot - obj.Offset);     
                        case '1-comp ind'
                            theta_mle = 1 ./ (obj.S./obj.N - obj.Offset)';      
                    end
                    [mle_llh, mle_llh_conds] = obj.LLH(modelName, theta_mle);   
                    std_err = 1./sqrt(obj.Ntot./(theta_mle.^2) + (obj.Stot-obj.Ntot*obj.Offset)./(1-theta_mle).^2);
                case {2,3}
                    % Offset corrected in LLH call for 2 and 3 component
                    % models
                    [~, mle_llh, ~, ~] = obj.fmincon_MLE(modelInfo);          
            end
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC(obj,constrainedP_bool,fixedP_vals)
            if nargin<3
                fixedP_vals = [];
            end
            if ischar(constrainedP_bool)
                modelName = constrainedP_bool;
            else
                modelName = obj.getModelName(constrainedP_bool);
            end
            modelInfo = obj.getModelInfo(modelName,fixedP_vals);
            [theta_mle,~,~] = obj.MLE(constrainedP_bool,fixedP_vals);
            modelInfo.MLE = theta_mle;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_sample(modelInfo) ;
        end
        %% LLH Calculation: 
        function [llh, llh_conds] = LLH(obj, modelName, theta)
            % [in]
            % model: string with name of model. taken from obj.MODELS
            % theta: size Nparams x Ntheta.  Ntheta is number of different
            %   thetas to evaluate (normally Ntheta=1).  Nparams is the
            %   number of parameters in theta for the given model type.
            % [out]
            %  llh - size: 1 x Ntheta.   This is the log-likelihood
            %  llh_conds - size Nconditions x Ntheta.  This is the
            %  log-likelihood for each conditions for each theta value.
            [alpha1,alpha2,alpha3,p1,p2,p3] = getValsFromTheta(obj,theta,modelName);
            [alpha1,alpha2,alpha3,p1,p2,p3] = thetasqareTheta(obj,alpha1,alpha2,alpha3,p1,p2,p3);
            %             p2 = n x k array. k is number of conditions. n is number of p2 values to return info for.
            Ncomp = modelNcomps(obj,modelName); %to speed things up could use another method to determine number of components (that doesn't require getModelInfo call)
            switch Ncomp
                case 1
                    llh_conds = zeros(size(p1'));
                    for k=1:obj.K
                        llh_conds(k,:) = obj.N(k) * log(p1(:,k)) + (obj.S(k)-obj.N(k)*obj.Offset-obj.N(k)) * log(1-p1(:,k));
                    end
                    llh = sum(llh_conds);
                case 2
                    llh_conds = zeros(size(p1'));
                    for k=1:obj.K
                        llh_conds(k,:) = obj.computeLLH_2comp(obj.D{k},alpha1(:,k), p1(:,k),p2(:,k),obj.Offset);
                    end
                    llh = sum(llh_conds,1);
                case 3
                    if isrow(theta)
                        theta=theta(:);
                    end
                    ntheta = size(theta,2);
                    llh_conds = zeros(obj.K, ntheta);
                    for k=1:obj.K
                        for n=1:ntheta
                            llh_conds(k,n) = obj.computeLLH_3comp(obj.D{k}, alpha1(:,k),alpha2(:,k),p1(:,k),p2(:,k),p3(:,k),obj.Offset);
                        end
                    end
                    llh = sum(llh_conds);
            end
        end 
        %% Model Selection
        function [ModelProbabilities,selectedModel,LogEvidence_Array,bestTheta,allsamples] = evidence_Test(obj,models)
            Nsamples = 20000;
            BurnFrac = 0.25;
            LogEvidence_Array = zeros(1,length(models));
            allsamples = cell(1,numel(models));
            for nn=1:length(models)
                [samples, ~, ~, ~] = obj.MCMC(models{nn});
                LogEvidence_Array(nn) = obj.calcLogEvidence(samples,models{nn});
                allsamples{nn} = samples;
            end
            ModelProbabilities=zeros(1,length(models));
            for nn=1:length(models)
                Tmp = LogEvidence_Array;
                LogNN=LogEvidence_Array(nn);
                Tmp(nn)=[];
                ModelProbabilities(nn)=1/(sum(exp(Tmp-LogNN))+1);
            end
            select_id = ModelProbabilities==max(ModelProbabilities);
            selectedModel = models{select_id};
            bestTheta = mean(allsamples{select_id});
            Nalpha = modelNalphas(obj,selectedModel);
            bestTheta(Nalpha+1:end) = obj.convert_prob2offrate(bestTheta(Nalpha+1:end));
            if obj.framerate>0
                bestTheta(Nalpha+1:end) = obj.covert_offrateFrames2Seconds(bestTheta(Nalpha+1:end),obj.framerate);
                display('bestTheta value as offrate in sec^-1, allsamples ps-values in frames');
            else
                display('bestTheta and allsamples ps-values in frames');
            end
        end
        function BIC_Test(obj)
        end
        function AIC_Test(obj)
        end
        function BayesFactor_Test(obj)
        end
        %% Simulate: input dependent on obj.framerate parameter
        function simulate(obj,alphas,ps,Nsim)
            % Use this to simulate data from Geometric Distribution.
            % alphas: k x 1 vector with k the number of simulated datasets
            %         alpha value represents the fraction of data for that component
            %         if k = 1 and obj.K > 1 then use same k param for all
            %         obj.K simulated datasets
            % ps: k x n matrix with k the number of simulated datasets and n the number of components for the geometric model
            %         p value represents probabilities in frames if default
            %         obj.tracklenghth = 0, or if obj.tracklength is set
            %         (~= 0) assumes input is in seconds^-1 and converts to
            %         frames for simulation using obj.tracklength
            %         ps must be given in decreasing value (fastest to slowest)
            %         if k = 1 and obj.K > 1 then use same k param for all
            % Important:
            % obj.K and obj.Offset must be defined: uses internal
            % simulateoutput function to generate data, removes tracks with
            % tracklength < Offset, and runs until Nsim number of
            % datapoints.         
            % - ie: 1 component model with 3 datasets:
            %    alphas = [];
            %    ps = [0.3;0.05;0.1];
            % - ie: 2 component model with 3 datasets:
            %    alphas = [0.4;0.5;0.1];
            %    ps = [0.3,0.05; 0.2,0.1; 0.5,0.03];
            % - ie: 3 component model with 3 datasets:
            %    alphas = [0.4,0.2;0.5,0.1;0.1,0.8];
            %    ps = [0.3,0.1,0.05; 0.2,0.1,0.01; 0.5,0.2,0.03];
            
            % because Offset value means that many simulated datasets won't meet the minimum
            % tracklength requirement, use while statement to ensure Nsim > obj.Offset
            % ------optional return of newD just gives what the actual outputs without the Offset
            if ~isempty(obj.D)
                display('(Warning: simulation results just overwrote previous obj.D values.)');
                obj.clearData();
            end
            if nargin < 4 % no Nsims given use default
                Nsim = 10000 * ones(1,obj.K);
            end
            if isscalar(Nsim)
                Nsim=Nsim * ones(1,obj.K);
            end
            % if obj.framerate is set then input is in sec^-1 so convert
            % offrate into frames for simulation
            if obj.framerate>0
                ps = obj.covert_offrateSeconds2Frames(ps,obj.framerate);
               display(['Assuming rate input given sec^-1 because obj.framerate is set.' ...
               ' Set obj.framerate = 0 to give input in frames.']);
            else
                display(['Assuming rate input given in frames. Set obj.frames to give input in sec^-1.']);
            end
            ps = obj.convert_offrate2prob(ps);
            addon = 0;
            in_Nsim = Nsim;
            while isempty(obj.D) || any(cellfun(@numel,obj.D) < Nsim)
                in_Nsim = in_Nsim + addon;
                newD = simulateoutput(obj,alphas,ps,in_Nsim);
                obj.setData(newD);
                addon = addon + Nsim;
            end
            setD = cellmap(@(x,y) x(randperm(size(x,1),y)),obj.D,num2cell(Nsim)); %ensure that only Nsim number of datapoints
            obj.setData(setD);
            function newD = simulateoutput(obj,alphas,ps,Nsim) % helper function
                newD=cell(1,obj.K);
                if isempty(alphas) % ---1 component model
                    assert(size(ps,2) == 1,'1 component models can only have a P value input (no alpha parameter)');
                    if size(ps,1) == 1 % if only one value input - use same value for all datasets
                        ps = repmat(ps,obj.K,1);
                    else
                        assert(numel(ps) == obj.K, 'if more than one sim value is given then there needs to be one input per simulated dataset');% if more than one sim value is given then there needs to be one input per simulated dataset
                    end
                    assert(all(all(0<ps)) & all(all(ps<1)),' make sure prob values are within range'); % make sure prob values are within range
                    newD = cellmap(@(k) 1+geornd(ps(k),1,Nsim(k)), 1:obj.K);
                elseif size(ps,2)== 2 % ---2 component model
                    assert(size(ps,1) == size(alphas,1),'make sure same # of inputs for ps and alphas') % make sure same # of inputs for ps and alphas
                    if size(alphas,1) == 1 % if only one value input - use same value for all datasets
                        alphas = repmat(alphas,obj.K,1);
                        ps = repmat(ps,obj.K,1);
                    end
                    assert(all(all(0<alphas) & all(alphas<1) & all(0<ps) & all(ps<1)),'make sure alpha and prob values are within range'); % make sure alpha and prob values are within range
                    assert(all(ps(:,1)>ps(:,2)),'make sure probs reported in decreasing value'); % make sure probs reported in decreasing value
                    for k=1:obj.K
                        pop1 = sum(rand(1,Nsim(k))<alphas(k));
                        newD{k} = 1+[geornd(ps(k,1),1,pop1), geornd(ps(k,2),1,Nsim(k)-pop1)];
                    end
                    
                elseif size(ps,2) == 3 % ---3 component model
                    assert(size(ps,1) == size(alphas,1),'make sure same # of inputs for ps and alphas') % make sure same # of inputs for ps and alphas
                    if size(ps,1) == 1 % if only one value input - use same value for all datasets
                        alphas = repmat(alphas,obj.K,1);
                        ps = repmat(ps,obj.K,1);
                    end
                    assert(all(all(0<alphas) & all(alphas<1)) && all(all(0<ps) & all(ps<1)),'Make sure alpha and prob values are within range'); % make sure alpha and prob values are within range
                    assert(all(alphas(:,1)+alphas(:,2) < 1), 'Check alphas do not add to more than 1'); %check alphas don't add to more than 1.
                    assert(all(ps(:,1)>=ps(:,2)) && all(ps(:,2) >= ps(:,3)),'Always should have Ps in decresaing values. (longer components last)'); %always should have ps in decresaing values. (longer components last)
                    for k=1:obj.K
                        samp = rand(1,Nsim(k));
                        pop1 = sum(samp<alphas(k,1));
                        pop2 = sum(alphas(k,1)<=samp & samp<alphas(k,1)+alphas(k,2));
                        pop3 = Nsim(k)-pop1-pop2;
                        newD{k} = 1+[geornd(ps(k,1),1,pop1), geornd(ps(k,2),1,pop2), geornd(ps(k,3),1,pop3)];
                    end
                end
            end
        end
        %% Plotting Helper Methods: units of plotted P values dependent on obj.framerate
        function plotLLH(obj,modelName,fixedP_vals)
            if nargin<3
                fixedP_vals = [];
            end
            if iscell(modelName)
                for n=1:numel(modelName)
                    obj.plotLLH_ndims(modelName{n},fixedP_vals)
                end
            else
                obj.plotLLH_ndims(modelName,fixedP_vals)
            end
        end
        function plotProfileLLH(obj,modelName)
            if nargin<3
                fixedP_vals = [];
            end
            if iscell(modelName)
                for n=1:numel(modelName)
                    obj.plotProfileLLH_ndims(modelName{n},fixedP_vals)
                end
            else
                obj.plotProfileLLH_ndims(modelName,fixedP_vals)
            end
        end
        function plotMCMC(obj,model,samples)
            % Plots histogram of each parameter values from chain
            % For p values, if obj.framerate is set uses sec^-1, otherwise
            % if obj.framerate = 0 reports values in frames
            Nparams = obj.modelNParam(model);
            desc = obj.modelParamDesc(model);
            Nalphas = obj.modelNalphas(model);
            [mean_vals, std_error, outsamples] = obj.calcMCMCstats(samples, model); 
            nrows = 3;
            ncols = ceil(Nparams/nrows);
            Nsamples = size(samples,1);
            f1 = figure;
            f1.Position = [10 10 900 600]; 
%             xmax = max(max(samples(:,Nalphas+1:end)));
%             xmin = min(min(samples(:,Nalphas+1:end)));
            for n = 1:Nparams
                subplot(nrows,ncols,n);
                %check if fixed
                if all(samples(:,n) == repmat(samples(1,n),Nsamples,1));
                    line(repmat(samples(1,n),Nsamples,1),(1:Nsamples));
                    leg_string = sprintf('Mean: %.3g\n',mean_vals(n));
                else
                histogram(outsamples(:,n));
%                 if n>Nalphas
%                 xlim([xmin xmax]);
%                 end
                leg_string = sprintf('Mean: %.3g\n 95CI: [%.3g - %.3g]\n',mean_vals(n), std_error(1,n),std_error(2,n));
                end
                legend(leg_string, 'FontSize',6);
                if n>Nalphas
                    addon = [' - ' obj.rateunits];
                else 
                    addon = '';
                end
                title([strrep(desc{n},'_','-') addon]);
            end
        end
        function plotMCMCchain(obj,model,samples, proposed_samples)
            Nparams = obj.modelNParam(model);
            Nalphas = obj.modelNalphas(model);
            desc = obj.modelParamDesc(model);
            f1 = figure;
            f1.Position = [15 15 900 600]; 
            nrows = 3;
            ncols = ceil(Nparams/nrows);
            for n = 1:Nparams
                subplot(nrows,ncols,n);
                if n<=Nalphas
                props = proposed_samples(:,n);
                samps = samples(:,n);
                addon = 'Fraction';
                else
                    props = obj.convert_prob2offrate(proposed_samples(:,n));
                    samps = obj.convert_prob2offrate(samples(:,n));
                    if obj.framerate>0
                        props = obj.covert_offrateFrames2Seconds(proposed_samples(:,n),obj.framerate);
                        samps = obj.covert_offrateFrames2Seconds(samples(:,n),obj.framerate);
                    end
                    addon = [' - ' obj.rateunits];
                end
                plot(props,'o'); hold on; plot(samps,'*');
                xlabel('Iteration'); legend({'Proposed Values', 'Accepted Values'});
                title([strrep(desc{n},'_','-') addon]);
            end
        end
        function [theta, std_er, samples] = calcMCMCstats(obj, samples, model,a)
            if nargin<4
                a = 0.682; % one standard error. change this to generate different statistics
            end
            Nalphas = obj.modelNalphas(model);
            Nsamples = size(samples,1);
            %---- remove if not used
            Nparams = size(samples,2);
            fixed = all(samples == repmat(samples(1,:),Nsamples,1)); % boolean row vector
            %----
            theta = mean(samples);
            lower = (1-a)/2; upper = 1-((1-a)/2);
            std_er = quantile(samples,[lower, upper],1);
            theta(Nalphas+1:end) = obj.convert_prob2offrate(theta(Nalphas+1:end));
            std_er(:,Nalphas+1:end) = obj.convert_prob2offrate(std_er(:,Nalphas+1:end)); 
            if obj.framerate>0
                theta(Nalphas+1:end) = obj.covert_offrateFrames2Seconds(theta(Nalphas+1:end),obj.framerate);
                std_er(:,Nalphas+1:end) = obj.covert_offrateFrames2Seconds(std_er(:,Nalphas+1:end),obj.framerate);
            end
        end
        function [theta, std_er, theta_llh] = calcMLEstats(obj,theta_mle,mle_llh,std_err,model)             
            Nalphas = obj.modelNalphas(model);
            assert(size(theta_mle,2) == size(std_err,2),'number of fit parameters must equal number of st_err values returned')
            theta = theta_mle;
            std_er = [theta - std_err; theta + std_err];
            theta_llh = mle_llh;
            theta(Nalphas+1:end) = obj.convert_prob2offrate(theta(Nalphas+1:end));
            std_er(:,Nalphas+1:end) = obj.convert_prob2offrate(std_er(:,Nalphas+1:end)); 
            if obj.framerate>0
                theta(Nalphas+1:end) = obj.covert_offrateFrames2Seconds(theta(Nalphas+1:end),obj.framerate);
                std_er(:,Nalphas+1:end) = obj.covert_offrateFrames2Seconds(std_er(:,Nalphas+1:end),obj.framerate);
            end
        end
        function CDF_out = computeCDF_2comp(obj,alpha,p1,p2)
            if numel(alpha == 3)
                p1 = alpha(2);
                p2 = alpha(3);
                alpha = alpha(1);
            end
            data = obj.D;
            offset = obj.Offset;
            CDF_out = cellmap(@(x) CDFhelper(obj,alpha,p1,p2,x,offset),data);
            function CDFhelp = CDFhelper(obj,alpha,p1,p2,dat,offs)
                originalCDF = CDF_2comp(obj,alpha,p1,p2,dat);
                CDF_offset = CDF_2comp(obj,alpha,p1,p2,offs);
                CDFhelp = (originalCDF - CDF_offset)./(1-CDF_offset);
            end
        end
        function CDF_out = computeCDF_3comp(obj,alpha1,alpha2,p1,p2,p3)
            if numel(alpha1 == 5)
                p1 = alpha1(3);
                p2 = alpha1(4);
                p3 = alpha1(5);
                alpha2 = alpha1(2);
                alpha1 = alpha1(1);
            end
            data = obj.D;
            offset = obj.Offset;
            CDF_out = cellmap(@(x) CDFhelper(obj,alpha1,alpha2,p1,p2,p3,x,offset),data);
            function CDFhelp = CDFhelper(obj,alpha1,alpha2,p1,p2,p3,dat,offs)
                originalCDF = CDF_2comp(obj,alpha1,alpha2,p1,p2,p3,dat);
                CDF_offset = CDF_2comp(obj,alpha1,alpha2,p1,p2,p3,offs);
                CDFhelp = (originalCDF - CDF_offset)./(1-CDF_offset);
            end
        end
        %% Analysis User Calls
        function [t, err_t,f] = fit_MLE(obj, constrainedP_bool,fixedP_vals)
            % Fit data using Maximum Likelihood Estimator
            % Results are returned in a table
            modelName = obj.getModelName(constrainedP_bool);
            if nargin<3
                fixedP_vals = [];
            end
            [theta_mle,mle_llh,std_err] = MLE(obj,modelName,fixedP_vals);
            [theta, std_err, ~] = calcMLEstats(obj,theta_mle,mle_llh,std_err,modelName);
            [t,err_t,f] = plotResults(obj,theta,std_err,modelName);
        end
        function [t, err_t,f] = fit_MCMC(obj,constrainedP_bool,fixedP_vals)
            % Fit data using Markov Chain Monte Carlo
            % Results are returned in a table
            modelName = obj.getModelName(constrainedP_bool);
            if nargin<3
                fixedP_vals = [];
            end
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC(obj,constrainedP_bool,fixedP_vals);
            [theta, std_er, samples] = calcMCMCstats(obj, samples, modelName);
            [t,err_t,f] = plotResults(obj,theta,std_er,modelName);
        end
        function [t,err_t,f] = plotResults(obj,theta,std_er,modelName)
        [t, err_t] = makeResultsTable(obj,theta,modelName,std_er);   
        % simulate data based on theta results, using obj info
        [alpha1, alpha2, alpha3, p1, p2, p3] = getValsFromTheta(obj,theta,modelName);
        [out_alpha1,out_alpha2,alpha3,out_p1,out_p2,out_p3] = obj.thetasqareTheta(alpha1,alpha2,alpha3,p1,p2,p3);
        alphas = [out_alpha1', out_alpha2'];
        ps = [out_p1',out_p2',out_p3'];
        Nsims = 20000;
        simobj = GeometricDist_MLEandMCMC();
        simobj.K = obj.K;
        simobj.setOffset(obj.Offset);
        simobj.setframerate(obj.framerate);
        simobj.simulate(alphas,ps,Nsims);
        simobj.setConditionNames(cellmap(@(x) [x ' FitResults'],obj.ConditionNames));
        % generate plot
        f = figure('units','normalized','outerposition',[0 0 0.9 0.9]);
        axH = axes();
        subplot(axH);
        hold('on');    
        data_linespecs = [];
        sim_linespecs.LineStyle = ':';
        sim_linespecs.Color = [0 0 0];
        % make CDF plot of results
        s1 = subplot(2,3,1);
        s1 = make_DataCDFplot(obj,data_linespecs,s1);
        hold on;
        fax = make_DataCDFplot(simobj,sim_linespecs,s1);
        s1.Title.String = sprintf('CDF of Track Lengths and Fit Using Model:  %s', modelName);  
        % make CDF residual plot 
        s2 = subplot(2,3,4);
        fh = obj.makeResidualPlot(fax,s2);
        s2.Title.String = sprintf('Residuals from [Data CDF - Model Fit CDF]');
        s2.YLabel.String = '[Data CDF - Model Fit CDF]';
        s2.XLabel.String = s1.XLabel.String;       
        % make PDF plot of results
        s2 = subplot(2,3,2);hold on;
        [outbins,s2] = make_DataPDFplot(obj,data_linespecs,s2);  
        [outbins,s2] = make_DataPDFplot(simobj,sim_linespecs,s2,outbins);
        s2.Title.String = sprintf('PDF of Track Lengths and Fit Using Model:  %s', modelName);  
        % make table of fit values
        s4 = subplot(2,3,3,'Parent',f);
        s4.Title.String = 'Fit Parameter Values';
        axis('off');
        pos = s4.Position;
        ut = obj.makeUITablefromTable(t,f);
        ut.Units = 'normalized';
        ut.Position = pos;
        ut.Position(3) = ut.Extent(3); 
%                   ut.Position(4) = ut.Extent(4); 
        % make table of fit value error bounds
        s5 = subplot(2,3,[5,6],'Parent',f);
        s5.Title.String = 'Parameter Confidence Interval:  [-\sigma, +\sigma]';
        axis('off');
        pos = s5.Position;
        ut = obj.makeUITablefromTable(err_t,f);
        ut.Units = 'normalized';
        ut.Position = pos;
        ut.Position(3) = ut.Extent(3);
%                   ut.Position(4) = ut.Extent(4);  
        end
        function [t, err_t] = makeResultsTable(obj,theta,modelName,std_err)
            [alpha1,alpha2,alpha3,p1,p2,p3] = getValsFromTheta(obj,theta,modelName);
            [alpha1,alpha2,alpha3,p1,p2,p3] = thetasqareTheta(obj,alpha1,alpha2,alpha3,p1,p2,p3);  
            [ste_alpha1,ste_alpha2,ste_alpha3,ste_p1,ste_p2,ste_p3] = getValsFromTheta(obj,std_err,modelName);
            [ste_alpha1,ste_alpha2,ste_alpha3,ste_p1,ste_p2,ste_p3] = thetasqareTheta(obj,ste_alpha1,ste_alpha2,ste_alpha3,ste_p1,ste_p2,ste_p3); 
            Ncomp = obj.modelNcomps(modelName);
            if ~isempty(obj.ConditionNames)
            DataConditions = obj.ConditionNames;
            else
               DataConditions = cellmap(@(k) sprintf('Condition #%i',k),1:obj.K);
            end
            switch Ncomp
                case 1
%                      t = table(p1,'RowNames',DataConditions,...
%                          'VariableNames',{sprintf('P\n%s', obj.rateunits)});
                     t = table(obj.N, p1','RowNames',DataConditions,...
                         'VariableNames',{'NumTracks','P1'});
                     t.Properties.VariableUnits = {'# Trajectories',obj.rateunits};
                     t.Properties.VariableDescriptions = {'# Tracks',...
                         ['P|' obj.rateunits_html]};
                      err_t = table();
                     err_t.N = obj.N;
                     err_t.PLow = ste_p1(1,:)';
                     err_t.PHigh = ste_p1(2,:)';
                     err_t.Properties.RowNames = DataConditions;
                     err_t.Properties.VariableUnits = {'# Trajectories',obj.rateunits,obj.rateunits};
                     err_t.Properties.VariableDescriptions = {'# Tracks',...
                         ['P|' obj.rateunits_html '|Lower Bound'],['P|' obj.rateunits_html '|Upper Bound']};
                case 2
                     t = table(obj.N,alpha1',p1',p2','RowNames',DataConditions,...
                         'VariableNames',{'NumTracks','Alpha','P1','P2'});
                     t.Properties.VariableUnits = {'# Trajectories','Fraction P1',obj.rateunits,obj.rateunits};
                     t.Properties.VariableDescriptions = {'# Tracks','<html>&alpha<sub>1</sub></htlml>|Fraction P1',...
                         ['P1|' obj.rateunits_html],['P2|' obj.rateunits_html]};
                     err_t = table();
                     err_t.N = obj.N;
                     err_t.ALow = ste_alpha1(1,:)';
                     err_t.AHigh = ste_alpha1(2,:)';
                     err_t.P1Low = ste_p1(1,:)';
                     err_t.P1High = ste_p1(2,:)';
                     err_t.P2Low = ste_p2(1,:)';
                     err_t.P2High = ste_p2(2,:)';
                     err_t.Properties.RowNames = DataConditions;
                     err_t.Properties.VariableUnits = {'# Trajectories','Fraction P1','Fraction P1',obj.rateunits,obj.rateunits,obj.rateunits,obj.rateunits};
                     err_t.Properties.VariableDescriptions = {'# Tracks','<html>&alpha<sub>1</sub></htlml>|Fraction P1|Lower Bound',...
                         '<html>&alpha<sub>1</sub></htlml>|Fraction P1|Upper Bound',...
                         ['P1|' obj.rateunits_html '|Lower Bound'],['P1|' obj.rateunits_html '|Upper Bound'],...
                         ['P2|' obj.rateunits_html '|Lower Bound'],['P2|' obj.rateunits_html '|Upper Bound']};
                case 3
                    t = table(obj.N,alpha1',alpha2',alpha3',p1',p2',p3','RowNames',DataConditions,...
                         'VariableNames',{'NumTracks','Alpha1','Alpha2','Alpha3','P1','P2','P3'});
                     t.Properties.VariableUnits = {'# Trajectories','Fraction P1','Fraction P2','Fraction P3',obj.rateunits,obj.rateunits,obj.rateunits};
                     t.Properties.VariableDescriptions = {'# Tracks','<html>&alpha<sub>1</sub></htlml>|Fraction P1','<html>&alpha<sub>2</sub></htlml>|Fraction P2',...
                         '<html>&alpha<sub>3</sub></htlml>|Fraction P3',...
                         ['P1|' obj.rateunits_html],['P2|' obj.rateunits_html],['P3|' obj.rateunits_html]};
                     err_t = table();
                     err_t.N = obj.N;
                     err_t.A1Low = ste_alpha1(1,:)';
                     err_t.A1High = ste_alpha1(2,:)';
                     err_t.A2Low = ste_alpha2(1,:)';
                     err_t.A2High = ste_alpha2(2,:)';
                     err_t.A3Low = ste_alpha3(1,:)';
                     err_t.A3High = ste_alpha3(2,:)';
                     err_t.P1Low = ste_p1(1,:)';
                     err_t.P1High = ste_p1(2,:)';
                     err_t.P2Low = ste_p2(1,:)';
                     err_t.P2High = ste_p2(2,:)';
                     err_t.P3Low = ste_p3(1,:)';
                     err_t.P3High = ste_p3(2,:)';
                     err_t.Properties.RowNames = DataConditions;
                     err_t.Properties.VariableUnits = {'# Trajectories','Fraction P1','Fraction P1','Fraction P2','Fraction P2',...
                         'Fraction P3','Fraction P3',...
                         obj.rateunits,obj.rateunits,obj.rateunits,obj.rateunits,obj.rateunits,obj.rateunits};
                     err_t.Properties.VariableDescriptions = {'# Tracks','<html>&alpha<sub>1</sub></htlml>|Fraction P1|Lower Bound',...
                         '<html>&alpha<sub>1</sub></htlml>|Fraction P1|Upper Bound',...
                         '<html>&alpha<sub>2</sub></htlml>|Fraction P2|Lower Bound',...
                         '<html>&alpha<sub>2</sub></htlml>|Fraction P2|Upper Bound',...
                          '<html>&alpha<sub>3</sub></htlml>|Fraction P3|Lower Bound',...
                         '<html>&alpha<sub>3</sub></htlml>|Fraction P3|Upper Bound',...
                         ['P1|' obj.rateunits_html '|Lower Bound'],['P1|' obj.rateunits_html '|Upper Bound'],...
                         ['P2|' obj.rateunits_html '|Lower Bound'],['P2|' obj.rateunits_html '|Upper Bound'],...
                         ['P3|' obj.rateunits_html '|Lower Bound'],['P3|' obj.rateunits_html '|Upper Bound']};
            end
            t.Properties.Description = ['Fit Results using Model: ' modelName];
            t.Properties.UserData = Ncomp;
            err_t.Properties.Description = ['Standard Error Bounds from Results using Model: ' modelName];
            err_t.Properties.UserData = Ncomp;
        end
        function in_fax = make_DataCDFplot(obj,linespecs,in_fax)
            % This method plots the CDF of the obj.D data
            % linespecs - a structure with fields containing any line information that can be
            % added to the CDF line objects. Otherwise a set of default values are used.
            % in_fax - input figure ax if plots are to be added to a specified figure axes.
            %        - if in_fax was given as input then output is the
            %          same. Otherwise it is the figure axes of the figure
            %          created within the method.
            % It uses the obj.rateunits to determine the x values,
            % obj.COLORS for the line colors and the obj.ConditionNames for the legend
            if nargin<3
                % no input figure handle
                fh = figure();
                in_fax = axes();
            end
            if obj.K>9
                cols = cell(1,obj.K);
                cm = lines(obj.K);
                for cc = 1:size(cm,1)
                    cols{cc} = cm(cc,:);
                end
            else
                cols = obj.COLORS;
            end
            %set default params
            lineinfo.LineWidth = 1.5;
            lineinfo.LineStyle = '-';
            lineinfo.Color = cols;
            if ~isempty(obj.ConditionNames)
                lineinfo.DisplayName = cellmap(@(x,y) strrep(x,'_','-'), obj.ConditionNames);
            end
            %adjust to lineinfo to linespecs input if any
            if nargin>1 && ~isempty(linespecs)
                assert(isstruct(linespecs),'linespecs must be a structure');
                fname = fieldnames(linespecs);
                for fieldID = 1:numel(fname)
                    lineinfo.(fname{fieldID}) = linespecs.(fname{fieldID});
                end
            end
            %make cdf plot
            data2plot = obj.convertDataUnits(obj.D);
            fnames = fieldnames(lineinfo);
            for cc = 1:obj.K
                [x,f] = ecdf(data2plot{cc});
                currline = plot(in_fax,f(2:end),x(2:end));
                hold('on');
                for fieldID = 1:numel(fnames)
                    lineproperty = lineinfo.(fnames{fieldID});
                    if iscell(lineproperty)
                        assert(numel(lineproperty) >= obj.K, 'Need to input a lineproperty for each obj.K');
                        currline.(fnames{fieldID}) = lineproperty{cc};
                    else
                        currline.(fnames{fieldID}) = lineproperty;
                    end
                end
            end
            ylabel('Cumulative Probability');
            xlabel(obj.rateunits);
            legend('location','best');
            temp_xlim = in_fax.XLim;
            in_fax.XScale='log';
            in_fax.XLim = temp_xlim;
            in_fax.YScale='linear';
        end
        function [in_bins,in_fax] = make_DataPDFplot(obj,linespecs,in_fax,in_bins)
            % This method plots the PDF of the obj.D data (if you want just
            % a histogram see matlab Histogram Method - or change the
            % histogram input within this method call.
            % linespecs - a structure with fields containing any line information that can be
            % added to the PDF line objects. Otherwise a set of default values are used.
            % in_fax - input figure ax if plots are to be added to a specified figure axes.
            %        - if in_fax was given as input then output is the
            %          same. Otherwise it is the figure axes of the figure
            %          created within the method.
            % in_bins - iput for histogram binning, if you want to manually
            %           set this. Useful when plotting into already
            %           existing plot.
            % It uses the obj.rateunits to determine the x values,
            % obj.COLORS for the line colors and the obj.ConditionNames for the legend
            if nargin<3
                % no input figure handle
                figure();
                in_fax = axes();
                in_bins = [];
            end
            if nargin<4
                in_bins = [];
            end
            if obj.K>9
                cols = cell(1,obj.K);
                cm = lines(obj.K);
                for cc = 1:size(cm,1)
                    cols{cc} = cm(cc,:);
                end
            else
                cols = obj.COLORS;
            end
            %set default params
            lineinfo.LineWidth = 1.5;
            lineinfo.LineStyle = '-';
            lineinfo.Color = cols;
            if ~isempty(obj.ConditionNames)
                lineinfo.DisplayName = cellmap(@(x,y) strrep(x,'_','-'), obj.ConditionNames);
            end
            %adjust to lineinfo to linespecs input if any
            if nargin>1 && ~isempty(linespecs)
                assert(isstruct(linespecs),'linespecs must be a structure');
                fname = fieldnames(linespecs);
                for fieldID = 1:numel(fname)
                    lineinfo.(fname{fieldID}) = linespecs.(fname{fieldID});
                end
            end
            %make cdf plot
            data2plot = obj.convertDataUnits(obj.D);
            fnames = fieldnames(lineinfo);
            for cc = 1:obj.K
                if isempty(in_bins)
                [N,in_bins] = histcounts(data2plot{cc},'Normalization','pdf');
                else
                    [N,in_bins] = histcounts(data2plot{cc},in_bins,'Normalization','pdf');
                end
                N(N == 0) = NaN;
                in_bins(N == 0) = NaN;
               
                binwidth = diff(in_bins(1:2));
                currline = plot(in_fax,in_bins(2:end)-binwidth,N);
                hold('on');
                for fieldID = 1:numel(fnames)
                    lineproperty = lineinfo.(fnames{fieldID});
                    if iscell(lineproperty)
                        assert(numel(lineproperty) >= obj.K, 'Need to input a lineproperty for each obj.K');
                        currline.(fnames{fieldID}) = lineproperty{cc};
                    else
                        currline.(fnames{fieldID}) = lineproperty;
                    end
                end
            end
            ylabel('Probility Density');
            xlabel(obj.rateunits);
            legend('location','best');
            temp_xlim = in_fax.XLim;
            in_fax.XScale='log';
            in_fax.XLim = temp_xlim;
            in_fax.YScale='linear';
        end
    end
    methods (Access = protected)
        % for MLE calculation
        function [llh, llh_conds] = fmincon_LLH(obj, modelName, theta_input)
            % this is the LLH function that is the objective function used for
            % fmincon in fmincon_MLE. 
            % theta_input must be a row vector and then this function
            % converts to column vector needed for fmincon optimization
            [llh, llh_conds] = LLH(obj, modelName, theta_input');
        end
        function [theta_mle, mle_llh, std_err] = fmincon_MLE(obj,modelInfo)
            max_iter=500;
            nParam = modelInfo.Nparam;
            sequence=zeros(nParam,max_iter+1);
            llh = zeros(obj.K,max_iter+1);
            sequence(:,1)=modelInfo.problem.x0(:);
            llh(1)=obj.fmincon_LLH(modelInfo.modelName,modelInfo.problem.x0);
            function stop = output(theta, opt, state)
                sequence(:,opt.iteration+1)=theta;
                llh(opt.iteration+1)=-opt.fval;
                stop=strcmp(state,'done');
            end
%             'interior-point' 'sqp'
            opts = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
            if obj.plotDiagnostics
                opts.Diagnostics='on';
                opts.Display='iter-detailed';
            else
                opts.Diagnostics='off';
                opts.Display='off';
            end
            opts.FunValCheck='on';
            opts.MaxIter=max_iter;
            opts.AlwaysHonorConstraints='bounds';
            opts.OutputFcn=@output;
            modelInfo.problem.solver = 'fmincon';
            modelInfo.problem.options = opts;
            if obj.Verbose
            fprintf('Optimizing Model: %s\n', modelInfo.modelName);
            end
            [theta_mle, mle_llh, exitFlag, out, ~, grad, hess] = fmincon(modelInfo.problem);
            %optimize with fmin con
            mle_llh = -mle_llh;
            if obj.Verbose
            fprintf('---> MLE_LLH: %g\n',mle_llh);
            end
            std_err = sqrt(diag(pinv(hess)))'; % this is observed information (hess is -hessian)
            niter = out.iterations;
            desc=modelInfo.modelParamDesc;
            if obj.plotDiagnostics
                figure();
                subplot(2,1,1);
                hold on;
                xs= 1:niter+1;
                for n=1:nParam
                    c=obj.COLORS{mod(n,numel(obj.COLORS))+1};
                    plot(xs,sequence(n,xs),'Color',c,'DisplayName',desc{n});
                end
                xlabel('Iteration');
                ylabel('Parameter Value');
                H=legend('Location','Best');
                H.Interpreter='None';
                title(sprintf('MLE Optimization Report: %s', modelInfo.modelName),'Interpreter','none')
                subplot(2,1,2);
                plot(xs,llh(xs),'Color','k','DisplayName','LLH');
                xlabel('Iteration');
                ylabel('LLH');
                H=legend('Location','Best');
                H.Interpreter='None';
            end
        end
        % for MCMC
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,modelInfo)
            Nparams = modelInfo.Nparam;
            Nalphas = modelInfo.Nalphas;
            Ncomp = modelInfo.Ncomp;
            Nsamples = obj.MCMCSampleNum;
            Ntot = ceil(Nsamples*(1+obj.MCMCBurnFrac)); %#ok<*PROP>
            samples = zeros(Ntot, Nparams);
            samples(1,:) = modelInfo.MLE;
            samplesLLH = zeros(Ntot,1);
            samplesLLH(1) = obj.LLH(modelInfo.modelName,samples(1,:));
            alpha_sigma = obj.MCMCSigma;
            [~, ~, ~, p1, p2, p3] = obj.getValsFromTheta(modelInfo.MLE,modelInfo.modelName);
            lambda_sigma = obj.MCMCSigma*[p1,p2,p3];
            proposed_samples = zeros(Ntot, Nparams);
            proposed_samples(1,:) = samples(1,:);
            proposed_samplesLLH = zeros(Ntot, 1);
            proposed_samplesLLH(1) = samplesLLH(1);
            Naccept = 0;
            Nreject = 0;
            Ninvalid = 0;
            % are some of the params fixed
            for k = 2:Ntot
                propTheta = samples(k-1,:);
                % modify propTheta: modify one variable at a time
                selected_param = mod(k-2,Nparams)+1;
                if ~isempty(modelInfo.fixedP_vals) && modelInfo.fixedP_vals(selected_param)
                    invalid = true; % this is to not make changes if that parameter is hard fixed. MCMCparams.fixed comes in as boolean.
                    proposed_samples(k,:) = propTheta;
                else % normal case where nothing is hard fixed
                    if selected_param<=Nalphas % is this an alpha
                        propTheta(selected_param) = propTheta(selected_param) + randn(1)*alpha_sigma;
                    else
                        propTheta(selected_param) = propTheta(selected_param) + randn(1)*lambda_sigma(selected_param-Nalphas);
                    end
                    proposed_samples(k,:) = propTheta;
                    
                    % check if new values are valid
                    [alpha1, alpha2, alpha3, lambda1, lambda2, lambda3] = obj.getValsFromTheta(propTheta,modelInfo.modelName);
                    if Ncomp == 1  % determine if lambdas are in correct order
                        lambda_invalid = false;
                    elseif Ncomp == 2
                        lambda_invalid = any(lambda1<lambda2);
                    elseif Ncomp == 3 %% hmmmmm is something wrong here/!!!???
                        lambda_invalid = any(lambda1<lambda2) || any(lambda2<lambda3);
                    end 
                    invalid =  lambda_invalid || any(propTheta<0) || any(propTheta(1:Nalphas)>1) ; % alpha must be betwee 0 and 1
                    % if 3 components - make sure alphas sum to less than or equal to 1
                    if Ncomp == 3
                        alpha_invalid = any(alpha1+alpha2 > 1);
                        invalid = invalid || alpha_invalid;
                    end
                end
                % now if not valid continue, otherwise if valid then calculate the LLH
                if invalid
                    proposed_samplesLLH(k) = -inf;
                    Ninvalid = Ninvalid+1;
                else
                    propThetaLLH = obj.LLH(modelInfo.modelName,propTheta);
                    proposed_samplesLLH(k) = propThetaLLH;
                end
                
                % given that it is valid, decided to accept or reject
                if ~invalid && (propThetaLLH>samplesLLH(k-1) || log(rand(1)) < propThetaLLH - samplesLLH(k-1)) % accept condition
                    samples(k,:)= propTheta;
                    if selected_param<=Nalphas
                        assert(samples(k,selected_param)<1,'P parameter must be <1');
                    end
                    samplesLLH(k) = propThetaLLH;
                    Naccept = Naccept+1;
                else % reject
                    samples(k,:) = samples(k-1,:);
                    samplesLLH(k) = samplesLLH(k-1);
                    Nreject = Nreject+1;
                end
            end
            
            % remove burn ins
            samples = samples(Ntot-Nsamples+1 : Ntot,:);
            samplesLLH = samplesLLH(Ntot-Nsamples+1 : Ntot);
            proposed_samples = proposed_samples(Ntot-Nsamples+1 : Ntot,:);
            proposed_samplesLLH = proposed_samplesLLH(Ntot-Nsamples+1 : Ntot);
            fprintf('\nAcceptance Ratio: %g\n#Accepted: %i\n#Rejected: %i\n', Naccept/Ntot,Naccept,Nreject);
            fprintf('\n# Invalids: %g\n', Ninvalid);
        end
        % for plotting functions
        function plotLLH_ndims(obj,modelName,fixedP_vals)
            M=1e3; %number of samples
            spacings = logspace(-3,0,M)';
            xs = [flip(-spacings);0;spacings];
            [mle, mle_llh] =  obj.MLE(modelName,fixedP_vals);
            nParam = obj.modelNParam(modelName);
            desc = obj.modelParamDesc(modelName);
            Ncomp = obj.modelNcomps(modelName);
            Nalphas = obj.modelNalphas(modelName);
            if ~Nalphas;
                figure; ax_p = gca; hold on;
            else
                figure;
                ax_alpha = subplot(2,1,1); hold on;
                ax_p = subplot(2,1,2); hold on;
            end
            for n=1:nParam
                thetas = repmat(mle,size(xs,1),1);
                thetas(:,n)= thetas(:,n) +  mle(:,n).*xs;
                [alpha1, alpha2, alpha3, p1, p2, p3] = getValsFromTheta(obj,thetas,modelName);
                plotxs = thetas(:,n);
                % do some checks to make sure these values are allowed 
                if Ncomp == 1
                    remove_bool = any(0>=p1 | p1>=1,2);
                    thetas(remove_bool,:) = [];
                    plotxs(remove_bool,:) = [];
                end
                if Ncomp == 2
                    remove_bool = any(p1<p2 | 0>=alpha1 | alpha1>=1 | 0>=p1 | p1>=1 | 0>=p2 | p2>=1,2);
                    thetas(remove_bool,:) = [];
                    plotxs(remove_bool,:) = [];
                end
                if Ncomp == 3 %% hmmmmm is something wrong here/!!!???
                        remove_bool = any(p1<p2 | p2<p3,2) | any(alpha1+alpha2 > 1,2);
                        thetas(remove_bool,:) = [];
                        plotxs(remove_bool,:) = [];
                end 
                [llh, ~] = obj.LLH(modelName, thetas);
                c = obj.COLORS{mod(n,numel(obj.COLORS))+1};
                name = sprintf('Theta(%i):%s',n,desc{n});
                if Nalphas && n<=Nalphas
                    plot(ax_alpha,plotxs', llh,'Color',c,'DisplayName',name); hold on;
                    name=sprintf('MLE(%i):%s =%.3g',n,desc{n},mle(n));
                    plot(ax_alpha,mle(n), mle_llh,'Color',c,'Marker','*','DisplayName',name);
                else
                    plotPs = obj.convert_prob2offrate(plotxs');
                    plot_mle = obj.convert_prob2offrate(mle(n));
                    if obj.framerate>0
                        plotPs = obj.covert_offrateFrames2Seconds(plotPs,obj.framerate);
                        plot_mle = obj.covert_offrateFrames2Seconds(plot_mle,obj.framerate);
                    end
                    plot(ax_p,plotPs, llh,'Color',c,'DisplayName',name);
                    name=sprintf('MLE(%i):%s =%.3g',n,desc{n},plot_mle);
                    plot(ax_p,plot_mle, mle_llh,'Color',c,'Marker','*','DisplayName',name);
                end
            end
            if Nalphas
                set(ax_alpha,'XScale','log','YScale','log');
                ylabel(ax_alpha,'LLH Alphas')
                xlabel(ax_alpha,'Alpha Parameter Value');
                title(ax_alpha,sprintf('alpha MLE Estimate: %s',modelName),'Interpreter','none');
                Ha=legend(ax_alpha,'Location','best');
                Ha.Interpreter='None';
            end
            set(ax_p,'XScale','log','YScale','log');
            xlabel(ax_p,[' P - ' obj.rateunits]);
            ylabel(ax_p,'LLH Ps') 
            title(ax_p,sprintf('P- MLE Estimate: %s',modelName),'Interpreter','none');
            Hp=legend(ax_p,'Location','best');
            Hp.Interpreter='None';
        end
        function plotProfileLLH_ndims(obj,modelName,fixedP_vals)
            M=1e1; %number of samples
            spacings = logspace(-3,0,M)';
            xs = [flip(-spacings);0;spacings];
            [mle, mle_llh] =  obj.MLE(modelName,fixedP_vals);
            nParam = obj.modelNParam(modelName);
            desc = obj.modelParamDesc(modelName);
            Ncomp = obj.modelNcomps(modelName);
            Nalphas = obj.modelNalphas(modelName);
            h = waitbar(0);
            if ~Nalphas;
                figure; ax_p = gca; hold on;
            else
                figure;
                ax_alpha = subplot(2,1,1); hold on;
                ax_p = subplot(2,1,2); hold on;
            end
            thetas = repmat(mle,size(xs,1),1);
            thetas = thetas + xs * mle;
            [alpha1, alpha2, alpha3, p1, p2, p3] = getValsFromTheta(obj,thetas,modelName);
            if Ncomp == 1
                remove_bool = any(0>=p1 | p1>=1,2);
                thetas(remove_bool,:) = [];
            end
            if Ncomp == 2
                remove_bool = any(p1<p2 | 0>=alpha1 | alpha1>=1 | 0>=p1 | p1>=1 | 0>=p2 | p2>=1,2);
                thetas(remove_bool,:) = [];
            end
            if Ncomp == 3 %% hmmmmm is something wrong here/!!!???
                remove_bool = any(p1<p2 | p2<p3,2) | any(alpha1+alpha2 > 1,2);
                thetas(remove_bool,:) = [];
            end
            for n=1:nParam
                profile_mle = zeros(size(thetas,1),1);
                for mm = 1:size(thetas,1)
                    InEqVec = zeros(1,size(thetas,2));
                    InEqVec(n) = thetas(mm,n);
                    profile_mle(mm) =  obj.MLE_profile(modelName,[],InEqVec);
                end
                % do some checks to make sure these values are allowed
                c = obj.COLORS{mod(n,numel(obj.COLORS))+1};
                name = sprintf('Theta(%i):%s',n,desc{n});
                if Nalphas && n<=Nalphas
                    plot(ax_alpha,thetas(:,n), profile_mle,'Color',c,'DisplayName',name); hold on;
                    name=sprintf('MLE(%i):%s =%.3g',n,desc{n},mle(n));
                    plot(ax_alpha,mle(n), mle_llh,'Color',c,'Marker','*','DisplayName',name);
                else
                    plotPs = obj.convert_prob2offrate(thetas(:,n));
                    plot_mle = obj.convert_prob2offrate(mle(n),obj.framerate);
                    if obj.framerate>0
                        plotPs = obj.covert_offrateFrames2Seconds(plotPs,obj.framerate);
                        plot_mle = obj.covert_offrateFrames2Seconds(plot_mle,obj.framerate);
                    end
                    plot(ax_p,plotPs, profile_mle,'Color',c,'DisplayName',name);
                    name=sprintf('MLE(%i):%s =%.3g',n,desc{n},plot_mle);
                    plot(ax_p,plot_mle, mle_llh,'Color',c,'Marker','*','DisplayName',name);
                end
                if Nalphas
                    set(ax_alpha,'XScale','log','YScale','log');
                    ylabel(ax_alpha,'MLE LLH Alphas')
                    xlabel(ax_alpha,'Alpha Parameter Value');
                    title(ax_alpha,sprintf('Alpha Profile Likelihood: %s',modelName),'Interpreter','none');
                    Ha=legend(ax_alpha,'Location','best');
                    Ha.Interpreter='None';
                end
                waitbar(n/nParam,h);
            end
            set(ax_p,'XScale','log','YScale','log');
            xlabel(ax_p,[' P - ' obj.rateunits]);
            ylabel(ax_p,'MLE LLH Ps')
            title(ax_p,sprintf('P- Profile Likelihood: %s',modelName),'Interpreter','none');
            Hp=legend(ax_p,'Location','best');
            Hp.Interpreter='None';
            close(h);
        end
        % for evidence calculation
        function LE=calcLogEvidence(obj,samples,model)
            %Calculate the evidence using Laplace approximation - from Keith
            % Tierney Kadane 1986
            % larger is better
            theta=mean(samples);
            goodids = logical(max(samples,[],1) - min(samples,[],1));
            A=inv(cov(samples(:,goodids)));
            LP=log(obj.calcPrior(theta,model));
            nn=size(samples,2);
            [LL, ~] = obj.LLH(model, theta);     
            LE=LL+ LP + nn/2*log(2*pi)-1/2*log(det(A));
        end  
        function p = calcPrior(obj,theta,model)
            na = obj.modelNalphas(model); %alpha prior =1 because uniform [0 1] so just 1 - don't include
            p = prod(1 ./ sqrt( theta(na+1:end) .* (1-theta(na+1:end)) ) ); %prior on all lambda Jeffrey's Prior
        end
        function data = convertDataUnits(obj,data_frames)
            % this method takes in data in frames (simulated or real) and
            % converts it into the units set by the class - seconds if
            % obj.framerate~=0 or frames if obj.framerate=0
            % this method is a helper method for plotting the results.
            % data_frames is 1 x obj.K cell array of tracklengths in frames
            if obj.framerate && obj.framerate ~= 0
                data = cellmap(@(x) obj.convert_frames2seconds(x,obj.framerate), data_frames);
            else
                data = data_frames;
            end
        end
        function CDF_out = CDF_2comp(obj,alpha,p1,p2,input)
            assert(nargin>3,'Inputs must include alpha,p1,p2,input');
            CDF_out = (alpha.*(1-(1-p1).^(input)) + (1-alpha).*(1-(1-p2).^(input)));
        end
        function CDF_out = CDF_3comp(obj,alpha1,alpha2,p1,p2,input)
            assert(nargin>3,'Inputs must include alpha1,alpha2,p1,p2,p3,input');
            alpha3 = 1-alpha1-alpha2;
            CDF_out = (alpha1.*(1-(1-p1).^(input)) + alpha2.*(1-(1-p2).^(input)) + alpha3.*(1-(1-p3).^(input)));
        end
    end
    methods (Static=true)
        %% LLH Calculations - static
        function modelName = getModelName(constrainedP_bool)
            % First find out what model is from the input ----------------------
            possible_modelstrNames = GeometricDist_MLEandMCMC.MODELS;
            if nargin == 0 % input was empty so just use default model
                display('No input model given. Defaulting to 1-component Geometric Distribution.');
                modelName = '1-comp ind';
            elseif ischar(constrainedP_bool) % input was a string with the model name
                matchname = cellfun(@(x) strcmp(x,constrainedP_bool), possible_modelstrNames);
                if ~any(matchname)
                    display('Input model string must be one of the following: 1-comp, 2-comp, 3-comp');
                end
                assert(any(matchname),'Input model string must be one of the following: 1-comp, 2-comp, 3-comp');
                modelName = constrainedP_bool;
            else % if input is not a string it should be vector, start by checking size of input vector
                assert(isrow(constrainedP_bool),'make sure constrainedP_bool is an input vector'); %make sure it is an input vector
                Np = numel(constrainedP_bool);
                if Np == 1
                    if constrainedP_bool
                        modelName = '1-comp';
                        if ~(constrainedP_bool == 1)
                            display('Warning: input expected to be boolean vector. Non-boolean values ignored. To fix values add a second input vector.')
                        end
                    else
                        modelName = '1-comp ind';
                    end
                elseif Np == 2
                    if constrainedP_bool(1) == 1
                        if constrainedP_bool(2) == 1
                            modelName = '2-comp constrain-P1P2';
                        else
                            modelName = '2-comp constrain-P1';
                        end
                    else
                        if constrainedP_bool(2) == 1
                            modelName = '2-comp constrain-P2';
                        else
                            modelName = '2-comp ind';
                        end
                    end
                elseif Np == 3
                    if constrainedP_bool(1)
                        if constrainedP_bool(2)
                            if constrainedP_bool(3)
                                modelName = '3-comp constrain-P1P2P3';
                            else
                                modelName = '3-comp constrain-P1P2';
                            end
                        else
                            if constrainedP_bool(3)
                                modelName = '3-comp constrain-P1P3';
                            else
                                modelName = '3-comp constrain-P1';
                            end
                        end
                    elseif constrainedP_bool(2)
                        if constrainedP_bool(3) == 1
                            modelName = '3-comp constrain-P2P3';
                        else
                            modelName = '3-comp constrain-P2';
                        end
                    elseif constrainedP_bool(3)
                        modelName = '3-comp constrain-P3';
                        
                    else
                        modelName = '3-comp ind'
                    end
                else
                    display('Invalid input: model input must be a boolean vector size (1 x n), with n value {1,2,3}');
                    return;
                end
            end
        end
        function constrainedP_bool = getConstrainedP_bool(modelName)
            switch modelName
                case '1-comp'
                    constrainedP_bool = 1;
                case '1-comp ind'
                    constrainedP_bool = 0;
                case '2-comp'
                    constrainedP_bool = [];
                    display(['There is no corresponding constrainedP_bool for this model.'...
                        'Use ''2-comp'' string to represent this on input.']);
                case '2-comp ind'
                    constrainedP_bool = [0 0];
                case '2-comp constrain-P1'
                    constrainedP_bool = [1 0];
                case '2-comp constrain-P2'
                    constrainedP_bool = [0 1];
                case '2-comp constrain-P1P2'
                    constrainedP_bool = [1 1];
                case '3-comp'
                    constrainedP_bool = [1 1 1];
                case '3-comp ind'
                    constrainedP_bool = [0 0 0];
                case '3-comp constrain-P1'
                    constrainedP_bool = [1 0 0];
                case '3-comp constrain-P2'
                    constrainedP_bool = [0 1 0];
                case '3-comp constrain-P3'
                    constrainedP_bool = [0 0 1];
                case '3-comp constrain-P1P2'
                    constrainedP_bool = [1 1 0];
                case '3-comp constrain-P1P3'
                    constrainedP_bool = [1 0 1];
                case '3-comp constrain-P2P3'
                    constrainedP_bool = [0 1 1];
                case '3-comp constrain-P1P2P3'
                    constrainedP_bool = [1 1 1];
            end
        end
        function llh = computeLLH_2comp(data,alpha,p1,p2,offset)
            % theta = [alpha, p1 (fast rate), p2 (slow rate)]
            % data: column vector of integers>=1
            % alpha, p1, p2: can be scaler or column vector (to test many values at once)
            %                   this is mainly for plotLLH visualization
            %                between 0 and 1
%             assert(all(p1>=p2));
%             assert(all(all(0<[alpha,p1,p2]) & all([alpha,p1,p2]<1)));
            assert(iscolumn(data), 'input data must be in column format')
            if size(alpha,1)>1 
                assert(isequal(size(alpha),size(p1),size(p2)),'alpha,p1,p2 must be same size')
%                 llh_c1 and llh_c2 should be column vectors
                llh_c1 = (repmat(log(alpha.*p1)',size(data,1),1) + (data-1)*log(1-p1'));
                llh_c2 = (repmat(log((1-alpha).*p2)',size(data,1),1) + (data-1)*log(1-p2'));
            else
                llh_c1 = log(alpha*p1) + (data-1)*log(1-p1);
                llh_c2 = log((1-alpha).*p2) + (data-1)*log(1-p2);
            end
            llh = sum(max(llh_c1,llh_c2) + log1p(exp(-abs(llh_c1-llh_c2)))); % if a>b then ln(a+b)==ln(a)+ln(1+exp(ln(b)-ln(a)))
            if nargin>4
                CDF_m = (alpha.*(1-(1-p1).^(offset)) + (1-alpha).*(1-(1-p2).^(offset)));
                llh = llh' - size(data,1)*log(1-CDF_m);
            end
        end
        function llh = computeLLH_3comp(data,alpha1,alpha2,p1,p2,p3,offset)
        % theta = [alpha, p1 (fast rate), p2 (slower rate), p3 (slowest rate)]
            % data: column vector of integers>=1
            % alpha, p1, p2, p3: can be scaler or column vector (to test many values at once)
            %                   this is mainly for plotLLH visualization
            %                between 0 and 1
        
%         ----------- copied from Geometric MixtureModelMLE: make sure is
%         properly updated: do we need to add a summation
         assert(iscolumn(data),'input data must be in column format')
            alpha3=max(1-alpha1-alpha2,0);
            llhM = zeros(numel(data),3);
            if size(alpha1,1)>1
               assert(isequal(size(alpha1),size(alpha2),size(p1),size(p2),size(p3)),'alpha1,alpha2,p1,p2,p3 must be same size')
%                 llh_c1 and llh_c2 should be column vectors
                llh_c1 = (repmat(log(alpha.*p1)',size(data,1),1) + (data-1)*log(1-p1'));
                llh_c2 = (repmat(log((1-alpha).*p2)',size(data,1),1) + (data-1)*log(1-p2'));
            else
                
                llhM(:,1) = log(alpha1*p1) + (data-1)*log(1-p1);
                llhM(:,2) = log(alpha2*p2) + (data-1)*log(1-p2);
                llhM(:,3) = log(alpha3*p3) + (data-1)*log(1-p3);
                llhM = sort(llhM,2,'descend');
                llh = llhM(:,1) + log1p(exp(llhM(:,2)-llhM(:,1)));
                llh = sum(llh + log1p(exp(llhM(:,3)-llh))); 
            end
            if nargin>6 % there's an offset
                CDF_m = (alpha1.*(1-(1-p1).^(offset)) + alpha2.*(1-(1-p2).^(offset)) + alpha3.*(1-(1-p3).^(offset)));
                 llh = llh' - size(data,1)*log(1-CDF_m);
            end
            
            %             ---------------
            
        end  
        function pdf = computePDF_2comp(data,alpha,p1,p2,offset)
        end
        function pdf = computePDF_3comp(data,alpha1,alpha2,p1,p2,p3,offset)    
        end  
        %% Helper Functions - static
        function offr = convert_prob2offrate(p)
           % this assumes offr is offrate in units of 1/frames
           offr = -log(1-p);   
        end
        function p = convert_offrate2prob(offr)
           % this assumes offr is offrate in units of 1/frames
           p = 1-exp(-offr);
        end
        function offr_frames = covert_offrateSeconds2Frames(offr_secs,framerate)
           % offr_secs is offrate in units of 1/sec
           % framerate is in units of frames/sec
           offr_frames = offr_secs./framerate;
        end
        function offr_secs = covert_offrateFrames2Seconds(offr_frames,framerate)
           % offr_secs is offrate in units of 1/sec
           % framerate is in units of frames/sec
           offr_secs = offr_frames*framerate;
        end
        function secs = convert_frames2seconds(frames,framerate)
            secs = frames/framerate;
        end
        function frames = convert_seconds2frames(secs,framerate)
            frames = secs*framerate;
        end
        function [ut,f] = makeUITablefromTable(t,f)
            if nargin<2 %no figure handle to make ui table given
             f = figure;
            end
            tc = table2cell(t);
            ut = uitable(f,'Data',tc);
            ut.ColumnName = t.Properties.VariableDescriptions;
            ut.RowName = t.Properties.RowNames;
            ut.Position(3) = ut.Extent(3);
        end
        function outax = makeResidualPlot(plotax,infax,ids1,ids2)
            % This method plots the residuals of the line objects YData values
            % within the plotax plot axes. It is mostly a obj.plotResults
            % helper.
            % plotax - the input plot axes to extract the line data from.
            %           this parameter is necessary.
            % infh - input figure handle for where to plot the residuals
            %           plot. if no input is given then makes a new figure;
            % ids1 - vector of indices from plotax.Children to use as first
            %           input for comparison
            % ids2 - vector of indices from plotax.Children to use as
            %           second input for comparison
            % Residual plot is line data from ids1 - line data from ids2
            % If no input is given then compares plotax.Children(end) - plotax.Children(#Children/2) etc.
            % Plot uses line coloring from ids1 line data

            if nargin<3
               numlines = numel(plotax.Children);
               assert(isEven(numlines),'Must having matching residual for each dataset - numlines must be even');
               ids1 = flip(numlines/2+1:numlines);
               ids2 = flip(1:numlines/2);
            end
            if nargin<2 || isempty(infax)
                figure(); outax = axes();
            else
                outax = infax;
            end  
            hold on;
            assert(numel(ids1) == numel(ids2),'Need same number of inputs for comparison');
            xs = cell(numel(ids1),1); ys = cell(numel(xs),1);
            plotcols = cell(numel(xs),1);
            displayname = cell(numel(xs),1);
            for i = 1:numel(ids1)
                id1Line = plotax.Children(ids1(i));
                id2Line = plotax.Children(ids2(i));
                plotcols{i} = id1Line.Color;
                xs{i} = id1Line.XData;
                ys{i} = id1Line.YData - interp1(id2Line.XData,id2Line.YData,id1Line.XData);
                displayname{i} = [id1Line.DisplayName ' - ' id2Line.DisplayName];
            end
            for p = 1:numel(xs)
                plot(outax,xs{p},ys{p},'Color',plotcols{p},'DisplayName',displayname{p});
            end
            legend('location','best');
            temp_xlim = outax.XLim;
            outax.XScale='log';
            outax.XLim = temp_xlim;
            outax.Title.String = 'Residual Plot';  
        end
    end
end
% add + 1 here because Matlab's geornd uses distribution that represents
% # of trials until failure and we use # of trials until success