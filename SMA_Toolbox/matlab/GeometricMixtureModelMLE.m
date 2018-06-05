 classdef GeometricMixtureModelMLE < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant=true)
        MODELS = {'1comp', '1comp_ind', '2comp', '2comp_ind', '2comp_ind_alpha', '3comp', '3comp_ind', '3comp_ind_alpha'};
    end
    properties (Constant=true, Hidden=true)
        COLORS = {[1,0,0],[0,1,0],[0,0,1],[1,0.5,0],[0.5,0,1],[0,0.5,1],...
                  [0, 0, 0], [0.5, 0, 0], [0.3, 0.7, 0.2]};
    end 
    properties
        K; % number of conditions
        Ntot; % Total number of samples
        N; % vector size:k Number of samples per-condition
        D; % data: cell array of vectors of track lengths. integer values >=1
        Stot; % Sum of the total data
        S; %vector size:k Sum of the data for each condition
        offset = 0; % all tracks are longer than this value (value in frames). ie: offset = 5 means minimun track length is 6 
        MCMCSampleNum = 10000;
        MCMCBurnFrac = 0.25;
    end

    properties (Hidden=true)
        plotDiagnostics = false;
    end
    
    methods
        function obj=GeometricMixtureModelMLE(varargin)
            if nargin>0
                obj.setData(varargin{:});
            end
        end
        
        function setData(obj,D)
            D = makecell(D);
            obj.K = numel(D);
            D=cellmap(@(d) d(:),D); %Make each D value a column vector
            allD= vertcat(D{:}); %all the data in one column
            assert(all(allD == fix(allD))); %make sure all data is an integer
            assert(all(allD>0)); % make sure all data is non-zero
            assert(all(isfinite(allD))); %make sure all data is finite
            % don't include any data that is shorter than offset
            if ~isempty(obj.offset)
                for ii = 1:numel(D)
                    D{ii}(D{ii}<=obj.offset) = [];
                end
            end
            %At this point data is OK
            obj.D = D;
            obj.N = cellfun(@numel,obj.D)';
            obj.Ntot = sum(obj.N);
            obj.S = cellfun(@sum,obj.D)';
            obj.Stot = sum(obj.S);
        end

        function nParams = modelNParam(obj, model)
            switch model
                case '1comp'
                    nParams = 1;
                case '1comp_ind'
                    nParams = obj.K;
                case '2comp'
                    nParams = 3;
                case '2comp_ind'
                    nParams = 3*obj.K;
                case '2comp_ind_alpha'
                    nParams = obj.K+2;
                case '2comp_ind_alpha_lambda_short' % alpha and lambda short free, lambda long fixed [alphaA,alphaB,lambda1A,lambda1B,lambda2]
                    nParams = 2*obj.K+1;
                case '2comp_ind_alpha_lambda_long' % alpha and lambda long free, lambda short fixed [alphaA,alphaB,lambda1,lambda2A,lambda2B]
                    nParams = 2*obj.K+1;
                case '3comp'
                    nParams = 5;
                case '3comp_ind'
                    nParams = 5*obj.K;
                case '3comp_ind_alpha'
                    nParams = 3+obj.K*2;
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        
        function nParams = modelNalphas(obj, model)
            switch model
                case '1comp'
                    nParams = 0;
                case '1comp_ind'
                    nParams = 0;
                case '2comp'
                    nParams = 1;
                case '2comp_ind'
                    nParams = obj.K;
                case '2comp_ind_alpha'
                    nParams = obj.K;
                case '2comp_ind_alpha_lambda_short'
                    nParams = obj.K;
                case '2comp_ind_alpha_lambda_long'
                    nParams = obj.K;
                case '3comp'
                    nParams = 2;
                case '3comp_ind'
                    nParams = obj.K*2;
                case '3comp_ind_alpha'
                    nParams = obj.K*2;
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        
        function nParams = modelNcomps(obj, model)
            switch model
                case '1comp'
                    nParams = 1;
                case '1comp_ind'
                    nParams = 1;
                case '2comp'
                    nParams = 2;
                case '2comp_ind'
                    nParams = 2;
                case '2comp_ind_alpha'
                    nParams = 2;
                case '2comp_ind_alpha_lambda_short'
                    nParams = 2;
                case '2comp_ind_alpha_lambda_long'
                    nParams = 2;
                case '3comp'
                    nParams = 3;
                case '3comp_ind'
                    nParams = 3;
                case '3comp_ind_alpha'
                    nParams = 3;
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        
        function desc = modelParamDesc(obj, model, paramNames)
            if nargin<3
                paramNames = {};
                for nn = 1:obj.K
                    paramNames{nn} = ' ';
                end
            end
            switch model
                case '1comp'
                    desc = cellmap(@(k) [sprintf('Lambda1_cond%i',k) '-' paramNames{k}],1:obj.K);
                case '1comp_ind'
                    desc = cellmap(@(k) [sprintf('Lambda1_cond%i',k) '-' paramNames{k}],1:obj.K);
                case '2comp'
                    desc = {'Alpha', 'Lambda1','Lambda2'};
                case '2comp_ind'
                    desc = [ cellmap(@(k) [sprintf('Alpha_cond%i',k) '-' paramNames{k}],1:obj.K),...
                             cellmap(@(k) [sprintf('Lambda1_cond%i',k) '-' paramNames{k}],1:obj.K),...
                             cellmap(@(k) [sprintf('Lambda2_cond%i',k) '-' paramNames{k}],1:obj.K)];
                case '2comp_ind_alpha'
                    desc = [ cellmap(@(k) [sprintf('Alpha_cond%i',k) '-' paramNames{k}],1:obj.K), 'Lambda1','Lambda2'];
                case '2comp_ind_alpha_lambda_short'
                    desc = [ cellmap(@(k) [sprintf('Alpha_cond%i',k) '-' paramNames{k}],1:obj.K),...
                        cellmap(@(k) [sprintf('Lambda1_cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'Lambda2'];
                case '2comp_ind_alpha_lambda_long'
                    desc = [ cellmap(@(k) [sprintf('Alpha_cond%i',k) '-' paramNames{k}],1:obj.K),...
                        'Lambda1',...
                        cellmap(@(k) [sprintf('Lambda2_cond%i',k) '-' paramNames{k}],1:obj.K)];
                case '3comp'
                    desc = {'Alpha1', 'Alpha2', 'Lambda1', 'Lambda2', 'Lambda3'};
                case '3comp_ind'
                    desc =  [ cellmap(@(k) [sprintf('Alpha1_cond%i',k) '-' paramNames{k}],1:obj.K),...
                              cellmap(@(k) [sprintf('Alpha2_cond%i',k) '-' paramNames{k}],1:obj.K),...
                              cellmap(@(k) [sprintf('Lambda1_cond%i',k) '-' paramNames{k}],1:obj.K),...
                              cellmap(@(k) [sprintf('Lambda2_cond%i',k) '-' paramNames{k}],1:obj.K),...
                              cellmap(@(k) [sprintf('Lambda3_cond%i',k) '-' paramNames{k}],1:obj.K)];
                case '3comp_ind_alpha'
                    desc =  [ cellmap(@(k) [sprintf('Alpha1_cond%i',k) '-' paramNames{k}],1:obj.K),...
                              cellmap(@(k) [sprintf('Alpha2_cond%i',k) '-' paramNames{k}],1:obj.K),...
                              'Lambda1','Lambda2','Lambda3'];
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        function tab = makeTable(obj,theta, model, tracklength,RowNames)
            if ~exist('RowNames','var')
                RowNames = cellmap(@(k) sprintf('Condition %i',k),1:obj.K);
            end
            switch model
                case '1comp'
                    VarNames = {'Lambda'};
                    Lambda = obj.convert_prob2lifetime(repmat(theta,[1 obj.K]),tracklength);
                    tab = table(Lambda','RowNames',RowNames','VariableNames',VarNames);
                case '1comp_ind'
                    VarNames = {'Lambda'};
                    Lambda = obj.convert_prob2lifetime(theta,tracklength);
                    tab = table(Lambda','RowNames',RowNames','VariableNames',VarNames);
                case '2comp'
                    VarNames = {'Alpha_Short';'Lambda_Short';'Lambda_Long'};
                    Alpha = repmat(theta(1),[1 obj.K]);     
                    Lambda1 = obj.convert_prob2lifetime(repmat(theta(2),[1 obj.K]),tracklength);
                    Lambda2 = obj.convert_prob2lifetime(repmat(theta(3),[1 obj.K]),tracklength);
                    tab = table(Alpha',Lambda1',Lambda2','RowNames',RowNames','VariableNames',VarNames);
                case '2comp_ind'
                    VarNames = {'Alpha_Short';'Lambda_Short';'Lambda_Long'};
                    Alpha = theta(1:obj.K);     
                    Lambda1 = obj.convert_prob2lifetime(theta(obj.K+1:obj.K*2),tracklength);
                    Lambda2 = obj.convert_prob2lifetime(theta(obj.K*2+1:end),tracklength);
                    tab = table(Alpha',Lambda1',Lambda2','RowNames',RowNames','VariableNames',VarNames);
                case '2comp_ind_alpha'
                    VarNames = {'Alpha_Short';'Lambda_Short';'Lambda_Long'};
                    Alpha = theta(1:obj.K);     
                    Lambda1 = obj.convert_prob2lifetime(repmat(theta(obj.K+1),[1 obj.K]),tracklength);
                    Lambda2 = obj.convert_prob2lifetime(repmat(theta(obj.K+2),[1 obj.K]),tracklength);
                    tab = table(Alpha',Lambda1',Lambda2','RowNames',RowNames','VariableNames',VarNames);
                case '3comp'
                    VarNames = {'Alpha_Short';'Alpha_Long1';'Lambda_Short';'Lambda_Long1';'Lambda_Long2'};
                    Alpha1 = repmat(theta(1),[1 obj.K]);
                    Alpha2 = repmat(theta(2),[1 obj.K]); 
                    Lambda1 = obj.convert_prob2lifetime(repmat(theta(3),[1 obj.K]),tracklength);
                    Lambda2 = obj.convert_prob2lifetime(repmat(theta(4),[1 obj.K]),tracklength);
                    Lambda3 = obj.convert_prob2lifetime(repmat(theta(5),[1 obj.K]),tracklength);
                    tab = table(Alpha1',Alpha2',Lambda1',Lambda2',Lambda3','RowNames',RowNames','VariableNames',VarNames);
                case '3comp_ind'
                    VarNames = {'Alpha_Short';'Alpha_Long1';'Lambda_Short';'Lambda_Long1';'Lambda_Long2'};
                    Alpha1 = theta(1:obj.K);
                    Alpha2 = theta(obj.K+1:obj.K*2);  
                    Lambda1 = obj.convert_prob2lifetime(theta(obj.K*2+1:obj.K*3),tracklength);
                    Lambda2 = obj.convert_prob2lifetime(theta(obj.K*3+1:obj.K*4),tracklength);
                    Lambda3 = obj.convert_prob2lifetime(theta(obj.K*4+1:obj.K*5),tracklength);
                    tab = table(Alpha1',Alpha2',Lambda1',Lambda2',Lambda3','RowNames',RowNames','VariableNames',VarNames);
                case '3comp_ind_alpha'
                    VarNames = {'Alpha_Short';'Alpha_Long1';'Lambda_Short';'Lambda_Long1';'Lambda_Long2'};
                    Alpha1 = theta(1:obj.K);
                    Alpha2 = theta(obj.K+1:obj.K*2);      
                    Lambda1 = obj.convert_prob2lifetime(repmat(theta(obj.K*2+1),[1 obj.K]),tracklength);
                    Lambda2 = obj.convert_prob2lifetime(repmat(theta(obj.K*2+2),[1 obj.K]),tracklength);
                    Lambda3 = obj.convert_prob2lifetime(repmat(theta(obj.K*2+3),[1 obj.K]),tracklength);
                    tab = table(Alpha1',Alpha2',Lambda1',Lambda2',Lambda3','RowNames',RowNames','VariableNames',VarNames);
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        function CDF = getCDF(obj,theta,model) %returns cell array of CDF values associated with each obj.D cell
            switch model
                case '1comp'
                    CDF = {};
                    p = round(repmat(theta,[1 obj.K]),3);
                    for ii = 1:obj.K
%                     CDF{ii} = 1-(1-p(ii)).^(obj.D{ii}+obj.offset);
                    CDF{ii} = 1-(1-p(ii)).^(obj.D{ii});
                    end
                case '1comp_ind'
                    CDF = {};
                    p = round(theta,3);
                    for ii = 1:obj.K
%                     CDF{ii} = 1-(1-p(ii)).^(obj.D{ii}+obj.offset);
                    CDF{ii} = 1-(1-p(ii)).^(obj.D{ii});
                    end
                case '2comp'
                    CDF = {};
                    A = repmat(round(theta(1),3),[1 obj.K]);
                    P1 = round(repmat(theta(2),[1 obj.K]),3);
                    P2 = round(repmat(theta(3),[1 obj.K]),3);
                    for ii = 1:obj.K
%                     CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset));
                    CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii})) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}));
                    end
                case {'2comp_ind','2comp_ind_alpha_lambda_short'}
                    CDF = {};
                    A = round(theta(1:obj.K),3);
                    P1 = round(theta(obj.K+1:obj.K*2),3);
                    P2 = round(theta(obj.K*2+1:end),3);
                    for ii = 1:obj.K
%                         CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset));
                        CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii})) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}));
%                         CDF{ii} = obj.compute_CDF_2comp(obj.D{ii},A(ii),P1(ii),P2(ii),obj.offset);
                    end
                case '2comp_ind_alpha' 
                    CDF = {};
                    A = round(theta(1:obj.K),3);     
                    P1 = round(repmat(theta(obj.K+1),[1 obj.K]),3);
                    P2 = round(repmat(theta(obj.K+2),[1 obj.K]),3);
                    for ii = 1:obj.K
                        CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii})) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}));
%                         CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset));
                    end
                    
                case '2comp_ind_alpha_lambda_long' 
                    CDF = {};
                    A = round(theta(1:obj.K),3);     
                    P1 = round(theta(obj.K+1),3);
                    P2 = round(theta(obj.K+2:end),3);
                    for ii = 1:obj.K
                        CDF{ii} = A(ii)*(1-(1-P1).^(obj.D{ii})) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}));
%                         CDF{ii} = A(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + (1-A(ii))*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset));
                    end
                case '3comp'
                    CDF = {};
                    A1 = repmat(round(theta(1),3),[1 obj.K]);
                    A2 = repmat(round(theta(2),3),[1 obj.K]);
                    A3 = 1-(A1+A2);
                    P1 = round(repmat(theta(3),[1 obj.K]),3);
                    P2 = round(repmat(theta(4),[1 obj.K]),3);
                    P3 = round(repmat(theta(5),[1 obj.K]),3);
                    for ii = 1:obj.K
                        %                         CDF{ii} = A1(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + A2(ii)*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset)) + ...
                        %                             A3(ii)*(1-(1-P3(ii)).^(obj.D{ii}+obj.offset));
                        CDF{ii} = A1(ii)*(1-(1-P1(ii)).^(obj.D{ii})) + A2(ii)*(1-(1-P2(ii)).^(obj.D{ii})) + ...
                            A3(ii)*(1-(1-P3(ii)).^(obj.D{ii}));
                    end
                case '3comp_ind'
                    CDF = {};
                    A1 = round(theta(1:obj.K),3);
                    A2 = round(theta(obj.K+1:obj.K*2),3);
                    A3 = 1-(A1+A2);
                    P1 = round(theta(obj.K*2+1:obj.K*3),3);
                    P2 = round(theta(obj.K*3+1:obj.K*4),3);
                    P3 = round(theta(obj.K*4+1:obj.K*5),3);
                    for ii = 1:obj.K
%                         CDF{ii} = A1(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + A2(ii)*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset)) + ...
%                             A3(ii)*(1-(1-P3(ii)).^(obj.D{ii}+obj.offset));
                        CDF{ii} = A1(ii)*(1-(1-P1(ii)).^(obj.D{ii})) + A2(ii)*(1-(1-P2(ii)).^(obj.D{ii})) + ...
                            A3(ii)*(1-(1-P3(ii)).^(obj.D{ii}));
                    end
                case '3comp_ind_alpha'
                    CDF = {};
                    A1 = round(theta(1:obj.K),3);
                    A2 = round(theta(obj.K+1:obj.K*2),3);
                    A3 = 1-(A1+A2);
                    P1 = round(repmat(theta(obj.K*2+1),[1 obj.K]),3);
                    P2 = round(repmat(theta(obj.K*2+2),[1 obj.K]),3);
                    P3 = round(repmat(theta(obj.K*2+3),[1 obj.K]),3);
                    for ii = 1:obj.K
%                         CDF{ii} = A1(ii)*(1-(1-P1(ii)).^(obj.D{ii}+obj.offset)) + A2(ii)*(1-(1-P2(ii)).^(obj.D{ii}+obj.offset)) + ...
%                             A3(ii)*(1-(1-P3(ii)).^(obj.D{ii}+obj.offset));
                           CDF{ii} = A1(ii)*(1-(1-P1(ii)).^(obj.D{ii})) + A2(ii)*(1-(1-P2(ii)).^(obj.D{ii})) + ...
                            A3(ii)*(1-(1-P3(ii)).^(obj.D{ii}));

                    end
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        
        function PDF = getPDF(obj,theta,model)
            PDF = cell(1,obj.K);
           [alpha1, alpha2,  alpha3, lambda1, lambda2, lambda3, VarNames] = getValsFromTheta(obj,theta,model);
           for ii = 1:obj.K
               if ~isempty(alpha2)
                   if ~isempty(alpha3)
              PDF{ii} = alpha1(ii)*lambda1(ii)*(1-lambda1(ii)).^(obj.D{ii}-1) + ...
                  alpha2(ii)*lambda2(ii)*(1-lambda2(ii)).^(obj.D{ii}-1) + ...
                  alpha3(ii)*lambda3(ii)*(1-lambda3(ii)).^(obj.D{ii}-1); 
                   else
                       PDF{ii} = alpha1(ii)*lambda1(ii)*(1-lambda1(ii)).^(obj.D{ii}-1) + ...
                  alpha2(ii)*lambda2(ii)*(1-lambda2(ii)).^(obj.D{ii}-1);
                   end
               else
                   PDF{ii} = alpha1(ii)*lambda1(ii)*(1-lambda1(ii)).^(obj.D{ii}-1);
               end
           end
        end
        
        function [alpha1, alpha2,  alpha3, lambda1, lambda2, lambda3, VarNames] = getValsFromTheta(obj,theta,model)
            VarNames = {'Alpha_1','Alpha_2','Alpha3','Lambda_1','Lambda_2','Lambda_3'};
            switch model
                case '1comp'
                    lambda1= theta;
                    lambda2 = [];
                    alpha1 = [];
                    alpha2 = [];
                    alpha3 = [];
                    lambda3 = [];
                case '1comp_ind'
                    lambda1= theta;
                    lambda2 = [];
                    alpha1 = [];
                    alpha2 = [];
                    alpha3 = [];
                    lambda3 = [];
                case '2comp'
                    alpha1 = theta(1);
                    alpha2 = 1-alpha1;
                    lambda1= theta(2);
                    lambda2 = theta(3); 
                    alpha3 = [];
                    lambda3 = [];
                case '2comp_ind'
                    alpha1= theta(1:obj.K);
                    alpha2 = 1-alpha1;
                    lambda1 = theta(obj.K+1:obj.K*2);
                    lambda2 = theta(obj.K*2+1:end);
                    alpha3 = [];
                    lambda3 = [];
                case '2comp_ind_alpha_lambda_short'
                    alpha1= theta(1:obj.K);
                    alpha2 = 1-alpha1;
                    lambda1 = theta(obj.K+1:obj.K*2);
                    lambda2 = theta(obj.K*2+1:end);
                    alpha3 = [];
                    lambda3 = [];
                case'2comp_ind_alpha_lambda_long'
                    alpha1= theta(1:obj.K);
                    alpha2 = 1-alpha1;
                    lambda1 = theta(obj.K+1);
                    lambda2 = theta(obj.K+2:end);
                    alpha3 = [];
                    lambda3 = [];
                case '2comp_ind_alpha'
                    alpha1 = theta(1:obj.K);
                    alpha2 = 1-alpha1;
                    lambda1 = theta(obj.K+1);
                    lambda2 = theta(obj.K+2);
                    alpha3 = [];
                    lambda3 = [];
                case '3comp'
                    alpha1 = theta(1);
                    alpha2 = theta(2);
                    alpha3 = 1-apha1-alpha2;
                    lambda1 = theta(3);
                    lambda2 = theta(4);
                    lambda3 = theta(5);
                case '3comp_ind'
                    alpha1 = theta(1:obj.K);
                    alpha2 = theta(obj.K+1:obj.K*2);
                    alpha3 = 1-apha1-alpha2;
                    lambda1 = theta(obj.K*2+1:obj.K*3);
                    lambda2 = theta(obj.K*3+1:obj.K*4);
                    lambda3 = theta(obj.K*4+1:obj.K*5);
                case '3comp_ind_alpha'
                    alpha1 = theta(1:obj.K);
                    alpha2 = theta(obj.K+1:obj.K*2);
                    alpha3 = 1-alpha1-alpha2;
                    lambda1 = theta(obj.K*2+1);
                    lambda2 = theta(obj.K*2+2);
                    lambda3 = theta(obj.K*2+3);
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        
        function [results, VarNames] = makeResultsMatrix(obj,theta,model,tracklength,std_err)
            % change to always have error
            [alpha1, alpha2, alpha3, lambda1, lambda2, lambda3, thetaNames] = getValsFromTheta(obj,theta,model);
            [errL_alpha1, errL_alpha2, errL_alpha3, errL_lambda1, errL_lambda2, errL_lambda3, thetaNames] = getValsFromTheta(obj,std_err(1,:),model);
            [errH_alpha1, errH_alpha2, errH_alpha3, errH_lambda1, errH_lambda2, errH_lambda3, thetaNames] = getValsFromTheta(obj,std_err(2,:),model);
            % assume values are 1xN
            switch model
                case '1comp'
                        VarNames = {['low_' thetaNames(4)],thetaNames(4),['high_' thetaNames(4)]};
                        % make lambdas
                        L1s = repmat([errL_lambda1;lambda1;errH_lambda1],[1 obj.K]);
                        results = round(obj.convert_prob2lifetime(L1s,tracklength),3)';
                case '1comp_ind'
                        VarNames = {['low_' thetaNames(4)],thetaNames(4),['high_' thetaNames(4)]};
                        preL1s = [errL_lambda1;lambda1;errH_lambda1];
                        L1s = round(obj.convert_prob2lifetime(L1s,tracklength),3);
                        results = L1s'
                case '2comp' %need to add in repmat
                        VarNames = {['low_' thetaNames(1)],thetaNames(1),['high_' thetaNames(1)],...
                            ['low_' thetaNames(2)],thetaNames(2),['high_' thetaNames(2)],...
                            ['low_' thetaNames(4)],thetaNames(4),['high_' thetaNames(4)],...
                            ['low_' thetaNames(5)],thetaNames(5),['high_' thetaNames(5)]};
                        a1s = [errL_alpha1;alpha1;errH_alpha1];
                        a2s = [errL_alpha2;alpha2;errH_alpha2];
                        preL1s = [errL_lambda1;lambda1;errH_lambda1];
                        preL2s = [errL_lambda2;lambda2;errH_lambda2];
                        L1s = round(obj.convert_prob2lifetime(preL1s,tracklength),3)
                        L2s = round(obj.convert_prob2lifetime(preL2s,tracklength),3)
                        results = [a1s',a2s',L1s',L2'];               
                case {'2comp_ind','2comp_ind_alpha_lambda_short','2comp_ind_alpha_lambda_long'}
                    VarNames = {['low_' thetaNames{1}],thetaNames{1},['high_' thetaNames{1}],...
                            ['low_' thetaNames{2}],thetaNames{2},['high_' thetaNames{2}],...
                            ['low_' thetaNames{4}],thetaNames{4},['high_' thetaNames{4}],...
                            ['low_' thetaNames{5}],thetaNames{5},['high_' thetaNames{5}]};
                        a1s = [errL_alpha1;alpha1;errH_alpha1];
                        a2s = [errL_alpha2;alpha2;errH_alpha2];
                        preL1s = [errL_lambda1;lambda1;errH_lambda1];
                        preL2s = [errL_lambda2;lambda2;errH_lambda2];
                        L1s = round(obj.convert_prob2lifetime(preL1s,tracklength),3)
                        L2s = round(obj.convert_prob2lifetime(preL2s,tracklength),3)
                        results = [a1s',a2s',L1s',L2s'];  
                    
                    
%                     if nargin>4
%                         VarNames = {'Alpha_Short';'+/-95%CI';'Lambda_Short';'+/-95%CI';'Lambda_Long';'+/-95%CI'};
%                         assert(numel(theta) == numel(err));
%                         Alpha = [round(theta(1:obj.K),3); round(err(1:obj.K),3)];
%                         Lambda1 = round(obj.convert_prob2lifetime(theta(obj.K+1:obj.K*2),tracklength,err(obj.K+1:obj.K*2)),3);
%                         Lambda2 = round(obj.convert_prob2lifetime(theta(obj.K*2+1:end),tracklength,err(obj.K*2+1:end)),3);
%                     else
%                         VarNames = {'Alpha_Short';'Lambda_Short';'Lambda_Long'};
%                         Alpha = round(theta(1:obj.K),3);
%                         Lambda1 = round(obj.convert_prob2lifetime(theta(obj.K+1:obj.K*2),tracklength),3);
%                         Lambda2 = round(obj.convert_prob2lifetime(theta(obj.K*2+1:end),tracklength),3);
%                     end
%                     results = [Alpha',Lambda1',Lambda2'];
                case '2comp_ind_alpha'
                     VarNames = {['low_' thetaNames{1}],thetaNames{1},['high_' thetaNames{1}],...
                            ['low_' thetaNames{2}],thetaNames{2},['high_' thetaNames{2}],...
                            ['low_' thetaNames{4}],thetaNames{4},['high_' thetaNames{4}],...
                            ['low_' thetaNames{5}],thetaNames{5},['high_' thetaNames{5}]};
                        a1s = [errL_alpha1;alpha1;errH_alpha1];
                        a2s = [errL_alpha2;alpha2;errH_alpha2];
                        preL1s = [errL_lambda1;lambda1;errH_lambda1];
                        preL2s = [errL_lambda2;lambda2;errH_lambda2];
                        L1s = repmat(round(obj.convert_prob2lifetime(preL1s,tracklength),3),1,size(a1s,2));
                        L2s = repmat(round(obj.convert_prob2lifetime(preL2s,tracklength),3),1,size(a1s,2));
                        results = [a1s',a2s',L1s',L2s']; 
%                     if nargin>4
%                         VarNames = {'Alpha_Short';'+/-95%CI';'Lambda_Short';'+/-95%CI';'Lambda_Long';'+/-95%CI'};
%                         assert(numel(theta) == numel(err));
%                         Alpha = [round(theta(1:obj.K),3); round(err(1:obj.K),3)];
%                         Lambda1 = obj.convert_prob2lifetime(repmat(theta(obj.K+1),[1 obj.K]),tracklength,repmat(err(obj.K+1),[1 obj.K]));
%                         Lambda2 = obj.convert_prob2lifetime(repmat(theta(obj.K+2),[1 obj.K]),tracklength,repmat(err(obj.K+2),[1 obj.K]));
%                     else
%                         VarNames = {'Alpha_Short';'Lambda_Short';'Lambda_Long'};
%                         Alpha = round(theta(1:obj.K),3);
%                         Lambda1 = obj.convert_prob2lifetime(repmat(theta(obj.K+1),[1 obj.K]),tracklength);
%                         Lambda2 = obj.convert_prob2lifetime(repmat(theta(obj.K+2),[1 obj.K]),tracklength);
%                     end
%                     results = [Alpha',Lambda1',Lambda2'];    
                case '3comp'
                    if nargin>4
                        VarNames = {'Alpha_1';'+/-95%CI';'Alpha_2';'+/-95%CI';'Alpha_3';'+/-95%CI';'+/-95%CI';'Lambda_1';'+/-95%CI';'Lambda_2';'+/-95%CI';'Lambda_3';'+/-95%CI'};
                        assert(numel(theta) == numel(err));
                        Alpha1 = [repmat(round(theta(1),3),[1 obj.K]); repmat(round(err(1),3),[1 obj.K])];
                        Alpha2 = [repmat(round(theta(2),3),[1 obj.K]); repmat(round(err(2),3),[1 obj.K])];
                        A3err = sqrt( repmat(round(err(1),3),[1 obj.K]).^2 + repmat(round(err(2),3),[1 obj.K]).^2 );
                        Alpha3 = [1-(Alpha1(1,:)+Alpha2(1,:)); A3err];
                        Lambda1 = round(obj.convert_prob2lifetime(repmat(theta(3),[1 obj.K]),tracklength,repmat(err(3),[1 obj.K])),3);
                        Lambda2 = round(obj.convert_prob2lifetime(repmat(theta(4),[1 obj.K]),tracklength,repmat(err(4),[1 obj.K])),3);
                        Lambda3 = round(obj.convert_prob2lifetime(repmat(theta(5),[1 obj.K]),tracklength,repmat(err(5),[1 obj.K])),3);
                    else
                        VarNames = {'Alpha_1';'Alpha_2';'Alpha_3';'Lambda_1';'Lambda_2';'Lambda_3'};
                        Alpha1 = repmat(round(theta(1),3),[1 obj.K]);
                        Alpha2 = repmat(round(theta(2),3),[1 obj.K]);
                        Alpha3 = 1-(Alpha1(1,:)+Alpha2(1,:));
                        Lambda1 = round(obj.convert_prob2lifetime(repmat(theta(3),[1 obj.K]),tracklength),3);
                        Lambda2 = round(obj.convert_prob2lifetime(repmat(theta(4),[1 obj.K]),tracklength),3);
                        Lambda3 = round(obj.convert_prob2lifetime(repmat(theta(5),[1 obj.K]),tracklength),3);
                    end
                    results = [Alpha1',Alpha2',Alpha3',Lambda1',Lambda2',Lambda3'];
                case '3comp_ind'
                    if nargin>4
                        VarNames = {'Alpha_1';'+/-95%CI';'Alpha_2';'+/-95%CI';'Alpha_3';'+/-95%CI';'+/-95%CI';'Lambda_1';'+/-95%CI';'Lambda_2';'+/-95%CI';'Lambda_3';'+/-95%CI'};
                        assert(numel(theta) == numel(err));
                        Alpha1 = round([theta(1:obj.K);err(1:obj.K)],3);
                        Alpha2 = round([theta(obj.K+1:obj.K*2);err(obj.K+1:obj.K*2)],3);
                        A3err = sqrt(err(1:obj.K).^2 + err(obj.K+1:obj.K*2).^2);
                        Alpha3 = [1-(Alpha1(1,:)+Alpha2(1,:)); A3err];
                        Lambda1 = round(obj.convert_prob2lifetime(theta(obj.K*2+1:obj.K*3),tracklength,err(obj.K*2+1:obj.K*3)),3);
                        Lambda2 = round(obj.convert_prob2lifetime(theta(obj.K*3+1:obj.K*4),tracklength,err(obj.K*3+1:obj.K*4)),3);
                        Lambda3 = round(obj.convert_prob2lifetime(theta(obj.K*4+1:obj.K*5),tracklength+err(obj.K*4+1:obj.K*5)),3);
                    else
                        VarNames = {'Alpha_1';'Alpha_2';'Alpha_3';'Lambda_1';'Lambda_2';'Lambda_3'};
                        Alpha1 = round(theta(1:obj.K),3);
                        Alpha2 = round(theta(obj.K+1:obj.K*2),3);
                        Alpha3 = 1-(Alpha1+Alpha2);
                        Lambda1 = round(obj.convert_prob2lifetime(theta(obj.K*2+1:obj.K*3),tracklength),3);
                        Lambda2 = round(obj.convert_prob2lifetime(theta(obj.K*3+1:obj.K*4),tracklength),3);
                        Lambda3 = round(obj.convert_prob2lifetime(theta(obj.K*4+1:obj.K*5),tracklength),3);
                    end
                    results = [Alpha1',Alpha2',Alpha3',Lambda1',Lambda2',Lambda3'];
                case '3comp_ind_alpha'
                    if nargin>4
                        VarNames = {'Alpha_1';'+/-95%CI';'Alpha_2';'+/-95%CI';'Alpha_3';'+/-95%CI';'Lambda_1';'+/-95%CI';'Lambda_2';'+/-95%CI';'Lambda_3';'+/-95%CI'};
                        assert(numel(theta) == numel(err));
                        Alpha1 = round([theta(1:obj.K); err(1:obj.K)],3);
                        Alpha2 = round([theta(obj.K+1:obj.K*2); err(obj.K+1:obj.K*2)],3);
                        A3err = sqrt(err(1:obj.K).^2 + err(obj.K+1:obj.K*2).^2);
                        Alpha3 = [1-(Alpha1(1,:)+Alpha2(1,:)); A3err];
                        Lambda1 = round(obj.convert_prob2lifetime(repmat(theta(obj.K*2+1),[1 obj.K]),tracklength,repmat(err(obj.K*2+1),[1 obj.K])),3);
                        Lambda2 = round(obj.convert_prob2lifetime(repmat(theta(obj.K*2+2),[1 obj.K]),tracklength,repmat(err(obj.K*2+2),[1 obj.K])),3);
                        Lambda3 = round(obj.convert_prob2lifetime(repmat(theta(obj.K*2+3),[1 obj.K]),tracklength,repmat(err(obj.K*2+3),[1 obj.K])),3);
                    else
                        VarNames = {'Alpha_1';'Alpha_2';'Alpha_3';'Lambda_1';'Lambda_2';'Lambda_3'};
                        Alpha1 = round(theta(1:obj.K),3);
                        Alpha2 = round(theta(obj.K+1:obj.K*2),3);
                        Alpha3 = 1-(Alpha1+Alpha2);
                        Lambda1 = round(obj.convert_prob2lifetime(repmat(theta(obj.K*2+1),[1 obj.K]),tracklength),3);
                        Lambda2 = round(obj.convert_prob2lifetime(repmat(theta(obj.K*2+2),[1 obj.K]),tracklength),3);
                        Lambda3 = round(obj.convert_prob2lifetime(repmat(theta(obj.K*2+3),[1 obj.K]),tracklength),3);   
                    end
                    results = [Alpha1',Alpha2',Alpha3',Lambda1',Lambda2',Lambda3'];
                otherwise
                    error('GeometricMixtureModelMLE:ModelName','Unknown model: %s', model);
            end
        end
        
        function uitab = makeUITable(obj,theta, model, tracklength,RowNames,err)
            if ~exist('RowNames','var')
                RowNames = cellmap(@(k) sprintf('Condition %i',k),1:obj.K);
            end
            
            if exist('err','var')
            [results, VarNames] = makeResultsMatrix(obj,theta,model,tracklength,err);
            else
                [results, VarNames] = makeResultsMatrix(obj,theta,model,tracklength);
            end
            uitab = uitable('Data',results,'RowName',RowNames','ColumnName',VarNames);
             tbpos=get(uitab,'Position');
             tbext=get(uitab,'Extent');
             set(uitab,'Position',[tbpos(1:2) tbext(3:4)]);
        end
        %% LLH computation
        
        function [llh, llh_conds] = LLH(obj, model, theta) %ADDED OFFSET HERE!!!
            % [in]
            % model: string with name of model. taken from obj.MODELS 
            % theta: size Nparams x Ntheta.  Ntheta is number of different
            %   thetas to evaluate (normally Ntheta=1).  Nparams is the
            %   number of parameters in theta for the given model type.
            % [out]
            %  llh - size: 1 x Ntheta.   This is the log-likelihood
            %  llh_conds - size Nconditions x Ntheta.  This is the
            %  log-likelihood for each conditions for each theta value.
            %  
            f = str2func(sprintf('LLH_%s',model));
            [llh, llh_conds] = f(obj,theta);
        end
        
        function [llh, llh_conds] = LLH_1comp(obj, theta)
            if numel(theta)==obj.K
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                llh_conds(k,:) = obj.N(k) * log(theta) + (obj.S(k)-obj.N(k)) * log(1-theta);
            end
            llh = sum(llh_conds);
        end
        
        function [llh, llh_conds] = LLH_1comp_ind(obj, theta)
            %theta should be a nParam x nTheta vector
            if numel(theta)==obj.K
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                llh_conds(k,:) = obj.N(k) * log(theta(k,:)) + (obj.S(k)-obj.N(k)) * log(1-theta(k,:));
            end
            llh = sum(llh_conds);
        end
        
        function [llh, llh_conds] = LLH_2comp(obj, theta)
            if isrow(theta) 
                theta=theta(:);
            end  
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_2comp(obj.D{k},theta(1,n), theta(2,n), theta(3,n),obj.offset);
                end
            end
            llh = sum(llh_conds); 
        end

        function [llh, llh_conds] = LLH_2comp_ind(obj, theta) %-ADDED OFFSET TO THIS ONE!!
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_2comp(obj.D{k},theta(k,n), theta(obj.K+k,n), theta(2*obj.K+k,n),obj.offset);
                end
            end
            llh = sum(llh_conds);            
        end
        function [llh, llh_conds] = LLH_2comp_ind_alpha(obj, theta)
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
%                     llh_conds(k,n) = obj.compute_LLH_2comp(obj.D{k},theta(k,n), theta(obj.K+1,n), theta(obj.K+2,n));
                    llh_conds(k,n) = obj.compute_LLH_2comp(obj.D{k},theta(k,n), theta(obj.K+1,n), theta(obj.K+2,n),obj.offset);
                end
            end
            llh = sum(llh_conds);            
        end
        
        function [llh, llh_conds] = LLH_2comp_ind_alpha_lambda_short(obj, theta) %[alphaA, alphaB, lambda1A, lambda1B, lambda2]
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_2comp(obj.D{k},theta(k,n), theta(obj.K+k,n), theta(2*obj.K+1,n),obj.offset);
                end
            end
            llh = sum(llh_conds);            
        end
        
        function [llh, llh_conds] = LLH_2comp_ind_alpha_lambda_long(obj, theta) %[alphaA, alphaB, lambda1, lambda2A, lambda2B]
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_2comp(obj.D{k},theta(k,n), theta(obj.K+1,n), theta(obj.K+1+k,n),obj.offset);
                end
            end
            llh = sum(llh_conds);            
        end
        
        function [llh, llh_conds] = LLH_3comp(obj, theta)
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_3comp(obj.D{k}, theta(1,n),theta(2,n),theta(3,n),theta(4,n),theta(5,n));                    
                end
            end
            llh = sum(llh_conds);           
        end

        function [llh, llh_conds] = LLH_3comp_ind(obj, theta)
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_3comp(obj.D{k}, theta(k,n),theta(obj.K+k,n),theta(obj.K*2+k,n),...
                                                                        theta(obj.K*3+k,n),theta(obj.K*4+k,n));                    
                end
            end
            llh = sum(llh_conds);            
        end

        function [llh, llh_conds] = LLH_3comp_ind_alpha(obj, theta) %-ADDED OFFSET TO THIS ONE!!
            if isrow(theta) 
                theta=theta(:);
            end
            ntheta = size(theta,2);
            llh_conds = zeros(obj.K, ntheta);
            for k=1:obj.K
                for n=1:ntheta
                    llh_conds(k,n) = obj.compute_LLH_3comp(obj.D{k}, theta(k,n),theta(obj.K+k,n),theta(obj.K*2+1,n),...
                                                                        theta(obj.K*2+2,n),theta(obj.K*2+3,n),obj.offset);    
                end
            end
            llh = sum(llh_conds);            
        end
        %% MCMC estimation
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC(obj,Nsamples,BurnFrac,model, fix1, fix2,fix3)
            if isempty(Nsamples)
                Nsamples = 10000;
            end
            if isempty(BurnFrac)
               BurnFrac = 0.25;
            end
            f = str2func(sprintf('MCMC_%s',model));
            if exist('fix3','var') && ~isempty(fix3) % just to make function call easier- could just set l1 and l2 to be [] and won't call fixed MCMC versions
                [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = f(obj,Nsamples,BurnFrac,fix1,fix2,fix3);
            elseif exist('fix2','var') && ~isempty(fix2) % just to make function call easier- could just set l1 and l2 to be [] and won't call fixed MCMC versions
                [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = f(obj,Nsamples,BurnFrac,fix1,fix2);
            elseif exist('fix1','var') && ~isempty(fix1)
                [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = f(obj,Nsamples,BurnFrac,fix1);
            else
                [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = f(obj,Nsamples,BurnFrac);
            end
        end
               
        function  [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_1comp(obj,Nsamples,BurnFrac) % -- how to calculate std_err?
            % first calculate init values from MLE
            theta_mle = obj.Ntot / obj.Stot;
            model='1comp';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle;
            MCMCparams.LLHfunc = @obj.LLH_1comp;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_1comp_ind(obj,Nsamples,BurnFrac) % -- how to calculate std_err?
             % first calculate init values from MLE
            theta_mle = obj.Ntot / obj.Stot;
            model='1comp_ind';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = repmat(theta_mle,1,obj.K);
            MCMCparams.LLHfunc = @obj.LLH_1comp_ind;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
       %----------- 2 component ----------------------
       %----------- ----------- ----------------------
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp(obj,Nsamples,BurnFrac,fixed_params,fixed_values)
            % fixed_params is a boolean length Nparams. 1 if fixed. 0 if free
            % fixed_values is a numeric array length Nparams. Values for each fixed param. 0 otherwise.
            
            % first calculate init values from MLE
            [theta_mle, ~, ~, ~] = MLE_2comp(obj);
            if nargin>3
               theta_mle(logical(fixed_params)) = fixed_values(logical(fixed_params)); 
               MCMCparams.fixed = fixed_params;
            end
            
            % now call MCMC
            model='2comp';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind(obj,Nsamples,BurnFrac,fixed_params,fixed_values)
            % first calculate init values from MLE
            [theta_mle, ~, ~, ~] = MLE_2comp_ind(obj);
            if nargin>3
                MCMCparams.fixed = fixed_params;
            end
            if nargin>4
               theta_mle(logical(fixed_params)) = fixed_values(logical(fixed_params)); 
            end
            % now call MCMC
            model='2comp_ind';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp_ind;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_alpha(obj,Nsamples,BurnFrac,fixed_params,fixed_values)
            % first calculate init values from MLE
            [theta_mle, ~, ~, ~] = MLE_2comp_ind_alpha(obj);
            if nargin>3
               theta_mle(logical(fixed_params)) = fixed_values(logical(fixed_params)); 
               MCMCparams.fixed = fixed_params;
            end
            % now call MCMC
            model='2comp_ind_alpha';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp_ind_alpha;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_alpha_lambda_short(obj,Nsamples,BurnFrac)
            Nparams = obj.modelNParam('2comp_ind');
            % first calculate init values from MLE
            [theta_mle_pre, ~, ~, ~] = MLE_2comp_ind_alpha_lambda_short(obj);
            lambdas_short = theta_mle_pre(obj.K+1:2*obj.K);
            lambda_long = theta_mle_pre(2*obj.K+1);
            theta_mle = [theta_mle_pre(1:obj.K);lambdas_short;lambda_long];
            % so need to repmat the lambda longs here - because just one value was optimized then input them as
            % fixed
%             fixed_params = false(1,Nparams);
%             fixed_params(2*obj.K+1:end) = true;
% %             fixed_values = zeros(1,Nparams);
% %             fixed_values(2*obj.K+1:end) = lambda_long;
%             MCMCparams.fixed = fixed_params;
             % now call MCMC
            model='2comp_ind_alpha_lambda_short';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp_ind_alpha_lambda_short;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_alpha_lambda_long(obj,Nsamples,BurnFrac)
            Nparams = obj.modelNParam('2comp_ind_alpha_lambda_long');
            % first calculate init values from MLE
            [theta_mle_pre, ~, ~, ~] = MLE_2comp_ind_alpha_lambda_long(obj);
            lambda_short = theta_mle_pre(obj.K+1);
            lambdas_long = theta_mle_pre(obj.K+2:end);
            theta_mle = [theta_mle_pre(1:obj.K);lambda_short;lambdas_long];
            % so need to repmat the lambda longs here - because just one value was optimized then input them as
            % fixed
%             fixed_params = false(1,Nparams);
%             fixed_params(obj.K+1:2*obj.K) = true;
% %             fixed_values = zeros(1,Nparams);
% %             fixed_values(obj.K+1:2*obj.K) = lambda_short;
%             MCMCparams.fixed = fixed_params;
            % now call MCMC
            model='2comp_ind_alpha_lambda_long';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp_ind_alpha_lambda_long;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
         %----------- 3 component ----------------------
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_3comp(obj,Nsamples,BurnFrac)
            % first calculate init values from MLE
            [theta_mle, ~, ~, ~] = MLE_3comp(obj);
            % now call MCMC
            model='3comp';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_3comp_ind(obj,Nsamples,BurnFrac,fixed_params,fixed_values)
            % first calculate init values from MLE
            [theta_mle, ~, ~, ~] = MLE_3comp_ind(obj);
            if nargin>3
               theta_mle(logical(fixed_params)) = fixed_values(logical(fixed_params)); 
               MCMCparams.fixed = fixed_params;
            end
            % now call MCMC
            model='3comp_ind';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_2comp_ind;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_3comp_ind_alpha(obj,Nsamples,BurnFrac,fixed_params,fixed_values)
            % first calculate init values from MLE
            [theta_mle, ~, ~, ~] = MLE_3comp_ind_alpha(obj);
            if nargin>3
                theta_mle(logical(fixed_params)) = fixed_values(logical(fixed_params));
                MCMCparams.fixed = fixed_params;
            end
            % now call MCMC
            model='3comp_ind_alpha';
            if nargin<1
                MCMCparams.Nsamples = 10000; %default
                MCMCparams.BurnFrac = 0.25;  %default
            else
                MCMCparams.Nsamples = Nsamples;
                MCMCparams.BurnFrac = BurnFrac;
            end
            MCMCparams.initTheta = theta_mle';
            MCMCparams.LLHfunc = @obj.LLH_3comp_ind_alpha;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams);
        end
        %----------- 2 component FIXED values 
        % fixed params are boolean with length = Nparams
        % fixed values are vector of param values size = size ( fixed params)
        % ------- 2 comp
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_fixShort(obj,Nsamples,BurnFrac,lambda)
            fixed_params = [0 1 0];
            fixed_values = [0 lambda 0];
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_fixLong(obj,Nsamples,BurnFrac,lambda)
            fixed_params = [0 0 1];
            fixed_values = [0 0 lambda];
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_fixBoth(obj,Nsamples,BurnFrac,lambda1,lambda2)
            % lambda1 = short component (larger value);
            fixed_params = [0 1 1];
            fixed_values = [0 lambda1 lambda2];
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp(Nsamples,BurnFrac,fixed_params,fixed_values);
        end 
        
        %-------- 2 comp independent
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_fixShort(obj,Nsamples,BurnFrac,lambda)
            Nparams = obj.modelNParam('2comp_ind');
            fixed_params = false(1,Nparams);
            fixed_params(obj.K+1:2*obj.K) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(obj.K+1:2*obj.K) = lambda;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp_ind(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_fixLong(obj,Nsamples,BurnFrac,lambda)
            Nparams = obj.modelNParam('2comp_ind');
            fixed_params = false(1,Nparams);
            fixed_params(2*obj.K+1:end) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(2*obj.K+1:end) = lambda;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp_ind(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_fixBoth(obj,Nsamples,BurnFrac,lambda1,lambda2)
            % lambda1 = short component (larger value);
            Nparams = obj.modelNParam('2comp_ind');
            fixed_params = false(1,Nparams);
            fixed_params(obj.K+1:end) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(obj.K+1:2*obj.K) = lambda1;
            fixed_values(2*obj.K+1:end) = lambda2;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp_ind(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        
           %-------- 2 comp independent alpha
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_alpha_fixShort(obj,Nsamples,BurnFrac,lambda)
            Nparams = obj.modelNParam('2comp_ind_alpha');
            fixed_params = false(1,Nparams);
            fixed_params(end-1) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(end-1) = lambda;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp_ind_alpha(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_alpha_fixLong(obj,Nsamples,BurnFrac,lambda)
            Nparams = obj.modelNParam('2comp_ind_alpha');
            fixed_params = false(1,Nparams);
            fixed_params(end) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(end) = lambda;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp_ind_alpha(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_2comp_ind_alpha_fixBoth(obj,Nsamples,BurnFrac,lambda1,lambda2)
            % lambda1 = short component (larger value);
            Nparams = obj.modelNParam('2comp_ind_alpha');
            fixed_params = false(1,Nparams);
            fixed_params(end-1:end) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(end-1) = lambda1;
            fixed_values(end) = lambda2;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_2comp_ind_alpha(Nsamples,BurnFrac,fixed_params,fixed_values);
        end
         %----------- 3 component FIXED values ----------------------
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_3comp_ind_alpha_fixBoth(obj,Nsamples,BurnFrac,lambda1,lambda2,lambda3)
            % lambda1, lambda2, lambda3 = short component first (larger value to smaller value);
            Nparams = obj.modelNParam('3comp_ind_alpha');
            fixed_params = false(1,Nparams);
            fixed_params(end-1:end) = true;
            fixed_values = zeros(1,Nparams);
            fixed_values(end-2) = lambda1;
            fixed_values(end-1) = lambda2;
            fixed_values(end) = lambda3;
            [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = obj.MCMC_3comp_ind_alpha(Nsamples,BurnFrac,fixed_params,fixed_values);          
        end
        %% MLE estimation
        %need to clean up std_err output. what is it that we  want?
        function varargout = MLE(obj, model)
            f = str2func(sprintf('MLE_%s',model));
            [varargout{1:nargout}] = f(obj);
        end
        
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_1comp(obj) % -- how to calculate std_err?
            theta_mle = obj.Ntot / obj.Stot;
            if nargout==2
                mle_llh = obj.LLH_1comp(theta_mle); % Return for convience if 2 ouputs asked for  
            end
            if nargout == 4
                mle_llh = obj.LLH_1comp(theta_mle); % Return for convience
                std_err = 1/sqrt(obj.Ntot/(theta_mle^2) + (obj.Stot-obj.Ntot)/(1-theta_mle)^2);
                std_err2 = std_err;
            end
        end
        
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_1comp_ind(obj) % -- how to calculate std_err?
            theta_mle = obj.N ./ obj.S;
            if nargout==2
                mle_llh = obj.LLH_1comp_ind(theta_mle); % Return for convience if 2 ouputs asked for
            end
            if nargout == 4
                mle_llh = obj.LLH_1comp_ind(theta_mle); % Return for convience
                std_err = 1./sqrt(obj.Ntot./(theta_mle.^2) + (obj.Stot-obj.Ntot)./(1-theta_mle).^2);
                std_err2 = std_err;
            end
        end


        function [theta_mle, mle_llh, std_err, std_err2] = MLE_2comp(obj)
            lambda1= 2*obj.Ntot / (obj.Ntot+obj.Stot);
            lambda2= (obj.Ntot/2) / (obj.Ntot/2+obj.Stot);
            theta_init=[0.5; lambda1;lambda2];
            model='2comp';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init ; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequalibn ty constraints force lambda1>=lambda2
            problem.Aineq = [0,-1,1];
            problem.bineq = 0;
            
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE('2comp',problem);
        end

        function [theta_mle, mle_llh, std_err, std_err2] = MLE_2comp_ind(obj)
            lambda1s= (2*obj.N) ./ (obj.N+obj.S);
            lambda2s= (obj.N/2) ./ (obj.N/2+obj.S);
            alphas = obj.S./obj.Stot;            
            theta_init=[alphas;lambda1s;lambda2s];
            model = '2comp_ind';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init ; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraints force lambda1>=lambda2
            A=zeros(obj.K, nParam);
            b=repmat(-1e-3,obj.K,1);
            for k=1:obj.K
                A(k,obj.K+k) = -1;
                A(k,2*obj.K+k) = 1;
            end
            problem.Aineq = A;
            problem.bineq = b;
                        
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end

        function [theta_mle, mle_llh, std_err, std_err2] = MLE_2comp_ind_alpha(obj)
            lambda1= 2*obj.Ntot / (obj.Ntot+obj.Stot);
            lambda2= (obj.Ntot/2) / (obj.Ntot/2+obj.Stot);
            alphas = obj.S./obj.Stot;
            theta_init=[alphas;lambda1;lambda2];
            model = '2comp_ind_alpha';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraint forces lambda1>=lambda2
            problem.Aineq = [zeros(1, nParam-2), -1, 1];
            problem.bineq = -1e-5;
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end
        
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_2comp_ind_alpha_lambda_long(obj)
            lambda1 = 2*obj.Ntot / (obj.Ntot+obj.Stot);
            lambda2s= (obj.N/2) ./ (obj.N/2+obj.S);
            alphas = obj.S./obj.Stot;
            theta_init=[alphas;lambda1;lambda2s];
            model = '2comp_ind_alpha_lambda_long';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraint forces lambda1>=lambda2
            A=zeros(obj.K, nParam);
            b=zeros(obj.K,1);
            for k=1:obj.K
                A(k,obj.K+1) = -1;
                A(k,obj.K+1+k) = 1;
            end
            % This is old stuff need to delete and update
%             A(k+1,obj.K+1:2*obj.K) = 1; % this makes sure that all lambda shorts are the same value?
%             b(obj.K+1,1) = obj.K;
             % for example. for obj.K = 2
%             A = [0 0 -1 0 1 0; 0 0 0 -1 0 1; 0 0 1 1 0 0]; b = [0 0 2]
%           [alphaA,alphaB,lambda1,lambda2A,lambda2B]
            problem.Aineq = A;
            problem.bineq = b;
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end
        
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_2comp_ind_alpha_lambda_short(obj)
            lambda1s= (2*obj.N) ./ (obj.N+obj.S);
            lambda2= (obj.Ntot/2) / (obj.Ntot/2+obj.Stot);%same value repeated obj.K times
            alphas = obj.S./obj.Stot;
            theta_init=[alphas;lambda1s;lambda2];
            model = '2comp_ind_alpha_lambda_short';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraint forces lambda1>=lambda2
            A=zeros(obj.K, nParam);
            b=zeros(obj.K,1);
            for k=1:obj.K
                A(k,obj.K+k) = -1;
                A(k,2*obj.K+1) = 1;
            end
             % This is old stuff need to delete and update
%             A(k+1,2*obj.K+1:end) = 1; % this makes sure that all lambda shorts are the same value?
%             b(obj.K+1,1) = obj.K;
            % for example. for obj.K = 2
%            A = [0 0 -1 0 1 0; 0 0 0 -1 0 1; 0 0 0 0 1 1]; b = [0 0 2]
            problem.Aineq = A;
            problem.bineq = b;
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end
          
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_3comp(obj)
            lambda1= 3*obj.Ntot / (2*obj.Ntot+obj.Stot);
            lambda2= 2*obj.Ntot / (obj.Ntot+obj.Stot);
            lambda3= (obj.Ntot/3) / (obj.Ntot/3+obj.Stot);
            assert(lambda1>lambda2 && lambda2>lambda3);
            theta_init=[0.3; 0.3; lambda1;lambda2;lambda3];
            model='3comp';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init ; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraints force:
            %  (1)  alpha1+alpha2<1
            %  (2)  -theta1+theta2<0
            %  (3)  -theta2+theta3<0            
            problem.Aineq = [1 1 0 0 0; 0 0 -1 1 0; 0 0 0 -1 1];
            problem.bineq = [1-eps;-eps;-eps];
           
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end
        
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_3comp_ind(obj)
            lambda1= 3*obj.N ./ (2*obj.N+obj.S);
            lambda2= 2*obj.N ./ (obj.N+obj.S);
            lambda3= (obj.N/3) ./ (obj.N/3+obj.S);
            alpha = 0.3*ones(2*obj.K,1);
            assert(all(lambda1>lambda2) && all(lambda2>lambda3));
            theta_init=[alpha; lambda1;lambda2;lambda3];
            model='3comp_ind';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init ; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraints force:
            %  (K)  alpha1+alpha2<1
            %  (K)  -theta1+theta2<0
            %  (K)  -theta2+theta3<0
            A=zeros(3*obj.K, nParam);
            b=zeros(3*obj.K,1);
            for k=1:obj.K %implement the alpha1+alpha2<1 constraints
                A(k,k) = 1;
                A(k,k+obj.K) = 1;
                b(k)=1-eps;
            end
            for k=1:obj.K %implement the lambda1>lambda2 constraints
                A(obj.K+k, 2*obj.K+k) = -1; %lambda1
                A(obj.K+k, 3*obj.K+k) = 1;  %lambda2
                b(obj.K+k)=-eps;
            end
            for k=1:obj.K %implement the lambda2>lambda3 constraints
                A(2*obj.K+k, 3*obj.K+k) = -1; %lambda2
                A(2*obj.K+k, 4*obj.K+k) = 1;  %lambda3
                b(2*obj.K+k)=-eps;
            end
            problem.Aineq = A;
            problem.bineq = b;
           
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end
        
        function [theta_mle, mle_llh, std_err, std_err2] = MLE_3comp_ind_alpha(obj)
            lambda1= 3*obj.Ntot / (2*obj.Ntot+obj.Stot);
            lambda2= 2*obj.Ntot / (obj.Ntot+obj.Stot);
            lambda3= (obj.Ntot/3) / (obj.Ntot/3+obj.Stot);
            alpha = 0.3*ones(2*obj.K,1);
            assert(lambda1>lambda2 && lambda2>lambda3);
            theta_init=[alpha; lambda1;lambda2;lambda3];
            model='3comp_ind_alpha';
            nParam = obj.modelNParam(model);
            problem.objective =  @(theta) -obj.LLH(model,theta);
            problem.x0 = theta_init ; %starting guess
            problem.lb = zeros(nParam,1)+eps;
            problem.ub = ones(nParam,1)-eps;
            %Linear inequality constraints force:
            %  (K)  alpha1+alpha2<1
            %  (1)  -theta1+theta2<0
            %  (1)  -theta2+theta3<0
            A=zeros(obj.K+2, nParam);
            b=zeros(obj.K+2,1);
            for k=1:obj.K %implement the alpha1+alpha2<1 constraints
                A(k,k) = 1;
                A(k,k+obj.K) = 1;
                b(k)=1-eps;
            end
            %implement the lambda1>lambda2 constraint
            A(end-1,end-2) = -1;
            A(end-1,end-1) = 1;
            b(end-1) = -eps;
            %implement the lambda2>lambda3 constraint
            A(end,end-1) = -1;
            A(end,end) = 1;
            b(end) = -eps;

            problem.Aineq = A;
            problem.bineq = b;
           
            [theta_mle, mle_llh, std_err, std_err2] = obj.fmincon_MLE(model,problem);
        end
        
        %% Model selction
        function [selected_model, theta_mle, model_aic] = AICTest(obj,models)
            if nargin<2
                models=obj.MODELS;
            end
            nModels = numel(models);
            assert(nModels>=2);
            theta_mle = cell(1,nModels);
            mle_llh = zeros(1,nModels);
            model_aic = zeros(1,nModels);
            for n=1:nModels
                [theta_mle{n}, mle_llh(n)] = obj.MLE(models{n});
                nParams = obj.modelNParam(models{n});
                model_aic(n) = 2*nParams - 2*mle_llh(n);
            end
            [~,selected] = min(model_aic);
            selected_model = models{selected};
            theta_mle = theta_mle{selected};
        end

        function [selected_model, theta_mle, model_bic, std_err, std_err2, alltheta] = BICTest(obj,models)
            if nargin<2
                models=obj.MODELS;
            end
            nModels = numel(models);
            assert(nModels>=2);
            alltheta = cell(1,nModels); %--- for returning all thetas later
            theta_mle = cell(1,nModels);
            std_err = cell(1,nModels); %-- should this be cell or array?
            std_err2 = cell(1,nModels);
            mle_llh = zeros(1,nModels);
            model_bic = zeros(1,nModels);
            for n=1:nModels
                [theta_mle{n}, mle_llh(n), std_err{n}, std_err2{n}] = obj.MLE(models{n});
                nParams = obj.modelNParam(models{n});
                model_bic(n) = -2*mle_llh(n) + nParams * log(obj.Ntot);
            end
            [~,selected] = min(model_bic);
            selected_model = models{selected};
            alltheta = theta_mle;
            theta_mle = theta_mle{selected}; %-- only returning best MLE theta from BIC?
        end
        
        function BF = BayesFactor(~,llh_m1, llh_m0, p_m1, p_m0)
            % BF1,0 = [p(M1|D)/p(M0|D]/[p(M1)/p(M0)]
            if nargin == 3
                % no priors
                p_m1 = 1; p_m0 = 1;
            end
            BF = exp(  llh_m1 - llh_m0 - log(p_m1) + log(p_m0) ); 
        end
        
        function [BF, llh_m1, llh_m0] = callBayesFactor(obj,models)
           assert(numel(models) == 2, 'Bayes Factor only compares between 2 models');
            Nsamples = 20000;
            BurnFrac = 0.25;
            [samples, ~, ~, ~] = obj.MCMC(Nsamples,BurnFrac,models{1});
            llh_m1 = obj.LLH(models{1}, mean(samples));
            [samples, ~, ~, ~] = obj.MCMC(Nsamples,BurnFrac,models{2});
            llh_m0 = obj.LLH(models{2}, mean(samples));
            
            BF = obj.BayesFactor(llh_m1, llh_m0);
        end
        
        function [ModelProbabilities,selectedModel,LogEvidence_Array,bestTheta,allsamples] = evidenceTest(obj,models)
            Nsamples = 20000;
            BurnFrac = 0.25;
            LogEvidence_Array = zeros(1,length(models));
            allsamples = cell(1,numel(models));
            for nn=1:length(models)
                [samples, ~, ~, ~] = obj.MCMC(Nsamples,BurnFrac,models{nn});
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
        end
        
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
             p = prod(1./theta(na+1:end)); %prior on all lambda is 1/lambda --- but this is probabilities so actually need to fix this !!!!!!
         end
        
        %% Plotting        
        function plotLLH(obj,model)
            if nargin<2
                nmodels = numel(obj.MODELS);
                for n=1:nmodels
                    obj.plotLLH_ndims(obj.MODELS{n})
                end
            else
                obj.plotLLH_ndims(model)
            end
        end
        



        %% Simulation
        % -- theta in the form: 
        %           [alpha1,alpha2,lambdashort1,lambdashort2,lambdalong1,lambdalong2]
        %    always should have lambdas in decresaing values. (longer components last)
        function simulate(obj, model, Nsim, theta)
            f = str2func(sprintf('simulate_%s',model));
            f(obj,Nsim,theta);
        end
        
        function simulate_1comp(obj, Nsim,theta)
            assert(0<theta && theta<1);
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            obj.setData(cellmap(@(n) geornd(theta,1,n)+1, Nsim));
        end
        
        function simulate_1comp_ind(obj, Nsim, theta)
            assert(numel(theta)==obj.K);
            assert(all(0<theta) && all(theta<1));
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            obj.setData(cellmap(@(k) geornd(theta(k),1,Nsim(k))+1, 1:obj.K));
        end
        
        function simulate_2comp(obj, Nsim, theta)
            assert(all(0<theta) && all(theta<1));
            assert(theta(2)>=theta(3)); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                pop1 = sum(rand(1,Nsim(k))<theta(1));
                newD{k} = 1+[geornd(theta(2),1,pop1), geornd(theta(3),1,Nsim(k)-pop1)];
            end
            obj.setData(newD);
        end

        function simulate_2comp_ind(obj, Nsim, theta)
            assert(all(0<theta) && all(theta<1));
            assert(numel(theta) == obj.modelNParam('2comp_ind'));
            assert(all(theta(obj.K+1:2*obj.K)>theta(2*obj.K+1:end))); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                pop1 = sum(rand(1,Nsim(k))<theta(k));
                newD{k} = 1+[geornd(theta(obj.K+k),1,pop1), geornd(theta(2*obj.K+k),1,Nsim(k)-pop1)];
            end
            
            obj.setData(newD);
        end

        function simulate_2comp_ind_alpha_lambda_short(obj, Nsim, theta)
            assert(all(0<theta) && all(theta<1));
            assert(numel(theta) == obj.modelNParam('2comp_ind'));
            assert(all(theta(obj.K+1:2*obj.K)>theta(2*obj.K+1:end))); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                pop1 = sum(rand(1,Nsim(k))<theta(k));
                newD{k} = 1+[geornd(theta(obj.K+k),1,pop1), geornd(theta(2*obj.K+k),1,Nsim(k)-pop1)];
            end
            obj.setData(newD);
        end
        
        function simulate_2comp_ind_alpha_lambda_long(obj, Nsim, theta)
            assert(all(0<theta) && all(theta<1));
            assert(numel(theta) == obj.modelNParam('2comp_ind_alpha_lambda_long'));
            assert(all(theta(obj.K+1)>theta(2*obj.K+1:end))); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                pop1 = sum(rand(1,Nsim(k))<theta(k));
                newD{k} = 1+[geornd(theta(obj.K+1),1,pop1), geornd(theta(obj.K+1+k),1,Nsim(k)-pop1)];
            end
            obj.setData(newD);
        end
        
        function simulate_2comp_ind_alpha(obj, Nsim, theta)
            assert(all(0<theta) && all(theta<1));
            assert(numel(theta) == obj.modelNParam('2comp_ind_alpha'));
            assert(theta(end-1)>=theta(end)); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                pop1 = sum(rand(1,Nsim(k))<theta(k));
                newD{k} = 1+[geornd(theta(obj.K+1),1,pop1), geornd(theta(obj.K+2),1,Nsim(k)-pop1)];
%                 newD{k}=newD{k}-obj.offset;
%                 newD{k}(newD{k}<1) = [];
            end
            
            obj.setData(newD);
        end

        function simulate_3comp(obj, Nsim, theta)
            %components are [Alpha1, Alpha2, Lambda1, Lambda2, Lambda3].
            %Lambdas should be in decreasing order.  Alpha1+Alpha2<1;
            assert(all(0<theta) && all(theta<1));
            assert(theta(1)+theta(2) < 1); %check alphas don't add to more than 1.
            assert(theta(3)>=theta(4) && theta(4) >= theta(5)); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                samp = rand(1,Nsim(k));
                pop1 = sum(samp<theta(1));
                pop2 = sum(theta(1)<=samp & samp< (theta(1)+theta(2)));
                pop3 = Nsim(k)-pop1-pop2;
                newD{k} = 1+[geornd(theta(3),1,pop1), geornd(theta(4),1,pop2), geornd(theta(5),1,pop3)];
            end
            obj.setData(newD);
        end
        
        function simulate_3comp_ind(obj, Nsim, theta)
            %components are [Alpha1_1,... Alpha1_k, Alpha2_1,... Alpha2_k, Lambda1_1,... Lambda1_k, Lambda2_1,... Lambda2_k,Lambda3_1,... Lambda3_k].
            %Lambdas should be in decreasing order for each K.  Alpha1_k+Alpha2_k<1;
            assert(all(0<theta) && all(theta<1));
            theta = reshape(theta,[],5); %now theta is obj.K x 5  and each row is a condition's theta
            assert(all(theta(:,1)+theta(:,2) < 1)); %check alphas don't add to more than 1.
            assert(all(theta(:,3)>=theta(:,4)) && all(theta(:,4)>=theta(:,5))); %always should have lambdas in decresaing values. (longer components last)
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                samp = rand(1,Nsim(k));
                pop1 = sum(samp<theta(k,1));
                pop2 = sum(theta(k,1)<=samp & samp<theta(k,1)+theta(k,2));
                pop3 = Nsim(k)-pop1-pop2;
                newD{k} = 1+[geornd(theta(k,3),1,pop1), geornd(theta(k,4),1,pop2), geornd(theta(k,5),1,pop3)];
            end
            obj.setData(newD);
        end
        
        function simulate_3comp_ind_alpha(obj, Nsim, theta)
            %components are [Alpha1_1,... Alpha1_k, Alpha2_1,... Alpha2_k, Lambda1, Lambda2, Lambda3].
            assert(all(0<theta) && all(theta<1));
            %Lambdas should be in decreasing order.  Alpha1_k+Alpha2_k<1;
            alphas = reshape(theta(1:2*obj.K),[],2); % size obj.Kx2 each row are alphas for that condition
            lambdas = theta(end-2:end); %all lambdas are the same
            assert(all(alphas(:,1)+alphas(:,2) < 1)); %check alphas don't add to more than 1.
            assert(lambdas(1)>=lambdas(2) && lambdas(2)>=lambdas(3));
            if isscalar(Nsim)
                Nsim=Nsim *ones(1,obj.K);
            end
            newD=cell(1,obj.K);
            for k=1:obj.K
                samp = rand(1,Nsim(k));
                pop1 = sum(samp<alphas(k,1));
                pop2 = sum(alphas(k,1)<=samp & samp<alphas(k,1)+alphas(k,2));
                pop3 = Nsim(k)-pop1-pop2;
                newD{k} = 1+[geornd(lambdas(1),1,pop1), geornd(lambdas(2),1,pop2), geornd(lambdas(3),1,pop3)];
            end
            obj.setData(newD);
        end
        
        %% --- plotting methods
        function plotMCMC(obj,model,samples, samplesLLH, proposed_samples, proposed_samplesLLH, framerate)
            Nparams = obj.modelNParam(model);
            desc = obj.modelParamDesc(model);
            Nalphas = obj.modelNalphas(model);
            [mean_vals, std_error] = obj.calcMCMCstats(samples, model, framerate); 
            nrows = 3;
            ncols = ceil(Nparams/nrows);
            Nsamples = size(samples,1);
            f1 = figure;
            f1.Position = [10 10 900 600]; 
            xmax = max(max(samples(:,Nalphas+1:end)));
            xmin = min(min(samples(:,Nalphas+1:end)));
            for n = 1:Nparams
                subplot(nrows,ncols,n);
                %check if fixed
                if all(samples(:,n) == repmat(samples(1,n),Nsamples,1));
                    line(repmat(samples(1,n),Nsamples,1),(1:Nsamples));
                    leg_string = sprintf('Mean: %.3g\n',mean_vals(n));
                else
                histogram(samples(:,n));
                if n>Nalphas
                xlim([xmin xmax]);
                end
                leg_string = sprintf('Mean: %.3g\n 95CI: [%.3g - %.3g]\n',mean_vals(n), std_error(1,n),std_error(2,n));
                end
                legend(leg_string, 'FontSize',6);
                title(strrep(desc{n},'_','-'));
            end
        end
        
        function plotMCMCchain(obj,model,samples, samplesLLH, proposed_samples, proposed_samplesLLH, framerate)
            Nparams = obj.modelNParam(model);
            desc = obj.modelParamDesc(model);
            f1 = figure;
            f1.Position = [15 15 900 600]; 
            nrows = 3;
            ncols = ceil(Nparams/nrows);
            for n = 1:Nparams
                subplot(nrows,ncols,n);
                plot(proposed_samples(:,n),'o'); hold on; plot(samples(:,n),'*');
                xlabel('Iteration'); ylabel('Sample'); legend({'Proposed Values', 'Accepted Values'});
                title(strrep(desc{n},'_','-'));
            end
        end
        
        function [theta, std_er] = calcMCMCstats(obj, samples, model, framerate)
            %delete out commented sections - old crap
            Nalphas = obj.modelNalphas(model);
            Nsamples = size(samples,1);
            Nparams = size(samples,2);
            fixed = all(samples == repmat(samples(1,:),Nsamples,1)); % boolean row vector
            theta = mean(samples);
            theta(Nalphas+1:end) = obj.convert_prob2lifetime(theta(Nalphas+1:end),framerate);
%             std_error = zeros(1, Nparams);
%             std_error(~fixed)= std(samples(:,~fixed));
            a = 0.682;
            lower = (1-a)/2; upper = 1-((1-a)/2);
            std_er = quantile(samples,[lower, upper],1);
            std_er(:,Nalphas+1:end) = obj.convert_prob2lifetime(std_er(:,Nalphas+1:end),framerate);
% %             [f,x] = arrayfun(@(x) ecdf(samples(:,x)), 1:size(samples,2),'UniformOutput',false);
% %             ci95_upper = cellfun(@(i,j) interp1(i,j,0.975), f,x);
% %             ci95_lower = cellfun(@(i,j) interp1(i,j,0.025), f,x);   
% %             ci95_upper = obj.convert_prob2lifetime(ci95(2,:),framerate);
% %             ci95_lower = obj.convert_prob2lifetime(ci95(1,:),framerate);
% %             ci95 = max([ci95_lower;ci95_upper], 0);
        end  
        
         function [theta, std_er] = calcMCMCstats_prob(obj, samples, model)
            %delete out commented sections - old crap
            Nalphas = obj.modelNalphas(model);
            Nsamples = size(samples,1);
            Nparams = size(samples,2);
            fixed = all(samples == repmat(samples(1,:),Nsamples,1)); % boolean row vector
            theta = mean(samples);
%             std_error = zeros(1, Nparams);
%             std_error(~fixed)= std(samples(:,~fixed));
            a = 0.682;
            lower = (1-a)/2; upper = 1-((1-a)/2);
            std_er = quantile(samples,[lower, upper],1);
%             [f,x] = arrayfun(@(x) ecdf(samples(:,x)), 1:size(samples,2),'UniformOutput',false);
%             ci95_upper = cellfun(@(i,j) interp1(i,j,0.975), f,x);
%             ci95_lower = cellfun(@(i,j) interp1(i,j,0.025), f,x);
%             ci95 = max([ci95_lower;ci95_upper], 0);
        end
     
        function bootStrap_singleDataset(obj,rat,numIter)
            assert(numel(obj.D) == 1, 'Obj can only have one data set for this call');
            newD = cell(numIter,1);
            numrand = numel(obj.D{1});
            for ii = 1:numIter
                rvals=rand(numrand,1);
                [~,ids]=sort(rvals);
                newids = ids(1:ceil(numrand*rat));
                newD{ii} = obj.D{1}(newids);
            end
            obj.setData(newD')
        end
        
        function fig = plot_Data_Fit(obj,theta,model,frameTime, condNames,err)
            % cdf_data = 1 - (1 - theta).^ obj.D{1};
            %  cdf_data = 1 - theta(1).*(1-theta(2)).^ obj.D{1} + (1-theta(1)).*(1-theta(3)).^ obj.D{1};
            if nargin<4
                frameTime = 1;
            end
            if nargin<5 %make empty conditionNames
                for cc = numel(obj.D)
                    condNames{cc} = '';
                end
            end
            fig = figure('Position',[-1700, 300, 2700 600]);
            axH = axes();
            if obj.K>9
                cols = cell(1,obj.K);
                cm = lines(obj.K);
                for cc = 1:size(cm,1)
                    cols{cc} = cm(cc,:);
                end
            else
                cols = obj.COLORS;
            end
            sim_gmm = GeometricMixtureModelMLE();
            sim_gmm.K = numel(obj.D);
            sim_gmm.Ntot = 50000;
            sim_gmm.offset = obj.offset;
            cdf_data = [];
            pdf_data = [];
%             cdf_data = obj.getCDF(theta,model);
%             pdf_data = obj.getPDF(theta,model);
            residx = {}; residy = {};
            modelstring = strrep(model,'_','-');
            legs = cellmap(@(x,y) strrep(x,'_','-'), obj.modelParamDesc(model,condNames));
%             if strcmp(model,'1comp') || strcmp(model,'2comp') || strcmp(model,'3comp')
%                 legs = repmat(legs,[numel(condNames) 1]);
%                 legs = cellmap(@(x,y) [x '-' y], legs,condNames');
%             end
            sim_gmm.simulate(model, sim_gmm.Ntot, theta);
%             --- delete this: fixed setData to include obj.offset
%             newD = sim_gmm.D;
%             for ii = 1:sim_gmm.K
%                 newD{ii} = newD{ii} - sim_gmm.offset;
%                 newD{ii}(newD{ii}<1)=[];
%             end
%             newD{k}=newD{k}-obj.offset;
%                 newD{k}(newD{k}<1) = [];
%             newdat = cellmap(@(x) x+sim_gmm.offset, sim_gmm.D)
%             sim_gmm.setData(newD);
%           --- delete above ----------------
            llval = sim_gmm.LLH(model,theta);
            % -------plot data over fit for CDF
            subplot(axH);
            s1 = subplot(1,5,1);
            hold('on');
            for rr = 1:obj.K
                [x,f] = ecdf((obj.D{rr}).*frameTime); %-- removed offset addition
                H = plot(f(2:end),x(2:end));
                hold('on');
                H.DisplayName = ['Data with ' num2str(numel(obj.D{rr})) ' tracks'];
                H.LineWidth = 1.5;
                H.Color = cols{rr};
                H.DisplayName = legs{rr}; 
                simdata =   (sim_gmm.D{rr})*frameTime; %-- removed offset addition
                [y,g] = ecdf(simdata);
                I = plot(g(2:end),y(2:end));
                I.DisplayName = 'Sim';
                I.LineWidth = 1.5;
                I.Color = [0 0 0];
                I.LineStyle = '--';
                I.DisplayName = [legs{rr} ' Sim Results'];
                if ~isempty(cdf_data) 
                    C = plot((obj.D{rr}).*frameTime',cdf_data{rr}', '.'); %-- removed offset addition
                    C.DisplayName = 'CDF';
                    C.LineWidth = 1;
                    C.MarkerSize = 5;
                    C.Color = [0 0 1];
                end
                % save values for residual plot
                residx{rr} = f(2:end);
                residy{rr} = x(2:end) - interp1(g(2:end),y(2:end),f(2:end));
            end
            ylabel('Cumulative Probability');
            if nargin>=4
                xlabel('Duration (Secs)');
            elseif nargin<4
                xlabel('Duration (frames)');
            end
            legend('location','best');
            temp_xlim = s1.XLim;
            s1.XScale='log';
            s1.XLim = temp_xlim;
            s1.YScale='linear';
            s1.Title.String = sprintf('CDF of Track Lengths. %s Model Sim. LLH =  %.3g', modelstring, llval);
            
             % -------plot data over fit for histogram
            s2 = subplot(1,5,2);
            hold('on');
            clear edges;
            for rr = 1:numel(obj.D)
                if rr == 1
                [N,edges] = histcounts((obj.D{rr}+obj.offset).*frameTime,'Normalization','probability');
                else
                    [N,edges] = histcounts((obj.D{rr}+obj.offset).*frameTime,edges,'Normalization','probability');
                end
                hold('on');
                binwidth = diff(edges(1:2));
                hst = plot(edges(2:end)-binwidth,N);
%                 HH.DisplayName = 'Data';
                hst.LineWidth = 1.5;
                hst.LineStyle = '-';
                hst.Color = cols{rr};
                hst.DisplayName = legs{rr};
                simdata =   (sim_gmm.D{rr}+obj.offset)*frameTime;
                [NI,edgesI] = histcounts(simdata,edges,'Normalization','probability');
                binwidthI = diff(edgesI(1:2));
                simhst = plot(edgesI(2:end)-binwidthI,NI);
%                 II.DisplayName = 'Sim';
                simhst.LineWidth = 1.5;
                simhst.Color = [0 0 0];%cols{rr};
                simhst.LineStyle = ':';
                simhst.DisplayName = [legs{rr} ' Sim Results'];
                 if ~isempty(pdf_data) 
                    C = plot((obj.D{rr}).*frameTime',pdf_data{rr}', '.'); %-- removed offset addition
                    C.DisplayName = 'PDF';
                    C.LineWidth = 1;
                    C.MarkerSize = 5;
                    C.Color = [0 0 1];
                end
            end
            ylabel('Probability Density');
            if nargin>=4
                xlabel('Duration (Secs)');
            elseif nargin<4
                xlabel('Duration (frames)');
            end
            legend('location','best');
            temp_xlim = s2.XLim;
            s2.XScale='log';
            s2.XLim = temp_xlim;
            s2.YScale='log';
            s2.Title.String = sprintf('PDF of Track Lengths. %s Model Sim. LLH =  %.3g', modelstring, llval);
            
            
            
            % ---------now plot residuals
            s3 = subplot(1,5,3);
            hold('on');
            for rr = 1:numel(obj.D)
                R = plot(residx{rr},residy{rr});
                R.Color = cols{rr};
                R.DisplayName = legs{rr};
            end
            ylabel('CDF Difference between Dat and Simulation');
            if nargin>=4
                xlabel('Duration (Secs)');
            elseif nargin<4
                xlabel('Duration (frames)');
            end
            legend('location','best');
            temp_xlim = s3.XLim;
            s3.XScale='log';
            s3.XLim = temp_xlim;

            s3.Title.String = sprintf('Residuals From %s Fit to Data. LLH =  %.3g', modelstring, llval);   
            % -----------now make uitable
            s4 = subplot(1,5,4);
            pos = get(subplot(1,5,4),'position');
            if nargin == 4
                uitab = obj.makeUITable(theta, model, frameTime);
            elseif nargin == 5
                uitab = obj.makeUITable(theta, model, frameTime,condNames);
            elseif nargin == 6
                uitab = obj.makeUITable(theta, model, frameTime,condNames,err);
            end
            set(uitab,'units','normalized');
            set(uitab,'position',pos);
            set(uitab,'FontSize',12);
            uitab.Position(3) = uitab.Extent(3); %make sure it's long enough in x
            s4.Title.String = 'Table of Fit Values (s)';
            axis off
        end
    end % public methods

    methods (Access=protected)
        function plotLLH_ndims(obj,model)
            M=1e4; %number of samples
            xs = logspace(-3,0,M);
            [mle, mle_llh] =  obj.MLE(model);
            nParam = obj.modelNParam(model);
            desc = obj.modelParamDesc(model);
            figure();
            hold on;
            for n=1:nParam
                thetas = repmat(mle,1,M);
                thetas(n,:)=xs;
                llh = obj.LLH(model, thetas);
                c = obj.COLORS{mod(n,numel(obj.COLORS))+1};
                name = sprintf('Theta(%i):%s',n,desc{n});
                plot(xs, llh,'Color',c,'DisplayName',name);
                name=sprintf('MLE(%i):%s =%.3g',n,desc{n},mle(n));
                plot(mle(n), mle_llh,'Color',c,'Marker','*','DisplayName',name);
            end
            set(gca(),'XScale','log','YScale','log');
            xlabel('Lambda1');
            ylabel('LLH');
            title(sprintf('MLE Estimate: %s',model),'Interpreter','none');
            H=legend('Location','best');
            H.Interpreter='None';
        end
        
        function [samples, samplesLLH, proposed_samples, proposed_samplesLLH] = MCMC_sample(obj,model,MCMCparams)
            Nparams = obj.modelNParam(model);
            Nalphas = obj.modelNalphas(model);
            Ncomp = obj.modelNcomps(model);
            Nsamples = MCMCparams.Nsamples;
            Ntot = ceil(Nsamples*(1+MCMCparams.BurnFrac)); %#ok<*PROP>
            samples = zeros(Ntot, Nparams);
            samples(1,:) = MCMCparams.initTheta;
            samplesLLH = zeros(Ntot,1);
            samplesLLH(1) = MCMCparams.LLHfunc(samples(1,:));
            alpha_sigma = 0.05;
            lambda_sigma = 0.05*MCMCparams.initTheta(Nalphas+1:Nparams);
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
                if isfield(MCMCparams,'fixed') && MCMCparams.fixed(selected_param)
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
                    [alpha1, alpha2, alpha3, lambda1, lambda2, lambda3] = obj.getValsFromTheta(propTheta,model);
%              delete       Nind_lambda = (Nparams-Nalphas)/Ncomp; % number of lamda_1
                    if Ncomp == 1  % determine if lambdas are in correct order
                        lambda_invalid = false;
                    elseif Ncomp == 2
%              delete:           lambda_invalid = any( propTheta(Nalphas+1:Nalphas+Nind_lambda) < propTheta(Nalphas+Nind_lambda+1:end) );
                        lambda_invalid = any(lambda1<lambda2);
                    elseif Ncomp == 3 %% hmmmmm is something wrong here/!!!???
%              delete:           lambda_invalid = any( propTheta(Nalphas+1:Nalphas+Nind_lambda) < propTheta(Nalphas+Nind_lambda+1:Nalphas+2*Nind_lambda) ) ||...
%                             any( propTheta(Nalphas+Nind_lambda+1:Nalphas+2*Nind_lambda)  < propTheta(Nalphas+2*Nind_lambda+1:Nalphas+3*Nind_lambda)) ;
                        lambda_invalid = any(lambda1<lambda2) || any(lambda2<lambda3);
                    end
                    
                    invalid =  lambda_invalid || any(propTheta<0) || any(propTheta(1:Nalphas)>1) ; % alpha must be betwee 0 and 1
                    
                    % if 3 components - make sure alphas sum to less than or equal to 1
                    if Ncomp == 3
                        alpha_invalid = any(propTheta(1:Nalphas/2)+propTheta(Nalphas/2+1:Nalphas) > 1);
                        invalid = invalid || alpha_invalid;
                    end
                end
                % now if not valid continue, otherwise if valid then calculate the LLH
                if invalid
                    proposed_samplesLLH(k) = -inf;
                    Ninvalid = Ninvalid+1;
                else
                    propThetaLLH = MCMCparams.LLHfunc(propTheta);
                    proposed_samplesLLH(k) = propThetaLLH;
                end
                
                % given that it is valid, decided to accept or reject
                if ~invalid && (propThetaLLH>samplesLLH(k-1) || log(rand(1)) < propThetaLLH - samplesLLH(k-1)) % accept condition
                    samples(k,:)= propTheta;
                    if selected_param<=Nalphas
                        assert(samples(k,selected_param)<1);
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

        function [theta_mle, mle_llh, std_err, std_err2] = fmincon_MLE(obj,model, problem)
            max_iter=500;
            nParam = obj.modelNParam(model);
            sequence=zeros(nParam,max_iter+1);
            llh = zeros(1,max_iter+1);
            sequence(:,1)=problem.x0;
            llh(1)=obj.LLH(model,problem.x0);
            function stop=output(theta, opt, state)
                sequence(:,opt.iteration+1)=theta;
                llh(opt.iteration+1)=-opt.fval;
                stop=strcmp(state,'done');
            end
            
            opts = optimoptions('fmincon');
            opts.Algorithm = 'interior-point';
%             opts.Algorithm = 'sqp';
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
            problem.solver = 'fmincon';
            problem.options = opts;
            fprintf('Optimizing Model: %s\n', model);
            [theta_mle, mle_llh, exitFlag, out, ~, grad, hess] = fmincon(problem);
            %optimize with fmin con
            mle_llh = -mle_llh;
            fprintf('---> MLE_LLH: %g\n',mle_llh);
%             fprintf('ExitFlag: %i\n',exitFlag);
            std_err = sqrt(diag(pinv(hess)));
            std_err2 = sqrt(1./diag(hess));
           % fprintf('SE(theta): %s',mat2str(sqrt(1./diag(hess)))); %-- old output, not as good 
            fprintf('SE(theta): %s',std_err2);
            niter = out.iterations;
            desc=obj.modelParamDesc(model);
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
                title(sprintf('MLE Optimization Report: %s', model),'Interpreter','none')
                subplot(2,1,2);
                plot(xs,llh(xs),'Color','k','DisplayName','LLH');
                xlabel('Iteration');
                ylabel('LLH');
                H=legend('Location','Best');
                H.Interpreter='None';
            end
        end
    end % protected methods

    methods (Static=true)
%         function llh = compute_LLH_2comp(data,alpha,p1,p2,offset)
%             % alpha, p1, p2 are scalars between 0 and 1
%             % data is column vector of integers>=1
% %             assert(nargin==5);
%             assert(iscolumn(data))
%             if nargin>4
%                 llh_c1 = log(alpha*p1) + (data+offset-1)*log(1-p1);
%                 llh_c2 = log((1-alpha)*p2) + (data+offset-1)*log(1-p2);
%                 % this is the probability that data<offset so = CDF
% %                 CDF = GeometricMixtureModelMLE.compute_CDF_2comp(data,alpha,p1,p2,offset);
%                 CDF = GeometricMixtureModelMLE.compute_CDF_2comp(offset,alpha,p1,p2);
%                 % this adds in a correction factor (LLH/(1-CDF(offset))) to account for the fact that only tracks > offset
%                 % are included in the data
%                 llh = sum(max(llh_c1,llh_c2) + log1p(exp(-abs(llh_c1-llh_c2))))-sum(numel(data)*log(1-CDF));% if a>b then ln(a+b)==ln(a)+ln(1+exp(ln(b)-ln(a)))
%             else
%                 llh_c1 = log(alpha*p1) + (data-1)*log(1-p1);
%                 llh_c2 = log((1-alpha)*p2) + (data-1)*log(1-p2);
%                 llh = sum(max(llh_c1,llh_c2) + log1p(exp(-abs(llh_c1-llh_c2))));% if a>b then ln(a+b)==ln(a)+ln(1+exp(ln(b)-ln(a)))
%             end
%         end
% new try --------------------------------
        function llh = compute_LLH_2comp(data,alpha,p1,p2,offset)
            % alpha, p1, p2 are scalars between 0 and 1
            % data is column vector of integers>=1
%             assert(nargin==5);
            assert(iscolumn(data))
            if nargin>4
                llh_c1 = log(alpha*p1) + (data-1)*log(1-p1);
                llh_c2 = log((1-alpha)*p2) + (data-1)*log(1-p2);
                % this is the probability that data<offset so = CDF
%                 CDF = GeometricMixtureModelMLE.compute_CDF_2comp(data,alpha,p1,p2,offset);
                CDF_m = (alpha*(1-(1-p1)^offset) + (1-alpha)*(1-(1-p2)^offset));
                % this adds in a correction factor (LLH/(1-CDF(offset))) to account for the fact that only tracks > offset
                % are included in the data
                llh = sum(max(llh_c1,llh_c2) + log1p(exp(-abs(llh_c1-llh_c2))))-sum(numel(data)*log(1-CDF_m));% if a>b then ln(a+b)==ln(a)+ln(1+exp(ln(b)-ln(a)))
            else
                llh_c1 = log(alpha*p1) + (data-1)*log(1-p1);
                llh_c2 = log((1-alpha)*p2) + (data-1)*log(1-p2);
                llh = sum(max(llh_c1,llh_c2) + log1p(exp(-abs(llh_c1-llh_c2))));% if a>b then ln(a+b)==ln(a)+ln(1+exp(ln(b)-ln(a)))
            end
        end
       % end new try -----------------------------------
        
        function CDF = compute_CDF_2comp(tracklength,alpha,p1,p2,offset) % -- need to figure out if this is actually needed anymore and fix
            % this calculates the CDF for a corrected tracklengh n, given offset and alpha,p1,p2
            % tracklength input is dataoffset corrected so that min tracklength =  1;
            assert(p1>=p2);
            if nargin>4 % there's an offset
                % the CDF value at the offset value - fraction of tracks that are cut off by offset
                cdf_ofs = (alpha*(1-(1-p1)^offset) + (1-alpha)*(1-(1-p2)^offset));
                % the true CDF at true tracklength values - where offset is added back so that the minimum tracklength is
                % offset + 1, ie: tracklength>offset
                originalCDF = alpha*(1-(1-p1).^(tracklength)) + (1-alpha)*(1-(1-p2).^(tracklength));
                % now shift and scale original CDF to correct for fraction removed by offset: divide by
                % 'complementary CDF' = prob that tracks are (strictly) longer than offset value
                CDF = (originalCDF - cdf_ofs)./(1-cdf_ofs);
            else
                CDF = alpha*(1-(1-p1).^tracklength) + (1-alpha)*(1-(1-p2).^tracklength);
            end
        end
        
        function CDF = compute_CDF_3comp(tracklength,a1,a2,a3,p1,p2,p3)
           CDF = a1*(1-(1-p1).^(tracklength)) + a2*(1-(1-p2).^(tracklength)) + ...
                            a3*(1-(1-p3).^(tracklength)); 
        end
        function llh = compute_LLH_3comp(data,alpha1,alpha2, lambda1,lambda2, lambda3,offset)
            % alpha, lambda1, lambda2 are scalars between 0 and 1
            % data is column vector of integers>=1
            assert(iscolumn(data))
            alpha3=max(1-alpha1-alpha2,0);
            llhM = zeros(numel(data),3);
            if nargin>6 % there's an offset
                llhM(:,1) = log(alpha1*lambda1) + (data+offset-1)*log(1-lambda1);
                llhM(:,2) = log(alpha2*lambda2) + (data+offset-1)*log(1-lambda2);
                llhM(:,3) = log(alpha3*lambda3) + (data+offset-1)*log(1-lambda3);
                llhM = sort(llhM,2,'descend');
                llh = llhM(:,1) + log1p(exp(llhM(:,2)-llhM(:,1)));
                CDF = GeometricMixtureModelMLE.compute_CDF_3comp(offset,alpha1,alpha2,alpha3,lambda1,lambda2,lambda3);
                llh = sum(llh + log1p(exp(llhM(:,3)-llh)))- sum(numel(data)*log(1-CDF));
            else
                llhM(:,1) = log(alpha1*lambda1) + (data-1)*log(1-lambda1);
                llhM(:,2) = log(alpha2*lambda2) + (data-1)*log(1-lambda2);
                llhM(:,3) = log(alpha3*lambda3) + (data-1)*log(1-lambda3);
                llhM = sort(llhM,2,'descend');
                llh = llhM(:,1) + log1p(exp(llhM(:,2)-llhM(:,1)));
                llh = sum(llh + log1p(exp(llhM(:,3)-llh))); 
            end
        end
        function out = convert_prob2lifetime(p,framerate,ci95p)
           lt = -framerate./log(1-p);   
           out = lt;
           if nargin>2
               ci95lt = abs(-framerate./log(1-(p+ci95p)) - lt);
               out = [lt; ci95lt];
           end
        end
        function p = convert_lifetime2prob(lt,framerate)
           p = 1-exp(-framerate./lt) ;
        end
        
        
        
    end % public static methods
    
end

