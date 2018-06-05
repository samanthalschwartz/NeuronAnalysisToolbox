classdef MLE < handle
    % A generic maximum likelihood evaluation and plotting class for multi-parameter models
    %

    properties (Constant=true)
        MIN_LLH=-4;
    end

    properties

        data;
        llhFunc;  % A function that takes (param, data) as parameters and returns a LLH

        params=struct('N',[],...
                      'names',[],...
                      'units',[],...
                      'init_val',[],...
                      'lowerBound',[],...
                      'upperBound',[]);


        % This structure contains the settings for the optimizer
        optimizer=struct('solver','fmincon',...
                         'algorithm','sqp',...
                         'display','final-detailed',...
                         'options',[]...
                         );

        %This structure contains all of the optimization information
        mle=struct('valid',false,... %Boolean if we were successful
                   'param_est',[],... %Parameter vector at MLE estimate
                   'llh',[],... %Parameter estimate
                   'grad',[],... %Estimated gradient at MLE
                   'hess',[],... %Estimated hessian at MLE
                   'fisher_info',[],... %Fisher information at MLE
                   'opt_exit_flag',[],... %Exit flag for optimization
                   'opt_output',[]... %optimization output structure
                 );
        
    end

    
    
    methods
        function obj=MLE(llhFunc, data, param_init_val)
            obj.llhFunc = llhFunc;
            obj.data = data;
            obj.params.init_val = param_init_val;
            obj.params.N = numel(param_init_val);
            obj.params.names = arrayfun(@(i) sprintf('Param%i',i), 1:obj.params.N, 'Uniform',0);
        end

        function estimate(obj)
            prob = obj.generateProblem(obj.optimizer);
            [obj.mle.param_est, obj.mle.llh, obj.mle.hess, obj.mle.grad] = obj.runOptimizer(prob);
            obj.mle.valid=true; % Check for exit flag next
        end


        



        function viewProfileLikelihood(obj)
            N = obj.params.N;
            figure();
            for n=1:N
                subplot(ceil(N/2),2,n);
                fprintf('Plotting param: %i\n',n);
                obj.plotProfileLikelihood(n);
            end
        end




    end 
    methods (Access=protected)
        function r = profileLikelihoodRange(obj,param_idx, valid_rng)
            mle_val = obj.mle.param_est(param_idx);
            objective = @(v) obj.profileLikelihood(param_idx,v)-obj.mle.llh-obj.MIN_LLH;
            opts.Display='iter';
            opts.TolX = 1e-6;
            lb = fzero(objective, [valid_rng(1), mle_val-eps] , opts);
            ub = fzero(objective, [mle_val+eps, valid_rng(2)], opts);
            assert(ub-lb>eps);
            r=[lb,ub];
        end

        function [vals,llh,mle_val,mle_llh] = estimateProfileLikelihood(obj,param_idx,num_points, valid_rng)
            mle_val = obj.mle.param_est(param_idx);
            mle_llh = obj.mle.llh;
            rng = obj.profileLikelihoodRange(param_idx, valid_rng);
            vals = linspace(rng(1),rng(2),num_points);
            llh = obj.profileLikelihood(param_idx, vals) - mle_llh;
        end


        function profileLLH = profileLikelihood(obj, param_idx, param_vals)
            N = numel(param_vals);
            profileLLH = zeros(1,N);
            prob = obj.generateProblem(obj.optimizer);
            prob.x0(param_idx) = [];
            prob.ub(param_idx) = [];
            prob.lb(param_idx) = [];            
            for n=1:N
                prob.objective = obj.profileObjective(param_idx, param_vals(n));
                [~, profileLLH(n)] = obj.runOptimizer(prob);
            end
        end

        function plotProfileLikelihood(obj, param_idx)
            num_points=30;
            prob = obj.generateProblem(obj.optimizer);
            valid_rng= [prob.lb(param_idx), prob.ub(param_idx)];
            [vals,llh,mle_val,mle_llh] = obj.estimateProfileLikelihood(param_idx,num_points,valid_rng);
            plot(vals,llh,'k-','DisplayName','Normalized Profile Likelihood');
            hold();
            s=sprintf('MLE Val:%.3g  (TrueLLH:%.3g)',mle_val,mle_llh);
            plot(mle_val, 0,'r*','DisplayName',s);
            ylim([obj.MIN_LLH,0]);
            xlim([vals(1),vals(end)]);
            legend('Location','best');
            xlabel(sprintf('%s',obj.params.names{param_idx}));
            ylabel('Normalized likelihood');
        end
        
        function opts = generateOpts(obj,Ostruct)
            opts = optimset(Ostruct.solver);
            opts.Algorithm = Ostruct.algorithm;
            opts.Diagonstics = 'off';
            opts.Display = 'off';
%             opts.Display = Ostruct.display;
            opts.HessianApproximation = 'bfgs';
        end

        function prob = generateProblem(obj,Ostruct)
            opts = obj.generateOpts(Ostruct);
            prob.objective = @(p) -obj.llhFunc(p,obj.data);
            prob.x0 = obj.params.init_val;
            prob.lb = obj.params.lowerBound;
            prob.ub = obj.params.upperBound;
            prob.options = opts;
            prob.solver = Ostruct.solver;
        end
        
        function F = profileObjective(obj, param_idx, param_val)
            F = @(p) -obj.llhFunc( [p(1:param_idx-1), param_val, p(param_idx:end)],obj.data);
        end

       

        function [mle,llh,hess, grad] = runOptimizer(obj, prob)
            optF=str2func(obj.optimizer.solver);
            [mle, fval, exit_flag, out, ~, grad, hess] = optF(prob);
            llh=-fval;
            grad=-grad;
            hess=-hess;
%             if exit_flag~=1
%                 fprintf('Exit Flag: %i\n',exit_flag);
%             end
        end
    end
end

