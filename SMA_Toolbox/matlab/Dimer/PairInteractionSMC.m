% File: PairInteractionSMC.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% 12/2016
% PairInteractionSMC 
%
% Runs a sequential monte carlo sampling (particle filter) for the pairinteraction model
%
% A high-speed C++ implementation

classdef PairInteractionSMC < IfaceMixin
    % PairInteractionSMC 
    %
    % Runs a sequential monte carlo sampling (particle filter) for the pairinteraction model
    %
    properties
        DefaultParams = struct(... %ErbanChapman with Diffusion
                            'DA',0.150,... %(um^2/s) particle A diffusion constant 
                            'DB',0.150,... %(um^2/s) particle B diffusion constant 
                            'DAB',0.035,... %(um^2/s) dimer diffusion constant 
                            'Dphi',1e1,... %(rad^2/s) dimer angle diffusion rate
                            'rho',0.050,... %(um) dimer bond distance between fluorophores
                            'sigma_rho',0.010,... %(um) dimer bond distance standard deviation
                            'gamma',1e2,... %(1/s) dimer bond distance relaxation rate
                            'rho_bind',0.050,... %(um) binding distance.  Binding is possible when distance < rho_bind
                            'rho_unbind',0.050,... %(um) unbinding distance.  
                            'lambda_bind',1e2,... % [1/s] Rate of binding when within rho_bind distance
                            'koff',1e0...%(1/s) (off) B->F rate.  (True kinetic rate)
                            );
        DefaultSimParams = struct(...
                                    'wmean',0.020, ... % [um] Mean localization error
                                    'wgamma',10., ... % [um] Mean localization error
                                    'wsigma',0.005 ... % [um] Mean localization error
                                   );
    end
   
    properties (Access = protected, Transient=true)
        initialized=false; %True once the object is correctly initialized by the initializeTrack method
    end

    methods
        function obj=PairInteractionSMC(varargin)
            % [in]
            %  parameters - struct
            %  prior - struct
            %  data - cell-array of observation matricies each [9
            %  Dmutual - Mutual diffusion constant (um^2/s)
            %  rho - Binding radius (um)
            %  lambda - Capture rate when within target (um^2/s)
            obj = obj@IfaceMixin(@PairInteractionSMC_Iface);
            if nargin>=1 && isstruct(varargin{1})
                params = varargin{1};
            else
                params = obj.DefaultParams;
            end
            if nargin==2
                data = varargin{2};
                obj.initialized = obj.openIface(params,data);
            else
                obj.initialized = obj.openIface(params);
            end
        end
        
        function Ndata = addData(obj,data)
            % Append new data to existing data
            % [in] data: cell-array of obersevation matricies
            % [out] Ndata: int - total number of observations in database
            Ndata = obj.call('addData',makecell(data));
        end

        function clearData(obj)
            % Clear all data
            obj.call('clearData');
        end

        function data = getData(obj)
            % Retrieve all data
            % [out] data: cell-array of obersevation matricies
            data = obj.call('getData');
        end
        
        function setParams(obj,params)
            % Set parameters
            % [in] params: struct of parameters to set
            obj.call('setParams',params);
        end

        function params = getParams(obj)
            % Get all parameters
            % [out] params: struct of parameters
            params = obj.call('getParams');
        end

        function setPrior(obj,prior)
            % [in] prior - a struct with the same format as PairInteractionMCMCModel.DefaultPrior
            if isfield(prior,'state')
                obj.call('setPriorState',prior.state);
            end
            if isfield(prior,'posShape')
                switch prior.posShape
                    case {'Gauss','Gaussian'}
                        obj.call('setPriorPositionGaussian',prior.posPoint, prior.posExtent);
                    case {'Rect','Rectangular'}
                        rect = repmat(prior.posPoint(:),2,1);
                        rect(3:4) = rect(3:4) + prior.posExtent;
                        obj.call('setPriorPositionRectangular',rect);
                    case {'Circ','Circular'}
                        obj.call('setPriorPositionCircular',prior.posPoint, prior.posExtent);                    
                end
            end
        end

        function success = runParticleFilter(obj,Nparticles, ProposalType)
            % Run the particle filter
            % [in] Nparticles- int - number of particles to simulates
            % [out] success - boolean - true if sampling was successful
            if nargin<3
                ProposalType='Transition';
            end
            success = obj.call('runParticleFilter',int32(Nparticles), ProposalType);
        end
        
        function obs_llh = obsLLH(obj)
            % Get the estimated LLH for each observation
            % [out] obsLLH - size:[Ndata] - log-likelihood of each observation 
            obs_llh = obj.call('obsLLH');
        end

        function [particles,llh] = sampleParticle(obj)
             % Sample a particle for each data and report LLH
             % [out] particles - cell-array size:[Ndata] each element is size:[5,Nobs] assigment to unknown variables 
             % [out] llh [optional] - size:[Ndata] log-likelihood of each particle 
             [particles,llh] = obj.call('sampleParticle');
        end

        function llh = sampleParticleLLH(obj)
             % Report SampleLLH of sampled particle for each data
             % [out] llh - size:[Ndata] log-likelihood of each sampled particle 
             llh = obj.call('sampleParticleLLH');
        end

        function [allParticles,weights,llh] = getAllParticles(obj,obsIdx)
            % Return all particles and the weights
            % [in] obsIdx - int [optional]  The index (1-based) of the observations to sample up-until
            %               [default = Nobs, i.e., to use all the data]  Set to a smaller value for
            %               debugging.
            % [out] allParticles - cell-array size:[Ndata] each element is [5,Nobs,Nparticles] 
            % [out] weights - size:[Nparticles, Ndata]
            % [out] llh - size:[Nparticles, Ndata]
            if nargin==2
                [allParticles,weights,llh] = obj.call('getAllParticles',int32(obsIdx)-1);
            else
                [allParticles,weights,llh] = obj.call('getAllParticles');
            end
        end

        function llh = computeLLH(obj,varargin)
            % Sample report LLH of a sampled particle for each data
            % [in] ObsData [optional] - cell-array of obersevation matricies [default:internal data]
            % [in] particle - cell-array of particle matricies
            % [out] llh - size:[Ndata] llh of particle for each data
            if nargin==2
                particles = varargin{1};
                if ~iscell(particles)
                    particles = {particles};
                end
                llh = obj.call('computeLLH',particles);
            else
                obsData = varargin{1};
                particles = varargin{2};

                if iscell(obsData) && iscell(particles)
                    llh = obj.call('computeLLH',obsData,particles);
                elseif ~iscell(obsData) && ~iscell(particles)
                    llh = obj.call('computeLLH',{obsData},{particles});
                    if ~isempty(llh)
                        llh = llh(1);
                    else
                        llh = NaN;
                    end
                end
            end
        end

        function llh = computeLLH_debug(obj,obsData,particles)
            % Sample report LLH of a sampled particle for each data
            % [in] ObsData - single obersevation matrix [9xNobs]
            % [in] particle - single particle matrix [5xNobs]
            % [out] llh - matrix [9xNobs] llh of particle for each variable
            llh = obj.call('computeLLH_debug',obsData,particles);
        end

        function  [obs,particles,llhAll,simParams] = simulate(obj,Nparticles,Nsteps,tmax,simType,simParams)
            % Use the DBN definition and the prior to simulate directly from the network
            % [in] Nparticles [uint64] - number of particles to simulate
            % [in] Nsteps [uint64]- number of steps to simulate
            % [in] tmax [double] - maximum simulation time.  Timesetep is interpolated based in this value
            % [in] simType [string] - 'simType': {"Ancestral","Blur","EC"} 
            % [in] simParams - struct of names doubles
            %                  'wmean':  [optional] mean value of error in observation (w) [Default: rho] 
            %                  'wgamma': [optional] relaxation rate of error in observation (w)  [Default: 10*dt]
            %                  'wsigma': [optional] standard deviation of error in observation [Default: .15*wmean]
            %                  'oversample': Required for "Blur".  Integer >1 number of additional samples to average for each measurment
            %                  'rho_unbind': Optional for "EC".  Unbinding radius.
            % [out] obs - cube size: [9,Nsteps,Nparticles] rows are [t, obs_xa, obs_xb, obs_ya, obs_xb, w_xa, w_xb, w_ya, w_yb]
            % [out] particles - cube size: [5,Nsteps,Nparticles] rows are [z, true_ax true_ay true_bx true_by]
            % [out] llhAll - matrix: [9,Nobs,Nparticles] llh of each hidden variable for each observation [z, yax yay ybx yby xax xay xbx xby]
            % [out] simParams - returned simulation parameters as a structure array of named doubles
            [obs,particles,llhAll,simParams]  = obj.call('simulate',uint32(Nparticles),uint32(Nsteps),tmax,simType,simParams);
        end
    end

    methods (Static=true)
        function [llhDelta, P] = testSimulation()
            smc = PairInteractionSMC();
            params = smc.DefaultSimParams;
            Nparticles = 10;
            Nsteps =100;
            tmax = 1;
            simType = 'Ancestral';
            [obs, particles,llhAll,simParams] = smc.simulate(Nparticles,Nsteps,tmax,simType,params);
            llhAllD = smc.computeLLH_debug(obs(:,:,1),particles(:,:,1));
            llhDelta = abs(llhAll(:,:,1) - llhAllD);
            llhDelta(llhDelta(:)<1e-7) = 0;
            llhDelta=max(llhDelta(:));
            P = particles(:,:,1);
        end
    end
end
