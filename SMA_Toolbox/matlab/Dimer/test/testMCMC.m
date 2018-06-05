
deltaT_sim = 1e-3;
nSteps = 1e3;
DA = 0.1;
DB = 0.15;
DAB = 0.075;
rho = 0.05; % [um]

NTs = 2; % Number of trajectories to simulate
frameT = 0.05; %[s] frame Time of observation


%% Simulate Trajectory Pair Observations with DimerSimulator

S = DimerSimulator(nSteps, deltaT_sim);
S.DA = DA;
S.DB = DB;
S.DAB = DAB;
S.dimerRigidLength = rho;

%Simulate trajectory pairs
Ts = S.simulateTrajectories(NTs);

%Observe trajectory pairs
Os = S.observeTrajectories(Ts, frameT);


%% Create a model

M = PairInteractionMCMCModel(frameT,Os);
M.initializeModelParams('Diffusive');

%% Evaluate log-likelihood of actual data
[LLH, llh_all_vars ] = M.computeLLH(M.particlesTrue);
M.plotLLHAllVars(llh_all_vars);

P_all_bound = M.particlesTrue;
P_all_free = M.particlesTrue;
for n=1:M.Ndata
    P_all_bound{n}(:,1)=1;
    P_all_free{n}(:,1)=0;
end

[LLH_bound, llh_all_vars_bound ] = M.computeLLH(P_all_bound);
M.plotLLHAllVars(llh_all_vars_bound);


[LLH_free, llh_all_vars_free ] = M.computeLLH(P_all_free);
M.plotLLHAllVars(llh_all_vars_free);
