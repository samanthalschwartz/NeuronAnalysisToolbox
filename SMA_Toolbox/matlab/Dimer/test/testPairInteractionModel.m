
%% Initialize simulation
Nsteps = 10000;
deltaT = 1e-3;
S= DimerSimulator(Nsteps, deltaT);
S.interactionType = 'Fixed';
S.interactionFixedTime = 1;
S.locErrorMean = 0.030;
S.dimerOrientation = 'DriftingRigid';
S.dimerOrientationDriftRate = 10;
Ts = S.simulateTrajectories(1);
S.browseTrajectories(Ts);
frameT = 1e-2;
Obs = S.observeTrajectories(Ts, frameT);
S.plotDistances(Ts{1},Obs{1});

%S.browseTrajectories(Ts);
%S.browseObservedTrajectories(Ts,Obs);
% S.plotDMLE(Ts,0);
S.plotDistances(Ts{1}, Obs{1});

%% Model
M = PairInteractionModel(FObs);
kon=1;
koff=5;
M.initializeModelParams('Distance', [S.DA,kon, koff]);
M.initializeModelParams('Diffusive', [S.DA, S.DB, S.DAB,kon, koff]);
M.calculateModel('Distance');
M.plotHMMObservation(1);

% 
% PairInteractionModel.plotHMM_Distance(trueState, pair, prior, params);
% PairInteractionModel.plotHMM_Sample(nSamples, trueState, pair, prior, params);

