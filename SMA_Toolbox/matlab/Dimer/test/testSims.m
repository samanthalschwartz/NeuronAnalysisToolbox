
Nsteps=1e3;
simTimeStep = 1e-3;
sampleTimeStep = 1e-2;
Ntraj=10;

DA = 0.050;  %um^2/s
DB = 0.050;  %um^2/s
DAB = 0.020; %um^2/s
Dtheta = 10*2*pi; % rad/s
rho = 0.050; %um
rho_sigma = 0.015; %um
gamma = 1e4; %[1/s] 

S=DimerSimulator(Nsteps,simTimeStep);
S.DA = DA;
S.DB = DB;
S.DAB = DAB;
S.dimerRigidLength = rho;
S.dimerLengthStddev = rho_sigma;
S.dimerRelaxationRate = gamma;
S.locErrorMean = 0.015; %[um]
S.overlayErrorStddev = 0.005; % micron

Ts=S.simulateTrajectories(Ntraj);
Os=S.observeTrajectories(Ts,sampleTimeStep);


S.browseObservedTrajectories(Ts,Os);
S.plotDistances(Ts{1}, Os{1})
