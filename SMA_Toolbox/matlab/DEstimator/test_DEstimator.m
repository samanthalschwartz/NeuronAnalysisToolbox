% test_DEstimator
% Test Script to Show Case the DEstimator class
%
% Peter K. Relich (physx.grad@gmail.com)
% Mark J. Olah (mjo@cs.unm.edu)

est = DEstimator(); % create empty instance of the DEstimator Object

% Simulation parameters
Dsim = 0.05; % Simulated diffusion constant (um^2/s)
N = 50; % Number of observations (trajectory length)
SE = 0.000; %  Mean standard error of observations (um). Simulated variances will follow a gamma distribution with this mean.
dT = 0.05; % Frame-to-frame time (s) (typically this is the same as exposure time)

% generate a 1D simulated trajectory with the above parameters
[ObSim, Tvec, SEvec] = DEstimator.simulate1D(Dsim,N,SE,dT,dT,'gamma');

% Initialize the simulated trajectory into the new 'est' DEstimator object
success = est.initializeTrack(ObSim, Tvec, SEvec, dT,'recursive');

% Find the maximum likelihood estimate for the diffusion constant of the trajectory
[D_mle, D_mle_llh] = est.MLE();

% sample some diffusion coefficients
min_D = 1e-4;
max_D = 3e-1;
num_samples = 100;
D_vals = linspace(min_D,max_D,num_samples);

% Look at the likelihood distribution for the 1-D trajectory
llh_vec = est.LLH(D_vals);

% convert the log likelihood into a re-scaled likelihood
lh_vec = exp(llh_vec - max(llh_vec));

%Plot the results.
plot(D_vals,lh_vec,'r-') %plot likelihood curve
hold('on');
plot(D_mle, exp(D_mle_llh-max(llh_vec)),'bs','MarkerFaceColor','k') %plot D_mle value location
legend('Relative Likelihood of D-value',sprintf('Maximum likelihood estimate: %f (um^2/s)',D_mle));
title('Diffusion constant likelihood');
xlabel('Diffusion constant (um^2/s)');
ylabel('Relative Likelihood');