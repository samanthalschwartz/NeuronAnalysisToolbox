% MLE_X_Comp

% Script to test the validity of returning the maximum likelihood positions
% with the Laplace method.
clear;
close all; % clean up the workspace

% First set up the simulation parameters
Dsim = 1; % Simulated diffusion parameter
dT = 1; % Simulated trajectory frame time
N = 100; % Numer of N observations per track
meanSE = 0.75; % The mean error term, gamma variance
exposureT = dT; % Make the exposure time 0 for now

% Generate a 1D simulated trajectory of true positions
X = DEstimator.simulate1DDiffusion(N+1, Dsim, dT); %make N+1 X values
% Convert true positions into observations
var_gen_method = 'gamma'; % localization variance

% Need to know mean realizations
SEvec = DEstimator.generateObsStandardError(meanSE,N,var_gen_method);

Tvec = (0:N-1)'.*dT;
alpha = exposureT/(2*dT);
Xbar_mean = (1-alpha)*X(1:end-1) + alpha*X(2:end);
Xbar_var = 2*Dsim*exposureT*(1/3-alpha/2);
Xbar_realize = Xbar_mean+randn(N,1).*sqrt(Xbar_var); % we want to know the realizations
ObSim = Xbar_realize + randn(N,1).*SEvec;

% Create an instance of the DEstimator class
est = DEstimator(ObSim, Tvec, SEvec, exposureT);

% Perform the MLE estimation of the true positions
[posMLE, posSE, MLE_LLH_X] = est.MLEPositions(Dsim);

% Plot overlay observations, true positions, estimated positions with standard error
figure;
h = gcf;
hold on;
plot(Tvec,ObSim,'bo','MarkerSize',6,'LineWidth',1);
plot(Tvec,Xbar_realize,'kx','MarkerSize',6,'LineWidth',2);
errorbar(Tvec,posMLE,posSE,'r+','MarkerSize',6,'LineWidth',1);
xlim([Tvec(1) Tvec(end)])
ha = gca;
set(ha,'FontSize',14);
legend({'Observations','True Motion Averaged Positions','Estimated True Positions'});
title('MLE Positions vs. True Positions','FontSize',14);
xlabel('Frame Number','FontSize',14);
ylabel('Position','FontSize',14);
hold off;

meansqError = mean((posMLE-Xbar_realize).^2);
obssqError = mean((ObSim - Xbar_realize).^2);

m1 = sprintf('The mean squared error of the MLE positions is %d',meansqError);
m2 = sprintf('The mean squared error of the observed positions is %d',obssqError);

% show the mean squared errors of the MLE and observations
disp(m1);
disp(m2);