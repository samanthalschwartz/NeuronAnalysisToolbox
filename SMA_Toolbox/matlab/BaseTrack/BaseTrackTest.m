%% BaseTrackTest

% Function to test the new BaseTrack Class!
clear;

%% Set up Simulation Parameters
D = 0.2; % diffusion coefficient
ROI = [128 128]; % Region of Interest in pixels
N = 2; % number of particles diffusing
T = 100; % Number of frames
Kon = 0.1; % blink on rate
Koff = 0.1; % blink off rate
Error = 0.1; % average localization error

%% simulate trajectories

Xtrue = 2*D*randn(T,N); 
Ytrue = 2*D*randn(T,N);

vEx = Error*rand(T,N); % this is the std, make sure to square it before putting into SPT
vEy = Error*rand(T,N);

Xtrue(1,:) = ROI(1)/4 *(2 + randn(1,N));
Ytrue(1,:) = ROI(2)/4 *(2 + randn(1,N));

Xtrue = cumsum(Xtrue);
Ytrue = cumsum(Ytrue);

OX = Xtrue+vEx.*randn(T,N);
OY = Ytrue+vEy.*randn(T,N);

Time = (1:T)';
Time = repmat(Time,1,N);

%% Add junk values for other parameters

% array of values
% 'x', 'y', 'lambda', 'I', 'std_x', 'std_y', 'std_lambda', 'std_I', 'frame'

% vectorize matrix values
X = OX(:);
Y = OY(:);

var_x = vEx(:).^2;
var_y = vEy(:).^2;

frame = Time(:);

% populate relevant components of input matrix
Matrix_In = zeros(length(X),11);

Matrix_In(:,1) = X;
Matrix_In(:,2) = Y;
Matrix_In(:,6) = var_x;
Matrix_In(:,7) = var_y;
Matrix_In(:,11) = frame+5; % see if this messes things up

% Randomize Matrix_In to verify link associations work!
matRand = rand(length(X),1);
[~, sortIdx] = sort(matRand);

Matrix_In = Matrix_In(sortIdx,:);

% add intermittencies!!!
intRand = rand(length(X),1);
Matrix_In2 = Matrix_In;
Matrix_In2(intRand < 0.2,:) = [];

% create the input ROI
ROIin = [1 ROI(1) 1 ROI(2) min(Matrix_In(:,11)) max(Matrix_In(:,11))];

%% run tracking to see if it works

A = BaseTrack(Matrix_In,ROIin);
A.D =D;
A.doLAP;

%% run tracking with intermittency to see if it works

B = BaseTrack(Matrix_In2,ROIin);
B.D = D;
B.doLAP;


