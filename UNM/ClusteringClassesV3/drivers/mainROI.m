clear all
close all

pixel2nm = 16000/150;

%load('../data/Fit4_ROI.mat');
load('../data/Fit3_ROIXAlb1.mat');

%SRtest.DriftCorrect = 1;
%SRtest.FrameConnect = 1;
%SRtest.Threshold;
%Sigma_Reg = std(SRD.DriftCorrect_XYShift) .* pixel2nm;  %nm
%if isnan(Sigma_Reg)
%   Sigma_Reg = [10, 10];
%end

X = double(RoI{1}.X);   % nm
Y = double(RoI{1}.Y);   % nm
X_STD = double(RoI{1}.X_STD);   % nm
Y_STD = double(RoI{1}.Y_STD);   % nm

SRc = SRcluster();
SRc.Cutoff = 10 : 10 : 1000;
SRc.LoS = 0.01;
SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
%cutoffFigs = SRc.cutoffPlots();
SRc.E = SRc.chooseCutoff();
SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();
