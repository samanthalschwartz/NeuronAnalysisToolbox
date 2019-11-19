clear all
close all

pixel2nm = 16000/150;

load('../data/Fit4.mat');

%SRtest.DriftCorrect = 1;
%SRtest.FrameConnect = 1;
%SRtest.Threshold;
Sigma_Reg = std(SRtest.DriftCorrect_XYShift) .* pixel2nm;  %nm
if isnan(Sigma_Reg)
   Sigma_Reg = [10, 10];
end

X = double(SRtest.Results.X) .* pixel2nm;   % nm
Y = double(SRtest.Results.Y) .* pixel2nm;   % nm
X_STD = double(SRtest.Results.X_STD) .* pixel2nm;   % nm
Y_STD = double(SRtest.Results.Y_STD) .* pixel2nm;   % nm

SRc = SRcluster();
SRc.PvalueStatistics = true;
%SRc.PlotFigures = false;
SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
%cutoffFigs = SRc.cutoffPlots();
SRc.E = SRc.chooseCutoff();
%SRc.Algorithm = 'DBSCAN_Daszykowski';   SRc.minPts = 3;   SRc.E = 30;
SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();
%saveas(SRclusterFig, 'A', 'fig');
%saveas(analysisFigs{1}, 'B', 'fig');
%close all
