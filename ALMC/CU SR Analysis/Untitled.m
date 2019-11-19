Sigma_Reg = [10, 10];
SRc = SRcluster();
%SRc.PvalueStatistics = true;
%SRc.PlotFigures = false;
% clusterSR can work in nD.
[xy_SR, sigma_SR, combined] = ...
    SRc.clusterSR([Xs_ch1, Ys_ch1], [prec_ch1, prec_ch1], Sigma_Reg);
% Functions below assume 2D data.
% cutoffFigs = SRc.cutoffPlots();
SRcollapseFig = SRc.plotSRcollapse();
SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();