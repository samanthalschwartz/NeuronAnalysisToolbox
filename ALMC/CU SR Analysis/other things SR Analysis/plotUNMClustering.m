function [c1,c2,f1,f2] = plotUNMClustering(Xs_ch1,Ys_ch1,Xs_ch2,Ys_ch2,input,savedir)
if nargin<2
    saveflag = false;
else
    saveflag = true;
end
if nargin>0
    E = input.E;
    minPts = input.minPts;
    algorithm = input.algorithm;
    savestr = input.savestr;
    
else
    E = 20;   % nm
    minPts = 10;
    algorithm = 'DBSCAN_Tran';
    savestr  = '';
end
%                DBSCAN_Daszykowski or DBSCAN
%                DBSCAN_Daszykowski_noE (computes its own value for E)
%                DBSCAN_Kovesi
%                DBSCAN_Pehlke
%                DBSCAN_Tran
%                Getis      (Carolyn Pehlke's Getis statistic based algorithm)
%                Hierarchal (Matlab's hierarchal clustering algorithm)
%                Voronoi    (Florian Levet et al's Voronoi based algorithm)
% channel 1
XY = [Xs_ch1,-Ys_ch1];
c1 = Clustering();
[nc1, C1, c1enters, ptsI] = c1.cluster(algorithm, XY, E, minPts);
fprintf('number of clusters = %d\n', nc1);
% results = c1.clusterStats(XY, c1, c1enters);
c1lusterFig = c1.plotClusters(XY, C1, c1enters, ptsI, algorithm);
f1 = figure(c1lusterFig);

% channel 2
XY = [Xs_ch2,-Ys_ch2];
c2 = Clustering();
[nc2, C2, c2enters, ptsI] = c2.cluster(algorithm, XY, E, minPts);
fprintf('number of clusters = %d\n', nc2);
% results = c2.clusterStats(XY, c2, c2enters);
c2lusterFig = c2.plotClusters(XY, C2, c2enters, ptsI, algorithm);
f2 = figure(c2lusterFig);

%-- saving if plotflag---
if saveflag
    save(fullfile(savedir,['cluster_ch1_' savestr '_E=' num2str(E) '_minPts=' num2str(minPts) '_' algorithm]),'c1'); 
    saveas(f1,fullfile(savedir,['cluster_ch1_' savestr '_E=' num2str(E) '_minPts=' num2str(minPts) '_' algorithm]));
    saveas(f1,fullfile(savedir,['cluster_ch1_' savestr '_E=' num2str(E) '_minPts=' num2str(minPts) '_' algorithm]),'png');
    
    save(fullfile(savedir,['cluster_ch2_' savestr '_E=' num2str(E) '_minPts=' num2str(minPts) '_' algorithm]),'c2'); 
    saveas(f2,fullfile(savedir,['cluster_ch2_' savestr '_E=' num2str(E) '_minPts=' num2str(minPts) '_' algorithm])); 
    saveas(f2,fullfile(savedir,['cluster_ch2_' savestr '_E=' num2str(E) '_minPts=' num2str(minPts) '_' algorithm]),'png'); 
end
end