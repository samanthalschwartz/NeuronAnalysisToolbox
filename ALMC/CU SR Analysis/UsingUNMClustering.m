% Example main program:
% 
%16000/150;
% 
%    load(...);
% 
%    Sigma_Reg = std(SRtest.DriftCorrect_XYShift) .* pixel2nm;  % nm
%    if isnan(Sigma_Reg)
%       Sigma_Reg = [10, 10];
%    end
% 
Sigma_Reg = [10, 10];
X = data.SML_data.ch1.position_x(1:1000);
Y = data.SML_data.ch1.position_y(1:1000);
X_STD = data.SML_data.ch1.precision;
Y_STD = data.SML_data.ch1.precision;
%    X = double(SRtest.Results.X) .* pixel2nm;   % nm
%    Y = double(SRtest.Results.Y) .* pixel2nm;   % nm
%    X_STD = double(SRtest.Results.X_STD) .* pixel2nm;   % nm
%    Y_STD = double(SRtest.Results.Y_STD) .* pixel2nm;   % nm
%
SRc = SRcluster();
%SRc.PvalueStatistics = true;
%SRc.PlotFigures = false;
% clusterSR can work in nD.
[xy_SR, sigma_SR, combined] = ...
    SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
% Functions below assume 2D data.
% cutoffFigs = SRc.cutoffPlots();
SRcollapseFig = SRc.plotSRcollapse();
SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();
%%
%  XY = load('../../analysis/methods/data/pts.csv');  % nm
% XY = SRc.XY;%[X,Y];
% XY = [X,Y];
% XY = [data.SML_data.ch1.position_x, data.SML_data.ch1.position_y];
%---for ch1 ----
XY = [Xs_ch1,-Ys_ch1];
%-- for ch2----
XY = [Xs_ch2,-Ys_ch2];
%----
E = 20;   % nm
minPts = 10;

c = Clustering();
%    algorithm = 'DBSCAN_Daszykowski';
algorithm = 'DBSCAN_Tran';
%               DBSCAN_Daszykowski or DBSCAN
%                DBSCAN_Daszykowski_noE (computes its own value for E)
%                DBSCAN_Kovesi
%                DBSCAN_Pehlke
%                DBSCAN_Tran
%                Getis      (Carolyn Pehlke's Getis statistic based algorithm)
%                Hierarchal (Matlab's hierarchal clustering algorithm)
%                Voronoi    (Florian Levet et al's Voronoi based algorithm)
[nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
fprintf('number of clusters = %d\n', nC);
results = c.clusterStats(XY, C, centers)
clusterFig = c.plotClusters(XY, C, centers, ptsI, algorithm);
%showm(clusterFig);
figure(clusterFig);

%    [x, y] = textread('../data/9021_5.txt',  '%*u %u %u %*u', ...
%                                             'headerlines', 1);
%    XY_5 =  [x, y];
%    [x, y] = textread('../data/9021_10.txt', '%*u %u %u %*u', ...
%                                             'headerlines', 1);
%    XY_10 = [x, y];
%    c.Results = 'results';
%    P{1} = XY_5  ./ 2.7559;
%    P{2} = XY_10 ./ 2.7559;
%    H_nm = 7400 / 2.7559;
%    V_nm = 6000 / 2.7559;
%    base_name = '9021';
%    particle_types = {'5', '10'};
%    c.cluster_stats('PairwiseDist',    P, base_name, particle_types, ...
%                                       H_nm, V_nm);
%    c.cluster_stats('Hopkins',         P, base_name, particle_types, ...
%                                       H_nm, V_nm);
%    c.cluster_stats('Ripley',          P, base_name, particle_types, ...
%                                       H_nm, V_nm);
%    c.cluster_stats('BivariateRipley', P, base_name, particle_types, ...
%                                       H_nm, V_nm);
%    c.cluster_stats('Dendrogram',      P, base_name, particle_types, ...
%                                       H_nm, V_nm);
