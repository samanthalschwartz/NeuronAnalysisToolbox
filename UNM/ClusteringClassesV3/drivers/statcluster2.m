clear all
close all

XY = load('../../analysis/methods/data/pts.csv');
E = 30;
minPts = 1;

c = Clustering();
algorithm = 'Hierarchical';
[nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
fprintf('number of clusters = %d\n', nC);
results = c.clusterStats(XY, C, centers)
clusterFig = c.plotClusters(XY, C, centers, ptsI, algorithm);
%showm(clusterFig);

%Cxy = cell(1, nC);
%for i = 1 : nC
%   Cxy{i} = XY(C{i}, :);
%end

[x, y] = textread('../data/9021_5.txt',  '%*u %u %u %*u', 'headerlines', 1);
XY_5 =  [x, y];
[x, y] = textread('../data/9021_10.txt', '%*u %u %u %*u', 'headerlines', 1);
XY_10 = [x, y];
c.Results = 'results';
P{1} = XY_5  ./ 2.7559;
P{2} = XY_10 ./ 2.7559;
H_nm = 7400 / 2.7559;
V_nm = 6000 / 2.7559;
base_name = '9021';
particle_types = {'5', '10'};
c.cluster_stats('PairwiseDist',    P, base_name, particle_types, H_nm, V_nm);
c.cluster_stats('PairwiseMutualDist', ...
                                   P, base_name, particle_types, H_nm, V_nm);
c.cluster_stats('Hopkins',         P, base_name, particle_types, H_nm, V_nm);
c.cluster_stats('Ripley',          P, base_name, particle_types, H_nm, V_nm);
c.cluster_stats('BivariateRipley', P, base_name, particle_types, H_nm, V_nm);
c.cluster_stats('Dendrogram',      P, base_name, particle_types, H_nm, V_nm);
