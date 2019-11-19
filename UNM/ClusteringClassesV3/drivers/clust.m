clear all
close all

XY = load('../../analysis/methods/data/pts.csv');
E = 30;
minPts = 1;

c = Clustering();
%c.Plotting = true;
c.Alpha = 1.1;
shrinkFactor = 0.5;
algorithms = {'DBSCAN_Daszykowski', 'Getis', 'Hierarchical', 'Voronoi'};
for i = 1 : length(algorithms);
   algorithm = algorithms{i};
   [nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
   fprintf('number of clusters = %d\n', nC);
   results = c.clusterStats(XY, C, centers, shrinkFactor)
   clusterFig = c.plotClusters(XY, C, centers, ptsI, algorithm, shrinkFactor);
   showm(clusterFig);
end
