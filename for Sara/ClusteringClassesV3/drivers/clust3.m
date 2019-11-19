clear all
close all

XY = load('../data/threeD.csv');
E = 500;
minPts = 3;

c = Clustering();
c.Results = 'results';
%c.Plotting = true;
c.Alpha = 1.1;
shrinkFactor = 0.5;
algorithms = {'DBSCAN_Daszykowski', 'Hierarchical', 'Voronoi'};
for i = 1 : length(algorithms);
   algorithm = algorithms{i};
   if strcmp(algorithm, 'Voronoi')
      E = 0;
   end
   [nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
   fprintf('number of clusters = %d\n', nC);
   results = c.clusterStats(XY, C, centers, shrinkFactor)
   clusterFig = c.plotClusters3(XY, C, centers, ptsI, algorithm, shrinkFactor);
   showm(clusterFig);
end
