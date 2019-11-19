clear all
close all

RT = ROITools();
c = Clustering();

RT.ROI_sizes = [2000, 2000];
E = 30;
minPts = 3;
options = 'O';
shrinkFactor = 0.5;
c.Plotting = true;

infile = '../data/Fit4.mat';

[n_ROIs, RoI, Sigma_Reg] = RT.getROI({infile});

basefile = regexprep(infile, '.mat', '');

%algorithms = {'DBSCAN_Daszykowski', 'Getis', 'Hierarchical', 'Voronoi'};
algorithms = {'Hierarchical'};
for j = 1 : n_ROIs
   fprintf('ROI %d: %d points\n', j, numel(RoI{j}.X{1}));
   for i = 1 : length(algorithms)
      X = RoI{j}.X{1};
      Y = RoI{j}.Y{1};
      XY = [X, Y];

      algorithm = algorithms{i};
      [nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
      fprintf('%s number of clusters = %d\n', algorithm, nC);
      results = c.clusterStats(XY, C, centers, shrinkFactor)
%     save(sprintf('%s_%s.mat', basefile, algorithm), 'results', 'c');
      clusterFig = c.plotClusters(XY, C, centers, ptsI, algorithm, ...
                                  shrinkFactor, options);
      showm(clusterFig);
%     saveas(clusterFig, sprintf('%s_%s.fig', basefile, algorithm));
   end
end
