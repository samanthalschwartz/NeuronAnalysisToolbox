function results_c = doClustering(n_ROIs, RoI, desc, results)
% Clustering for each label in each ROI. 

   E = 30;       % epsilon: max distance between 2 adjacent points in a cluster
   minPts = 5;   % minimum number of points in a cluster

   c = Clustering();
   c.Plotting = true;
   c.Alpha = 2;
   c.Valgorithm = 2;
   algorithm = 'DBSCAN';   % clustering algorithm
   shrinkFactor = 0.5;     % used to make cluster boundaries convex or concave
   options = 'O';

   XY = cell(2, 1);
   results_c = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      XY{1} = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      XY{2} = [ RoI{i}.X{2}, RoI{i}.Y{2} ];

      for j = 1 : 2
         [nC, C, centers, ptsI] = c.cluster(algorithm, XY{j}, E, minPts);
         fprintf('%s number of clusters ROI %d label %d = %d\n', ...
                 algorithm, i, j, nC);

         results_c{i}{j} = c.clusterStats(XY{j}, C, centers, shrinkFactor);
         results_c{i}{j}.C = C;
         results_c{i}{j}.centers = centers;
         results_c{i}{j}.ptsI = ptsI;

         txt = sprintf('%s ROI%d L%d', desc, i, j);
         clusterFig = c.plotClusters(XY{j}, C, centers, ptsI, txt, ...
                                     shrinkFactor, options);
         %showm(clusterFig);
         txt = sprintf('%s_ROI%d_L%d_%s', desc, i, j, algorithm);
         saveas(clusterFig, fullfile(results, sprintf('%s.png', txt)));
      end
   end

end
