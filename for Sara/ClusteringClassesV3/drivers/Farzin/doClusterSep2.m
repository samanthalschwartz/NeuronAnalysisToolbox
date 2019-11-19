function results_cs = doClusterSep2(n_ROIs, results_c)
% Find the nearest neighbor of each label 2 cluster to each label 1 cluster
% using center-to-center distances.

   results_cs = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      centers1 = results_c{i}{1}.centers';
      centers2 = results_c{i}{2}.centers';
      [indx, dist] = knnsearch(centers1, centers2);
      results_cs{i}.indx = indx;
      results_cs{i}.dist = dist;
   end

end
