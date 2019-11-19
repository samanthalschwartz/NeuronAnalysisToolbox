function doBiStats(n_ROIs, RoI, desc, results)
% Pairwise mutual distances and bivariate Ripley's statistics for each ROI.

   c = Clustering();
   c.Fig_ext = 'png';
   c.Results = results;
   particle_types = {'L1', 'L2'};

   results_birip = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      txt = sprintf('%s_ROI%d', desc, i);
      XY1 = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      XY2 = [ RoI{i}.X{2}, RoI{i}.Y{2} ];

      P{1} = XY1;
      P{2} = XY2;
      H_nm = RoI{i}.ROI(2) - RoI{i}.ROI(1);
      V_nm = RoI{i}.ROI(4) - RoI{i}.ROI(3);
      txt = sprintf('%s_ROI%d', desc, i);
      c.cluster_stats('PairwiseMutualDist', ...
                      P, txt, particle_types, H_nm, V_nm);
      c.cluster_stats('BivariateRipley', P, txt, particle_types, H_nm, V_nm);
   end

end
