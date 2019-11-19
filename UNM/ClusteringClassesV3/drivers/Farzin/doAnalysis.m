function [results_pcc, resultsRC_pcc, results_c, results_cs] = ...
   doAnalysis(data_dir, results_dir, desc, pixel2nm)

   [X1, Y1, X2, Y2] = getData(data_dir);
   % double is needed below as the boundary command used with clustering only
   % works correctly with doubles.
   XY1 = double([X1, Y1] .* pixel2nm);
   XY2 = double([X2, Y2] .* pixel2nm);

   load(fullfile(results_dir, 'ROIs.mat'));

   doBiStats(n_ROIs, RoI, desc, results_dir);
   [results_pcc, resultsRC_pcc] = doPairCorr(n_ROIs, RoI, desc, results_dir);
   results_c  = doClustering(n_ROIs, RoI, desc, results_dir);
   results_cs = doClusterSep2(n_ROIs, results_c);

   save(fullfile(results_dir, sprintf('%s_results.mat', desc)), 'n_ROIs', ...
        'RoI', 'results_pcc', 'resultsRC_pcc', 'results_c', 'results_cs');

end
