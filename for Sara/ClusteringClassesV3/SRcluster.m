classdef SRcluster < handle

% SRcluster class written by Michael Wester and Keith Lidke (7/21/2018)
%    <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2018 by Michael J. Wester and Keith A. Lidke
%
% Example main program:
%
%    pixel2nm = 16000/150;
%
%    load(...);
%
%    Sigma_Reg = std(SRtest.DriftCorrect_XYShift) .* pixel2nm;  % nm
%    if isnan(Sigma_Reg)
%       Sigma_Reg = [10, 10];
%    end
%
%    X = double(SRtest.Results.X) .* pixel2nm;   % nm
%    Y = double(SRtest.Results.Y) .* pixel2nm;   % nm
%    X_STD = double(SRtest.Results.X_STD) .* pixel2nm;   % nm
%    Y_STD = double(SRtest.Results.Y_STD) .* pixel2nm;   % nm
%
%    SRc = SRcluster();
%    %SRc.PvalueStatistics = true;
%    %SRc.PlotFigures = false;
%    % clusterSR can work in nD.
%    [xy_SR, sigma_SR, combined] = ...
%       SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
%    % Functions below assume 2D data.
%    cutoffFigs = SRc.cutoffPlots();
%    SRcollapseFig = SRc.plotSRcollapse();
%    SRclusterFig = SRc.plotSRclusters();
%    [results, analysisFigs] = SRc.analyzeSRclusters();

% =============================================================================
properties
% =============================================================================

   Method = 'hierarchal_singlelabel';  % H-SET collapse method
   %Method = 'keith_HC';
   %Method = 'trivial';
   Algorithm = 'hierarchalSimple';     minPts = 1;  % default clustering setup
   %Algorithm = 'DBSCAN_Daszykowski';   minPts = 3;
   Dim = 2;   % allowed values are 2 (2D) and 3 (3D)
   PlotFigures = true;         % plot various cluster related figures
   Printing = true;            % print various cluster statistics
   ProduceLegend = true;       % produce legends for some cluster figures
   % If true, produce various Pvalue plots at the expense of slowing down the
   % collapsing algorithm
   PvalueStatistics = false;
   Timing = true;              % produce timings for clustering
   Verbose = false;            % print cluster sizes as found

   % range of cutoff distances (nm) for analyses with cutoff as a parameter
   Cutoff = 10 : 10 : 5000;
   E = 30;       % cutoff distance (nm)
   LoS = 0.01;   % level of significance
   A_ROI = 1;    % area of the ROI (nm^2)
   ShrinkFactor = 0.5;   % boundary shrink factor (0 = convex hull)

   % experimental collapse parameters
   MaxLD = inf;   % maximum linkage distance allowed
   MaxLN = inf;   % maximum leaf nodes allowed

   % reject properties
   Histograms = false;           % reject histograms
   DetailedStatistics = false;   % reject detailed statistics

   % Voronoi properties
   Alpha      = 2;       % Ratio of local density / overall density for a
                         % point's Voronoi region to be considered sufficiently
                         % dense for clustering purposes.
   Valgorithm = 2;       % Voronoi algorithm to apply:
                         % 1   [0] calculations consider Voronoi regions only
                         % 2   [1] calculations consider Voronoi regions and
                         %         their adjacent neighbors
                         % 3   [1M] consider the median density of each cell
                         %          and its neighbors
   Plotting   = false;   % Produce Voronoi plots.
   PtIDs      = false;   % Label the points in the plots.

% =============================================================================
end

properties(SetAccess = protected)
% =============================================================================

   Sigma_Reg;        % sigma registration
   XY_orig;          % input (x, y)
   Sigma_orig;       % input sigma
   C_orig;           % initial clusters at cutoff E
   XY;               % collapsed (x, y)
   Sigma;            % collapsed sigma
   C;                % clusters at cutoff E of collapsed data
   Nodes_combined;   % multiple nodes collapsed into single nodes
   PDfig = [];       % Pvalue vs (linkage distance, # leaf nodes) figures
   PD = [];          % Pvalue vs (linkage distance, # leaf nodes) data

% =============================================================================
end % properties(SetAccess = protected)

methods
% =============================================================================

function [XY, sigma, combined] = clusterSR(obj, XY_orig, sigma_orig, Sigma_Reg)
% Combine multiple clustered points into single localizations when appropriate
% via a top-down descent through a hierarchal dendrogram relationship between
% points.
% n is the original number of points and m is the dimension.
% n' is the final number of points after combinations have occurred.
%
% Inputs:
%    XY_orig      n x m matrix of coordinates
%    sigma_orig   n x m matrix of position uncertainties (1 standard deviation)
%    Sigma_Reg    1 x m array of registration error (1 standard deviation)
%
% Outputs:
%    XY           n' x m final coordinate matrix
%    sigma        n' x m final position uncertainty matrix
%    combined     cell array of indices of combined points per cluster

   obj.XY_orig = XY_orig;
   obj.Sigma_orig = sigma_orig;
   obj.Sigma_Reg = Sigma_Reg;

   % Find clusters in the data.
   if obj.Timing
      tic
   end
   switch obj.Method
   case 'trivial'
      % For testing purposes only.
      obj.XY = XY_orig;
      obj.Sigma = sigma_orig;

   case 'hierarchal_singlelabel'
      % Collapse multiple emitters into single emitters.
      if obj.PvalueStatistics
         [obj.XY, obj.Sigma, obj.Nodes_combined] =          ...
            obj.hierarchalSingleLabelP(XY_orig, sigma_orig, ...
                                       Sigma_Reg, obj.LoS);
      else
         [obj.XY, obj.Sigma, obj.Nodes_combined] =               ...
            SRcluster.hierarchalSingleLabel(XY_orig, sigma_orig, ...
                                            Sigma_Reg, obj.LoS);
      end

   case 'keith_HC'
      % Collapse multiple emitters into single emitters.
      [obj.XY, obj.Sigma, obj.Nodes_combined] = ...
         SRcluster.keith_HC(XY_orig, sigma_orig, Sigma_Reg, obj.LoS);

   otherwise
      error('Unknown method: %s\n', obj.Method);
   end
   if obj.Timing
      toc
   end

   XY = obj.XY;
   sigma = obj.Sigma;
   combined = obj.Nodes_combined;

end

% =============================================================================

function set_XY_sigma(obj, Sigma_Reg, XY_orig, sigma_orig, XY, sigma, combined)
% Copy inputs to the SRcluster object.

   obj.Sigma_Reg = Sigma_Reg;
   obj.XY_orig = XY_orig;
   obj.Sigma_orig = sigma_orig;
   obj.XY = XY;
   obj.Sigma = sigma;
   obj.Nodes_combined = combined;

end

% =============================================================================

function [cutoffFigs, M] = cutoffPlots(obj)
% Produce a series of plots by varying the cutoff distance in a hierarchal
% clustering algorithm.
%
% Outputs:
%    cutoffFigs   figure handles for the following plots:
%                 # of clusters vs. cutoff
%                 # of objects per clusters vs. cutoff
%                 mean cluster area vs. cutoff
%                 mean nearest neighbor center-to-center distances vs. cutoff
%    M            # of clusters with > 2 objects per cutoff distance
%                 [1 x Ncutoff]

   if size(obj.XY, 1) <= 1
      warning('cutoffPlots: 0 or 1 points provided!');
      return
   end

   Z = linkage(obj.XY, 'single');

   n_cutoffs = length(obj.Cutoff);
   N = zeros(1, n_cutoffs);    % # of clusters
   M = zeros(1, n_cutoffs);    % # of clusters with > 2 objects
   A = zeros(1, n_cutoffs);    % mean area of clusters with > 2 objects
   Oa = zeros(1, n_cutoffs);   % mean # of objects per cluster
   Om = zeros(1, n_cutoffs);   % mean # of objects per cluster with > 2 objects
   for i = 1 : n_cutoffs
      T = cluster(Z, 'Cutoff', obj.Cutoff(i), 'Criterion', 'distance', ...
                  'Depth', 2);
      nC = max(T);

      % Find the indices for each cluster.
      C = [];
      for j = 1 : nC
         c = find(T == j);
         C{j} = c';
      end

      centers = zeros(2, nC);
      for j = 1 : nC
         centers(:, j) = mean(obj.XY(C{j}, :), 1);
      end

      % # of clusters
      N(i) = nC;

      % Nm: # of clusters with multiple objects (> 2)
      % Am: area of clusters with multiple objects (> 2)
      % Oa: # of objects per cluster
      % Om: # of objects per cluster with multiple objects (> 2)
      Nm = 0;
      Am = 0;
      for j = 1 : nC
         C_j = C{j};
         nC_j = length(C_j);
         Oa(i) = Oa(i) + nC_j;
         if nC_j > 2
            Nm = Nm + 1;
            Om(i) = Om(i) + nC_j;
            %[~, area] = convhull(obj.XY(C_j, 1), obj.XY(C_j, 2));
            [~, area] = boundary(obj.XY(C_j, 1), obj.XY(C_j, 2), ...
                                 obj.ShrinkFactor);
            Am = Am + area;
         end
      end
      M(i) = Nm;
      A(i) = Am / Nm;
      Oa(i) = Oa(i) / nC;
      Om(i) = Om(i) / Nm;

      % mean nearest neighbor center-to-center distances
      %min_c2c_dists = min(squareform(pdist(centers')) + 1.0e+10 * eye(nC));
      min_c2c_dists = SRcluster.nn_distances(centers');
      if isempty(min_c2c_dists)   % 1 cluster
         D(i) = 0;
      else                        % > 1 cluster
         D(i) = mean(min_c2c_dists);
      end
   end

   n_figs = 0;

   % # of clusters vs. cutoff
   n_figs = n_figs + 1;
   if obj.PlotFigures
      cutoffFigs(n_figs) = figure();
   else
      cutoffFigs(n_figs) = figure('Visible', 'off');
   end
   hold on
   plot(obj.Cutoff, N, 'k-',  'LineWidth', 3);
   plot(obj.Cutoff, M, 'm--', 'LineWidth', 3);
   xlabel('cutoff (nm)');
   ylabel('# of clusters');
   legend('total', '# objs > 2', 'Location', 'NorthEast');
   hold off

   % # of objects per clusters vs. cutoff
   n_figs = n_figs + 1;
   if obj.PlotFigures
      cutoffFigs(n_figs) = figure();
   else
      cutoffFigs(n_figs) = figure('Visible', 'off');
   end
   hold on
   plot(obj.Cutoff, Oa, 'k-',  'LineWidth', 3);
   plot(obj.Cutoff, Om, 'm--', 'LineWidth', 3);
   xlabel('cutoff (nm)');
   ylabel('# of objects per cluster');
   legend('all', '# objs > 2', 'Location', 'SouthEast');
   hold off

   % mean cluster area vs. cutoff
   n_figs = n_figs + 1;
   if obj.PlotFigures
      cutoffFigs(n_figs) = figure();
   else
      cutoffFigs(n_figs) = figure('Visible', 'off');
   end
   hold on
   plot(obj.Cutoff, A, 'k-', 'LineWidth', 3);
   xlabel('cutoff (nm)');
   ylabel('mean cluster area (nm^2)');
   hold off

   % mean nearest neighbor center-to-center distances vs. cutoff
   n_figs = n_figs + 1;
   if obj.PlotFigures
      cutoffFigs(n_figs) = figure();
   else
      cutoffFigs(n_figs) = figure('Visible', 'off');
   end
   hold on
   plot(obj.Cutoff, D, 'k-', 'LineWidth', 3);
   xlabel('cutoff (nm)');
   ylabel('mean nearest neighbor center-to-center distances (nm)');
   hold off

end

% =============================================================================

function E = chooseCutoff(obj)
% Choose a cutoff value from the peak of the curve of the number of clusters
% with multiple objects vs. cutoff distance.  If more than one cutoff value is
% at the peak, choose the largest.
%
% Output:
%    E   largest cutoff at the peak of the curve of the number of clusters with
%        multiple objects vs. cutoff distance

   if size(obj.XY, 1) <= 1
      warning('chooseCutoff: 0 or 1 points provided!');
      return
   end

   Z = linkage(obj.XY, 'single');

   n_cutoffs = length(obj.Cutoff);
   M = zeros(1, n_cutoffs);    % # of clusters with > 2 objects
   for i = 1 : n_cutoffs
      T = cluster(Z, 'Cutoff', obj.Cutoff(i), 'Criterion', 'distance', ...
                  'Depth', 2);
      nC = max(T);

      % Find the indices for each cluster.
      C = [];
      for j = 1 : nC
         c = find(T == j);
         C{j} = c';
      end

      % Nm: # of clusters with multiple objects (> 2)
      Nm = 0;
      for j = 1 : nC
         C_j = C{j};
         nC_j = length(C_j);
         if nC_j > 2
            Nm = Nm + 1;
         end
      end
      M(i) = Nm;
   end

   %[Mmax, k] = max(M);
   Mmax = max(M);
   k = find(M == Mmax);
   E = obj.Cutoff(k);
   E = max(E);
   fprintf('\ncutoff = %.0f nm (%d multiple object clusters)\n', E, Mmax);
   if length(k) > 1
      fprintf('--- chosen from cutoffs =');
      fprintf(' %.0f', obj.Cutoff(k));
      fprintf('\n');
   end

end

% =============================================================================

function SRcollapseFig = plotSRcollapse(obj, xy_region)
% Plot original and collapsed data, using cluster boundary outlines to indicate
% which original observations were collapsed into a single localization.
%
% Input:
%    xy_region       (OPTIONAL) coordinates of the domain vertices [N x 2]
% Output:
%    SRcollapseFig   figure handle

   if obj.PlotFigures
      SRcollapseFig = figure();
   else
      SRcollapseFig = figure('Visible', 'off');
   end
   hold on
   % Fake points used to make the legend.
   plot(-1e+10, -1e+10, 'c.', 'MarkerSize', 10);
   plot(-1e+10, -1e+10, 'b.', 'MarkerSize', 10);
   plot(-1e+10, -1e+10, 'k-');
   % Plot original data.
   plot(obj.XY_orig(:, 1), obj.XY_orig(:, 2), 'c.', 'MarkerSize', 15);
   % Plot processed data.
   plot(obj.XY(:, 1), obj.XY(:, 2), 'b.', 'MarkerSize', 10);
   % Plot convex hull of combined nodes.
   n_eliminated = 0;
   n_collapses = length(obj.Nodes_combined);
   for i = 1 : n_collapses
      i_collapsed = obj.Nodes_combined{i};
      n_collapsed = length(i_collapsed);
      n_eliminated = n_eliminated + n_collapsed - 1;
      if n_collapsed == 2
         plot(obj.XY_orig(i_collapsed, 1), obj.XY_orig(i_collapsed, 2), 'k-');
      else
         %k = convhull(obj.XY_orig(i_collapsed, 1), ...
         %             obj.XY_orig(i_collapsed, 2));
         k = boundary(obj.XY_orig(i_collapsed, 1), ...
                      obj.XY_orig(i_collapsed, 2), obj.ShrinkFactor);
         k = i_collapsed(k);
         plot(obj.XY_orig(k, 1), obj.XY_orig(k, 2), 'k-');
      end
   end

   if exist('xy_region', 'var')
      plot(xy_region(:, 1), xy_region(:, 2), 'k-', 'LineWidth', 3);

      xmin = min([obj.XY_orig(:, 1); xy_region(:, 1)]);
      xmax = max([obj.XY_orig(:, 1); xy_region(:, 1)]);
      ymin = min([obj.XY_orig(:, 2); xy_region(:, 2)]);
      ymax = max([obj.XY_orig(:, 2); xy_region(:, 2)]);
   else
      xmin = min(obj.XY_orig(:, 1));
      xmax = max(obj.XY_orig(:, 1));
      ymin = min(obj.XY_orig(:, 2));
      ymax = max(obj.XY_orig(:, 2));
   end
   dx = xmax - xmin;
   dy = ymax - ymin;
   if isempty(dx) | dx == 0
      dx = 1;
   end
   if isempty(dy) | dy == 0
      dy = 1;
   end
   axis([xmin - 0.05*dx, xmax + 0.05*dx, ymin - 0.05*dy, ymax + 0.05*dy]);
   xlabel('x (nm)');
   ylabel('y (nm)');
   title(sprintf('collapses = %d, eliminated = %d (P > %.3f)', ...
                 n_collapses, n_eliminated, obj.LoS));
   legend('original objects', 'collapsed objects', 'collapse boundaries', ...
          'Location', 'Best');
   hold off

end

% =============================================================================

function SRcollapseFig = plotSRcollapse3(obj, xy_region)
% Plot original and collapsed data, using cluster boundary outlines to indicate
% which original observations were collapsed into a single localization.
%
% Input:
%    xy_region       (OPTIONAL) coordinates of the domain vertices [N x 3]
% Output:
%    SRcollapseFig   figure handle

   if obj.PlotFigures
      SRcollapseFig = figure();
   else
      SRcollapseFig = figure('Visible', 'off');
   end
   hold on
   % Fake points used to make the legend.
   plot3(-1e+10, -1e+10, -1e+10, 'c.', 'MarkerSize', 10);
   plot3(-1e+10, -1e+10, -1e+10, 'b.', 'MarkerSize', 10);
   plot3(-1e+10, -1e+10, -1e+10, 'k-');
   % Plot original data.
   plot3(obj.XY_orig(:, 1), obj.XY_orig(:, 2), obj.XY_orig(:, 3), ...
         'c.', 'MarkerSize', 15);
   % Plot processed data.
   plot3(obj.XY(:, 1), obj.XY(:, 2), obj.XY(:, 3), 'b.', 'MarkerSize', 10);
   % Plot convex hull of combined nodes.
   n_eliminated = 0;
   n_collapses = length(obj.Nodes_combined);
   for i = 1 : n_collapses
      i_collapsed = obj.Nodes_combined{i};
      n_collapsed = length(i_collapsed);
      n_eliminated = n_eliminated + n_collapsed - 1;
      if n_collapsed == 2
         plot3(obj.XY_orig(i_collapsed, 1), obj.XY_orig(i_collapsed, 2), ...
               obj.XY_orig(i_collapsed, 3), 'k-');
      else
         %k = convhull(obj.XY_orig(i_collapsed, 1), ...
         %             obj.XY_orig(i_collapsed, 2));
         k = boundary(obj.XY_orig(i_collapsed, 1), ...
                      obj.XY_orig(i_collapsed, 2), ...
                      obj.XY_orig(i_collapsed, 3), obj.ShrinkFactor);
         k = i_collapsed(k);
         plot3(obj.XY_orig(k, 1), obj.XY_orig(k, 2), obj.XY_orig(k, 3), 'k-');
      end
   end

   if exist('xy_region', 'var')
      plot3(xy_region(:, 1), xy_region(:, 2), xy_region(:, 3), ...
            'k-', 'LineWidth', 3);

      xmin = min([obj.XY_orig(:, 1); xy_region(:, 1)]);
      xmax = max([obj.XY_orig(:, 1); xy_region(:, 1)]);
      ymin = min([obj.XY_orig(:, 2); xy_region(:, 2)]);
      ymax = max([obj.XY_orig(:, 2); xy_region(:, 2)]);
      zmin = min([obj.XY_orig(:, 3); xy_region(:, 3)]);
      zmax = max([obj.XY_orig(:, 3); xy_region(:, 3)]);
   else
      xmin = min(obj.XY_orig(:, 1));
      xmax = max(obj.XY_orig(:, 1));
      ymin = min(obj.XY_orig(:, 2));
      ymax = max(obj.XY_orig(:, 2));
      zmin = min(obj.XY_orig(:, 3));
      zmax = max(obj.XY_orig(:, 3));
   end
   dx = xmax - xmin;
   dy = ymax - ymin;
   dz = zmax - zmin;
   if isempty(dx) | dx == 0
      dx = 1;
   end
   if isempty(dy) | dy == 0
      dy = 1;
   end
   if isempty(dz) | dz == 0
      dz = 1;
   end
   axis([xmin - 0.05*dx, xmax + 0.05*dx, ymin - 0.05*dy, ymax + 0.05*dy, ...
         zmin - 0.05*dz, zmax + 0.05*dz]);
   xlabel('x (nm)');
   ylabel('y (nm)');
   zlabel('z (nm)');
   title(sprintf('collapses = %d, eliminated = %d (P > %.3f)', ...
                 n_collapses, n_eliminated, obj.LoS));
   legend('original objects', 'collapsed objects', 'collapse boundaries', ...
          'Location', 'Best');
   hold off

end

% =============================================================================

function SRclusterFig = plotSRclusters(obj)
% Plot original observations, collapsed localizations, clusters of collapsed
% points that pass the single label test so really should also be collapsed
% (pass 2 of the algorithm), and clusters of collapsed points that are
% considered to be true clusters.
%
% Output:
%    SRclusterFig   figure handle

   if isempty(obj.XY)
      warning('No points to cluster!');
      return
   end

   switch obj.Algorithm
   case 'hierarchalSimple'
      obj.C = SRcluster.hierarchalSimple(obj.XY, obj.E);
      nC = length(obj.C);
   otherwise
      c = Clustering();
      if strcmp(obj.Algorithm, 'Voronoi')
         c.Alpha      = obj.Alpha;
         c.Valgorithm = obj.Valgorithm;
         c.Plotting   = obj.Plotting;
         c.PtIDs      = obj.PtIDs;
      end
      [nC, obj.C, ~, ptsI] = ...
         c.cluster(obj.Algorithm, obj.XY, obj.E, obj.minPts);

      for i = 1 : length(ptsI)
         obj.C{nC + i} = ptsI(i);
      end
      nC = nC + length(ptsI);
   end

   collapsed = false;
   clustered = false;

   % Check each cluster in turn to see if it is likely that it really
   % represents a single label.
   if obj.PlotFigures
      SRclusterFig = figure();
   else
      SRclusterFig = figure('Visible', 'off');
   end
   hold on
   % Fake points used to make the legend.
   plot(-1e+10, -1e+10, 'c.', 'MarkerSize', 10);
   plot(-1e+10, -1e+10, 'b.', 'MarkerSize', 10);
   plot(-1e+10, -1e+10, 'g-', 'MarkerSize', 10);
   plot(-1e+10, -1e+10, 'r-', 'MarkerSize', 10);
   % Plot original data.
   plot(obj.XY_orig(:, 1), obj.XY_orig(:, 2), 'c.', 'MarkerSize', 10);
   % Look for clusters.
   xmin = 1.0e+10;
   ymin = 1.0e+10;
   xmax = -1.0e+10;
   ymax = -1.0e+10;
   for i = 1 : nC
      C_i = obj.C{i};
      xmin = min(xmin, min(obj.XY(C_i, 1)));
      ymin = min(ymin, min(obj.XY(C_i, 2)));
      xmax = max(xmax, max(obj.XY(C_i, 1)));
      ymax = max(ymax, max(obj.XY(C_i, 2)));
      nC_i = length(C_i);
      [Pvalue, XY_wm, sigma_wm] =                                     ...
         SRcluster.singleLabelTest(obj.XY(C_i, :), obj.Sigma(C_i, :), ...
                                   obj.Sigma_Reg);
      if obj.Verbose
         fprintf('[%3d] Cluster of %3d: Pvalue = %f', i, nC_i, Pvalue);
      end
      % A cluster of 1, so already a single label.
      if nC_i == 1
         if obj.Verbose
            fprintf(' *');
         end
      % Collapsed to a single label.
      elseif Pvalue > obj.LoS
         collapsed = true;
         if obj.Verbose
            fprintf(' *');
         end
         % Draw outline.
         if nC_i == 2
            plot(obj.XY(C_i, 1), obj.XY(C_i, 2), 'r-', 'LineWidth', 3);
         elseif nC_i > 2
            %k = convhull(obj.XY(C_i, 1), obj.XY(C_i, 2));
            k = boundary(obj.XY(C_i, 1), obj.XY(C_i, 2), obj.ShrinkFactor);
            plot(obj.XY(C_i(k), 1), obj.XY(C_i(k), 2), 'r-', 'LineWidth', 3);
         end
      % Cluster containing multiple labels.
      else
         clustered = true;
         % Draw outline.
         if nC_i == 2
            plot(obj.XY(C_i, 1), obj.XY(C_i, 2), 'g-', 'LineWidth', 3);
         elseif nC_i > 2
            %k = convhull(obj.XY(C_i, 1), obj.XY(C_i, 2));
            k = boundary(obj.XY(C_i, 1), obj.XY(C_i, 2), obj.ShrinkFactor);
            plot(obj.XY(C_i(k), 1), obj.XY(C_i(k), 2), 'g-', 'LineWidth', 3);
         end
      end
      if obj.Verbose
         fprintf('\n');
      end
   end
   % Plot processed data.
   plot(obj.XY(:, 1), obj.XY(:, 2), 'b.', 'MarkerSize', 10);
   dx = xmax - xmin;
   dy = ymax - ymin;
   % Degenerate cases that show up when nC == 1, for example.
   if dx == 0
      dx = 1;
   end
   if dy == 0
      dy = 1;
   end
   axis([xmin - 0.05*dx, xmax + 0.05*dx, ymin - 0.05*dy, ymax + 0.05*dy]);
   xlabel('x (nm)');
   ylabel('y (nm)');
   title(sprintf('(P > %.3f, cutoff = %.0f)', obj.LoS, obj.E));
   if collapsed
      legend('original objects', 'collapsed objects',          ...
             'multiple label cluster', 'single label cluster', ...
             'Location', 'Best');
   elseif clustered
      legend('original objects', 'collapsed objects',          ...
             'multiple label cluster', 'Location', 'Best');
   else
      legend('original objects', 'collapsed objects', 'Location', 'Best');
   end
   hold off

end

% -----------------------------------------------------------------------------

function SRclusterFig = plotSRclusters3(obj)
% Plot original observations, collapsed localizations, clusters of collapsed
% points that pass the single label test so really should also be collapsed
% (pass 2 of the algorithm), and clusters of collapsed points that are
% considered to be true clusters.
%
% Output:
%    SRclusterFig   figure handle

   if isempty(obj.XY)
      warning('No points to cluster!');
      return
   end

   switch obj.Algorithm
   case 'hierarchalSimple'
      obj.C = SRcluster.hierarchalSimple(obj.XY, obj.E);
      nC = length(obj.C);
   otherwise
      c = Clustering();
      if strcmp(obj.Algorithm, 'Voronoi')
         c.Alpha      = obj.Alpha;
         c.Valgorithm = obj.Valgorithm;
         c.Plotting   = obj.Plotting;
         c.PtIDs      = obj.PtIDs;
      end
      [nC, obj.C, ~, ptsI] = ...
         c.cluster(obj.Algorithm, obj.XY, obj.E, obj.minPts);

      for i = 1 : length(ptsI)
         obj.C{nC + i} = ptsI(i);
      end
      nC = nC + length(ptsI);
   end

   collapsed = false;
   clustered = false;

   % Check each cluster in turn to see if it is likely that it really
   % represents a single label.
   if obj.PlotFigures
      SRclusterFig = figure();
   else
      SRclusterFig = figure('Visible', 'off');
   end
   hold on
   % Fake points used to make the legend.
   plot3(-1e+10, -1e+10, -1e+10, 'c.', 'MarkerSize', 10);
   plot3(-1e+10, -1e+10, -1e+10, 'b.', 'MarkerSize', 10);
   plot3(-1e+10, -1e+10, -1e+10, 'g-', 'MarkerSize', 10);
   plot3(-1e+10, -1e+10, -1e+10, 'r-', 'MarkerSize', 10);
   % Plot original data.
   plot3(obj.XY_orig(:, 1), obj.XY_orig(:, 2), obj.XY_orig(:, 3), ...
         'c.', 'MarkerSize', 10);
   % Look for clusters.
   xmin = 1.0e+10;
   ymin = 1.0e+10;
   zmin = 1.0e+10;
   xmax = -1.0e+10;
   ymax = -1.0e+10;
   zmax = -1.0e+10;
   for i = 1 : nC
      C_i = obj.C{i};
      xmin = min(xmin, min(obj.XY(C_i, 1)));
      ymin = min(ymin, min(obj.XY(C_i, 2)));
      zmin = min(zmin, min(obj.XY(C_i, 3)));
      xmax = max(xmax, max(obj.XY(C_i, 1)));
      ymax = max(ymax, max(obj.XY(C_i, 2)));
      zmax = max(zmax, max(obj.XY(C_i, 3)));
      nC_i = length(C_i);
      [Pvalue, XY_wm, sigma_wm] =                                     ...
         SRcluster.singleLabelTest(obj.XY(C_i, :), obj.Sigma(C_i, :), ...
                                   obj.Sigma_Reg);
      if obj.Verbose
         fprintf('[%3d] Cluster of %3d: Pvalue = %f', i, nC_i, Pvalue);
      end
      % A cluster of 1, so already a single label.
      if nC_i == 1
         if obj.Verbose
            fprintf(' *');
         end
      % Collapsed to a single label.
      elseif Pvalue > obj.LoS
         collapsed = true;
         if obj.Verbose
            fprintf(' *');
         end
         % Draw outline.
         if nC_i == 2
            plot3(obj.XY(C_i, 1), obj.XY(C_i, 2), obj.XY(C_i, 3), ...
                  'r-', 'LineWidth', 3);
         elseif nC_i > 2
            %k = convhull(obj.XY(C_i, 1), obj.XY(C_i, 2));
            k = boundary(obj.XY(C_i, 1), obj.XY(C_i, 2), obj.XY(C_i, 3), ...
                         obj.ShrinkFactor);
            plot3(obj.XY(C_i(k), 1), obj.XY(C_i(k), 2), obj.XY(C_i(k), 3), ...
                  'r-', 'LineWidth', 3);
         end
      % Cluster containing multiple labels.
      else
         clustered = true;
         % Draw outline.
         if nC_i == 2
            plot3(obj.XY(C_i, 1), obj.XY(C_i, 2), obj.XY(C_i, 3), ...
                  'g-', 'LineWidth', 3);
         elseif nC_i > 2
            %k = convhull(obj.XY(C_i, 1), obj.XY(C_i, 2));
            k = boundary(obj.XY(C_i, 1), obj.XY(C_i, 2), obj.XY(C_i, 3), ...
                         obj.ShrinkFactor);
            plot3(obj.XY(C_i(k), 1), obj.XY(C_i(k), 2), obj.XY(C_i(k), 3), ...
                  'g-', 'LineWidth', 3);
         end
      end
      if obj.Verbose
         fprintf('\n');
      end
   end
   % Plot processed data.
   plot3(obj.XY(:, 1), obj.XY(:, 2), obj.XY(:, 3), 'b.', 'MarkerSize', 10);
   dx = xmax - xmin;
   dy = ymax - ymin;
   dz = zmax - zmin;
   axis([xmin - 0.05*dx, xmax + 0.05*dx, ymin - 0.05*dy, ymax + 0.05*dy, ...
         zmin - 0.05*dz, zmax + 0.05*dz]);
   xlabel('x (nm)');
   ylabel('y (nm)');
   zlabel('z (nm)');
   title(sprintf('(P > %.3f, cutoff = %.0f)', obj.LoS, obj.E));
   if collapsed
      legend('original objects', 'collapsed objects',          ...
             'multiple label cluster', 'single label cluster', ...
             'Location', 'Best');
   elseif clustered
      legend('original objects', 'collapsed objects',          ...
             'multiple label cluster', 'Location', 'Best');
   else
      legend('original objects', 'collapsed objects', 'Location', 'Best');
   end
   hold off

end

% =============================================================================

function [results, analysisFigs] = analyzeSRclusters(obj)
% Compute statistics for original observations and collapsed localizations.
% Pass 2 of the algorithm (collapsing any clusters detected after the initial
% collapse that pass the single label test) is performed with the results given
% in the final analysis.
%
% Output:
%   results        structure with the fields below
%   analysisFigs   various handles of figures analyzing the SRcollapse results
%
% numclust_orig(1,2,3)       number of singles, doubles, multiples in original
%                            data (observations)
% singles_per_total_clusters_orig
%                            fraction of singles as clusters
% numobjs_orig(1,2,3)        number of objects in single, double, multiple
%                            clusters
% numobjs_per_multiple_cluster_orig
%                            average number of objects per multiple (>= 3)
%                            cluster
% n_collapses_pass1          number of clusters that were collapsed by the
%                            hierarchal top-down algorithm (pass 1)
% numobjs_collapsed_pass1    number of objects in each cluster before it
%                            was collapsed _including_ singlets
% n_objs_collapsed_pass1     sum(numobjs_collapsed_pass1(1 : n_collapses))
% NOTE: Therefore, the number of objects that were eliminated = 
%          n_objs_collapsed_pass1 - n_collapses_pass1
% n_objs_per_cluster_orig    number of objects per original cluster
% n_objs_per_cluster         number of objects per collapsed cluster
% n_collapses_pass2          number of clusters collapsed in pass 2
% n_objs_collapsed_pass2     total number of objects in the clusters before
%                            they were collapsed
% NOTE: Thus, the number of objects eliminated =
%          n_objs_collapsed_pass2 - n_collapsed_pass2
% numclust(1,2,3)            number of singles, doubles, multiples in the
%                            collapsed data (localizations)
% singles_per_total_clusters fraction of singles as clusters
% numdensity(1,2,3)          numclust(1,2,3) / ROI area
% numobjs(1,2,3)             number of objects in single, double, multiple
%                            clusters
% numobjs_per_multiple_cluster
%                            number of objects per multiple (>= 3) cluster
% cluster_numobjs            number of objects per cluster for all clusters
% cluster_areas              area per cluster for all clusters
% cluster_radii              equivalent radii per cluster for all clusters
% cluster_perimeters         perimeter per cluster for all clusters
% cluster_compactness        4 pi area / perimeter^2 per cluster for all
%                            clusters
% area_per_multiple_cluster  average area per multiple (>= 3) cluster
% radius_per_2orMore_cluster average radius per double and multiple cluster
% numobjs_per_area           number of objects in each multiple cluster
%                            divided by its area
% centers                    centers of each final cluster
% min_c2c_dists              nearest neighbor center-to-center distances
%                            between all cluster centers
% min_c2c_dists_ge3          nearest neighbor center-to-center distances
%                            between centers of clusters of 3 or more objects
% min_c2c_dist(1,2,3,4)      minimum nearest neighbor center-to-center
%                            distances between centers of clusters of
%                            (1) singles, (2) doubles, (3) multiples, (4) all
% singlets                   indices of the original objects forming each
%                            singlet
% multiples                  indices of the original objects forming each
%                            multiple
% singlets_always            fraction of singlets that were always singlets

   if isempty(obj.XY)
      warning('No points to cluster!');
      return
   end

   switch obj.Algorithm
   case 'hierarchalSimple'
      obj.C_orig = SRcluster.hierarchalSimple(obj.XY_orig, obj.E);
      obj.C =      SRcluster.hierarchalSimple(obj.XY,      obj.E);
      nC_orig = length(obj.C_orig);
      nC      = length(obj.C);

      % Find the center of each final cluster.
      centers = zeros(obj.Dim, nC);
      for i = 1 : nC
         centers(:, i) = mean(obj.XY(obj.C{i}, :), 1);
      end

   otherwise
      c = Clustering();
      if strcmp(obj.Algorithm, 'Voronoi')
         c.Alpha      = obj.Alpha;
         c.Valgorithm = obj.Valgorithm;
         c.Plotting   = obj.Plotting;
         c.PtIDs      = obj.PtIDs;
      end
      [nC, obj.C, centers, ptsI]          = ...
         c.cluster(obj.Algorithm, obj.XY,      obj.E, obj.minPts);
      if strcmp(obj.Algorithm, 'DBSCAN_Daszykowski_noE')
         obj.E = c.E;
      end
      [nC_orig, obj.C_orig, ~, ptsI_orig] = ...
         c.cluster(obj.Algorithm, obj.XY_orig, obj.E, obj.minPts);

      for i = 1 : length(ptsI_orig)
         obj.C_orig{nC_orig + i} = ptsI_orig(i);
      end
      nC_orig = nC_orig + length(ptsI_orig);
      for i = 1 : length(ptsI)
         obj.C{nC + i} = ptsI(i);
         centers(:, nC + i) = obj.XY(ptsI(i), :);
      end
      nC = nC + length(ptsI);
   end

   % Compute the number of objects in each cluster.
   nCi_orig = arrayfun(@(x) length(cell2mat(x)), obj.C_orig);
   nCi      = arrayfun(@(x) length(cell2mat(x)), obj.C);

   % (1) single, (2) double, (>2) multiple:
   %
   % Number of clusters containing 1/2/>2 objects.
   numclust_orig = zeros(1, 3);
   numclust      = zeros(1, 3);
   % Total number of objects in clusters containing 1/2/>2 objects.
   numobjs_orig  = zeros(1, 3);
   numobjs       = zeros(1, 3);
   % Number of objects per cluster.
   n_objs = [];
   % Area per cluster.
   areas = [];
   % Equivalent radii per cluster.
   radii = [];
   % Perimeter per cluster.
   perims = [];
   % Compactness per cluster.
   compactness = [];
   % Total area of multiple clusters.
   area = 0;
   % Total equivalent radii of double and multiple clusters.
   radius = 0;
   % Minimum center-to-center distances between objects for clusters containing
   %   1/2/>2/all objects.
   min_c2c_dist  = zeros(1, 4);
   % Number of objects in a cluster per area of the cluster for each multiple
   %    cluster.
   numobjs_per_area = [];

   % --------------------------------------------------------------------------

   % Analyze the original clusters.
   for i = 1 : nC_orig
      nC_i = length(obj.C_orig{i});
      if nC_i == 1
         numclust_orig(1) = numclust_orig(1) + 1;
         numobjs_orig(1)  = numobjs_orig(1) + 1;
      elseif nC_i == 2
         numclust_orig(2) = numclust_orig(2) + 1;
         numobjs_orig(2)  = numobjs_orig(2)  + 2;
      else
         numclust_orig(3) = numclust_orig(3) + 1;
         numobjs_orig(3)  = numobjs_orig(3)  + nC_i;
      end
   end

   % --------------------------------------------------------------------------

   if obj.Printing
      fprintf('\n');
      fprintf('Original Analysis:\n\n');
   end

   results.numclust_orig = numclust_orig;
   if obj.Printing
      fprintf('# of clusters  = %d\n', nC_orig);
      fprintf('# of singles   = %d\n', numclust_orig(1));
      fprintf('# of doubles   = %d\n', numclust_orig(2));
      fprintf('# of multiples = %d\n', numclust_orig(3));
   end

   results.singles_per_total_clusters_orig = numclust_orig(1) / nC_orig;
   if obj.Printing
      fprintf('singles / total clusters = %.3f\n', ...
              results.singles_per_total_clusters_orig);
   end

   results.numobjs_orig = numobjs_orig;
   if numclust_orig(3) ~= 0
      results.numobjs_per_multiple_cluster_orig = ...
         numobjs_orig(3) / numclust_orig(3);
   else
      results.numobjs_per_multiple_cluster_orig = [];
   end
   if obj.Printing
      fprintf('# objects (TOTAL)     = %d\n', sum(numobjs_orig));
      fprintf('# objects (singles)   = %d\n', numobjs_orig(1));
      fprintf('# objects (doubles)   = %d\n', numobjs_orig(2));
      fprintf('# objects (multiples) = %d\n', numobjs_orig(3));
      fprintf('# objects / multiple cluster = %.3f\n', ...
              results.numobjs_per_multiple_cluster_orig);
   end

   % --------------------------------------------------------------------------

   % Analyze the combined nodes (collapse analysis).
   %
   % n_collapses       is the number of clusters that were collapsed by the
   %                   hierarchal top-down algorithm (pass 1)
   % numobjs_collapsed is the number of objects in each cluster before it was
   %                   collapsed _including_ singlets
   % Therefore, the number of objects that were eliminated = 
   %   sum(numobjs_collapsed(1 : n_collapses)) - n_collapses.
   [analysisFigs{1}, n_collapses, numobjs_collapsed] = ...
      obj.collapseHistogram(numobjs_orig);
   results.n_collapses_pass1       = n_collapses;
   results.numobjs_collapsed_pass1 = numobjs_collapsed;
   results.n_objs_collapsed_pass1  = sum(numobjs_collapsed);

   % --------------------------------------------------------------------------

   singlets  = [];
   multiples = [];
   n_singlets  = 0;
   n_multiples = 0;
   % Find the singlets that were originally singlets.
   sing = setdiff(1 : size(obj.XY_orig, 1), [obj.Nodes_combined{:}]);

   n_combined = length(obj.Nodes_combined);
   n_collapsed = 0;
   n_objs_collapsed = 0;
   % Check each final cluster in turn to see if it is likely that it really
   % represents a single label and the analyze the results.
   %
   % n_collapsed      is the number of clusters collapsed in pass 2
   % n_objs_collapsed is the total number of objects in the clusters before
   %                  they were collapsed
   % Thus, the number of objects eliminated = n_objs_collapsed - n_collapsed.
   % singlets         are the indices of the original objects forming each
   %                  singlet
   % multiples        are the indices of the original objects forming each
   %                  multiple
   % singlets_always  fraction of singlets that were always singlets
   for i = 1 : nC
      C_i = obj.C{i};
      nC_i = length(C_i);
      [Pvalue, XY_wm, sigma_wm] =                                     ...
         SRcluster.singleLabelTest(obj.XY(C_i, :), obj.Sigma(C_i, :), ...
                                   obj.Sigma_Reg);
      % A cluster of 1, so already a single label.
      if nC_i == 1
         numclust(1) = numclust(1) + 1;
         numobjs(1)  = numobjs(1) + 1;
         n_singlets = n_singlets + 1;
         if C_i <= n_combined
            % The indices of the original uncombined objects forming a singlet.
            singlets{n_singlets} = obj.Nodes_combined{C_i};
         else
            % This is a bit of a fudge, but the counts are correct.
            singlets{n_singlets} = sing(C_i - n_combined);
         end
      % Collapsed to a single label.
      elseif Pvalue > obj.LoS
         n_collapsed = n_collapsed + 1;
         n_objs_collapsed = n_objs_collapsed + nC_i;
         centers(:, i) = XY_wm';  % update the label's center
         numclust(1) = numclust(1) + 1;
         numobjs(1)  = numobjs(1) + 1;
         n_singlets = n_singlets + 1;
         singlets{n_singlets} = [];
         for j = 1 : nC_i
            if C_i(j) <= n_combined
               singlets{n_singlets} = ...
                  [singlets{n_singlets}, obj.Nodes_combined{C_i(j)}];
            else
               singlets{n_singlets} = ...
                  [singlets{n_singlets}, sing(C_i(j) - n_combined)];
            end
         end
      % Cluster containing multiple labels.
      else
         for j = 1 : nC_i
            n_multiples = n_multiples + 1;
            if C_i(j) <= n_combined
               multiples{n_multiples} = obj.Nodes_combined{C_i(j)};
            else
               multiples{n_multiples} = sing(C_i(j) - n_combined);
            end
         end
         % 2 labels.
         if nC_i == 2
            numclust(2) = numclust(2) + 1;
            numobjs(2)  = numobjs(2)  + 2;
            radius = radius + pdist(obj.XY(C_i, :)) / 2;
         % > 2 labels.
         else
            numclust(3) = numclust(3) + 1;
            numobjs(3)  = numobjs(3)  + nC_i;
            n_objs = [n_objs, nC_i];
            %[~, A] = convhull(obj.XY(C_i, 1), obj.XY(C_i, 2));
            [k, A] = boundary(obj.XY(C_i, 1), obj.XY(C_i, 2), ...
                              obj.ShrinkFactor);
            R = sqrt(A / pi);
            areas = [areas, A];
            radii = [radii, R];
            perim = 0;
            for j = 1 : length(k) - 1
               perim = perim + pdist([obj.XY(C_i(k(j)), :); ...
                                      obj.XY(C_i(k(j + 1)), :)]);
            end
            perims = [perims, perim];
            compactness = [compactness, 4*pi*A / perim^2];
            area = area + A;
            radius = radius + R;
            numobjs_per_area = [numobjs_per_area, nC_i / A];
         end
      end
   end
   results.singlets  = singlets;
   results.multiples = multiples;
   if isempty(singlets)
      results.singlets_always = NaN;
   else
      results.singlets_always = ...
         numel(find(cellfun(@numel, singlets) == 1)) / numel(singlets);
   end

   % Analyze the final singlet and multiple nodes (collapse analysis)
   % like as was done for the pass1 combined nodes.
   [analysisFigs{2}, n_collapses_s, numobjs_collapsed_s] = ...
      obj.collapseHistogramSM(singlets, 'singlet');
   results.n_collapses_s       = n_collapses_s;
   results.numobjs_collapsed_s = numobjs_collapsed_s;
   results.n_objs_collapsed_s  = sum(numobjs_collapsed_s);

   [analysisFigs{3}, n_collapses_m, numobjs_collapsed_m] = ...
      obj.collapseHistogramSM(multiples, 'multiple');
   results.n_collapses_m       = n_collapses_m;
   results.numobjs_collapsed_m = numobjs_collapsed_m;
   results.n_objs_collapsed_m  = sum(numobjs_collapsed_m);

   analysisFigs = analysisFigs(~cellfun('isempty', analysisFigs));

   % --------------------------------------------------------------------------

   if obj.Printing
      fprintf('\n');
      fprintf('Final Analysis:\n\n');
   end

   results.n_objs_per_cluster_orig = nCi_orig;
   results.n_objs_per_cluster      = nCi;

   results.n_collapses_pass2 = n_collapsed;
   results.n_objs_collapsed_pass2 = n_objs_collapsed;
   if obj.Printing
      fprintf('# of clusters collapsed at this level = %d\n', n_collapsed);
      fprintf('# of objects eliminated at this level = %d\n', ...
              n_objs_collapsed - n_collapsed);
   end

   results.numclust = numclust;
   if obj.Printing
      fprintf('# of clusters  = %d\n', nC);
      fprintf('# of singles   = %d\n', numclust(1));
      fprintf('# of doubles   = %d\n', numclust(2));
      fprintf('# of multiples = %d\n', numclust(3));
   end

   results.singles_per_total_clusters = numclust(1) / nC;
   if obj.Printing
      fprintf('singles / total clusters = %.3f\n', ...
              results.singles_per_total_clusters);
   end

   results.numdensity = numclust / obj.A_ROI;
   if obj.Printing
      fprintf('# density (singles)   (#/nm^2) = %g\n', results.numdensity(1));
      fprintf('# density (doubles)   (#/nm^2) = %g\n', results.numdensity(2));
      fprintf('# density (multiples) (#/nm^2) = %g\n', results.numdensity(3));
   end

   results.numobjs = numobjs;
   if numclust(3) ~= 0
      results.numobjs_per_multiple_cluster = numobjs(3) / numclust(3);
   else
      results.numobjs_per_multiple_cluster = [];
   end
   if obj.Printing
      fprintf('# objects (TOTAL)     = %d\n', sum(numobjs));
      fprintf('# objects (singles)   = %d\n', numobjs(1));
      fprintf('# objects (doubles)   = %d\n', numobjs(2));
      fprintf('# objects (multiples) = %d\n', numobjs(3));
      fprintf('# objects / multiple cluster = %.3f\n', ...
              results.numobjs_per_multiple_cluster);
   end

   results.cluster_numobjs = n_objs;
   results.cluster_areas = areas;
   results.cluster_radii = radii;
   results.cluster_perimeters  = perims;
   results.cluster_compactness = compactness;
   if numclust(3) ~= 0
      results.area_per_multiple_cluster = area / numclust(3);
   else
      results.area_per_multiple_cluster = [];
   end
   if numclust(2) + numclust(3) ~= 0
      results.radius_per_2orMore_cluster = ...
         radius / (numclust(2) + numclust(3));
   else
      results.radius_per_2orMore_cluster = [];
   end
   if obj.Printing
      fprintf('mean area / multiple cluster (nm^2) = %f\n', ...
              results.area_per_multiple_cluster);
      fprintf('mean radius / double or multiple cluster (nm) = %f\n', ...
              results.radius_per_2orMore_cluster);
   end

   % List, by cluster, of the number of objects in each multiple cluster
   % divided by its area.
   results.numobjs_per_area = numobjs_per_area;
   if length(numobjs_per_area) > 0
      mean_numobjs_per_area = mean(numobjs_per_area);
   else
      mean_numobjs_per_area = [];
   end
   if obj.Printing
      fprintf('mean(# objects / area per multiple cluster) (nm^-2) = %f\n', ...
              mean_numobjs_per_area);
   end

   results.centers = centers;
   n = size(centers, 2);
   if n > 1
      %min_c2c_dists = min(squareform(pdist(centers')) + 1.0e+10 * eye(nC));
      min_c2c_dists = SRcluster.nn_distances(centers');
      min_c2c_dist(4) = min(min_c2c_dists);
      results.min_c2c_dists = min_c2c_dists;
   else
      results.min_c2c_dists = [];
   end
   c = find(nCi == 1);
   n = length(c);
   if n > 1
      %min_c2c_dists = ...
      %   min(squareform(pdist(centers(:, c)')) + 1.0e+10 * eye(n));
      min_c2c_dists = SRcluster.nn_distances(centers(:, c)');
      min_c2c_dist(1) = min(min_c2c_dists);
   end
   c = find(nCi == 2);
   n = length(c);
   if n > 1
      %min_c2c_dists = ...
      %   min(squareform(pdist(centers(:, c)')) + 1.0e+10 * eye(n));
      min_c2c_dists = SRcluster.nn_distances(centers(:, c)');
      min_c2c_dist(2) = min(min_c2c_dists);
   end
   c = find(nCi >= 3);
   n = length(c);
   if n > 1
      %min_c2c_dists = ...
      %   min(squareform(pdist(centers(:, c)')) + 1.0e+10 * eye(n));
      min_c2c_dists = SRcluster.nn_distances(centers(:, c)');
      min_c2c_dist(3) = min(min_c2c_dists);
      results.min_c2c_dists_ge3 = min_c2c_dists;
   else
      results.min_c2c_dists_ge3 = [];
   end
   results.min_c2c_dist = min_c2c_dist;

   if obj.Printing
      fprintf('minimum c2c distance (TOTAL)     (nm) = %f\n', min_c2c_dist(4));
      fprintf('minimum c2c distance (singles)   (nm) = %f\n', min_c2c_dist(1));
      fprintf('minimum c2c distance (doubles)   (nm) = %f\n', min_c2c_dist(2));
      fprintf('minimum c2c distance (multiples) (nm) = %f\n', min_c2c_dist(3));
   end

end

% =============================================================================

function [collapseHistFig, n_collapses, numobjs_collapsed] = ...
   collapseHistogram(obj, numobjs_orig)
% Produce a histogram of the number of observations collapsed per localization.
%
% Inputs:
%    numobjs_orig   number of objects in each cluster pre-collapse
% Outputs:
%    collapseHistFig   figure handle
%    see below

   % Analyze the combined nodes.
   %
   % n_collapses_1     is the number of clusters that were collapsed by the
   %                   hierarchal top-down algorithm (pass 1) _excluding_
   %                   singlets
   % n_collapses       is the number of clusters that were collapsed by the
   %                   hierarchal top-down algorithm (pass 1) _including_
   %                   singlets
   % numobjs_collapsed is the number of objects in each cluster before it was
   %                   collapsed _including_ singlets
   % Therefore, the number of objects that were eliminated = 
   %   sum(numobjs_collapsed(1 : n_collapses_1)) - n_collapses_1.
   % The number of singlet objects that remained singlets =
   %   sum(numobjs_orig) - (n_eliminated + n_collapses_1)
   % The number of collapsed (>= 2 objects) clusters + number of singlets =
   %   number of original objects - number of objects that were eliminated
   n_collapses_1 = length(obj.Nodes_combined);
   n_collapses = sum(numobjs_orig)                                          ...
      - (sum(arrayfun(@(i) length(obj.Nodes_combined{i}), 1:n_collapses_1)) ...
         - n_collapses_1);
   numobjs_collapsed = ones(1, n_collapses);
   for i = 1 : n_collapses_1
      numobjs_collapsed(i) = length(obj.Nodes_combined{i});
   end

   % --------------------------------------------------------------------------

   if obj.Printing
      fprintf('\n');
      fprintf('Collapse Analysis:\n\n');

      fprintf('# of collapses = %d\n', n_collapses);
      fprintf('# objects eliminated = %d\n', ...
              sum(numobjs_collapsed(1 : n_collapses)) - n_collapses);
   end

   if ~isempty(numobjs_collapsed)
      if obj.Printing
         fprintf('mean # of objects collapsed into 1 = %.3f +- %.3f\n', ...
                 mean(numobjs_collapsed), std(numobjs_collapsed));
      end
      if obj.PlotFigures
         collapseHistFig = figure();
      else
         collapseHistFig = figure('Visible', 'off');
      end
      hold on
      collapsed = sort(numobjs_collapsed);
      hist(collapsed, 1 : max(collapsed));
      xlabel('# of objects collapsed to 1 in a cluster');
      ylabel('# of clusters collapsed');
      title(sprintf('mean # of objects collapsed into 1 = %.3f +- %.3f\n', ...
                    mean(numobjs_collapsed), std(numobjs_collapsed)));
      hold off
   else
      if obj.Printing
         fprintf('mean # of objects collapsed into 1 = 0 +- 0\n');
      end
      collapseHistFig = [];
   end

end

% =============================================================================

function [collapseHistFig, n_collapses, numobjs_collapsed] = ...
   collapseHistogramSM(obj, sm, type)
% Produce a histogram of the number of observations collapsed per localization.
% Here, the histograms are for observations collapsed to either a final singlet
% localization or collapsed to a final cluster of multiple localizations.
%
% Inputs:
%    sm     final singlet or multiplet nodes
%    type   'singlets' or 'multiplets'
% Outputs:
%    collapseHistFig   figure handle
%    see below

   % Analyze the combined nodes.
   %
   % n_collapses       is the number of clusters that were collapsed by the
   %                   hierarchal top-down algorithm (pass 1 and 2)
   % numobjs_collapsed is the number of objects in each cluster before it was
   %                   collapsed (pass 1 and 2)
   % Therefore, the number of objects that were eliminated = 
   %   sum(numobjs_collapsed(1 : n_collapses)) - n_collapses.
   n_collapses = length(sm);
   numobjs_collapsed = ones(1, n_collapses);
   for i = 1 : n_collapses
      numobjs_collapsed(i) = length(sm{i});
   end

   % --------------------------------------------------------------------------

   if obj.Printing
      fprintf('\n');
      fprintf('%s Collapse Analysis:\n\n', type);

      fprintf('# of collapses = %d\n', n_collapses);
      fprintf('# objects eliminated = %d\n', ...
              sum(numobjs_collapsed(1 : n_collapses)) - n_collapses);
   end

   if ~isempty(numobjs_collapsed)
      if obj.Printing
         fprintf('mean # of objects collapsed into 1 %s = %.3f +- %.3f\n', ...
                 type, mean(numobjs_collapsed), std(numobjs_collapsed));
      end
      if obj.PlotFigures
         collapseHistFig = figure();
      else
         collapseHistFig = figure('Visible', 'off');
      end
      hold on
      collapsed = sort(numobjs_collapsed);
      hist(collapsed, 1 : max(collapsed));
      xlabel('# of objects collapsed to 1 in a cluster');
      ylabel('# of clusters collapsed');
      title(                                                               ...
         sprintf('mean # of objects collapsed into 1 %s = %.3f +- %.3f\n', ...
                 type, mean(numobjs_collapsed), std(numobjs_collapsed)));
      hold off
   else
      if obj.Printing
         fprintf('mean # of objects collapsed into 1 %s = 0 +- 0\n', type);
      end
      collapseHistFig = [];
   end

end

% =============================================================================

function [results, SRclusterFig] = ...
   analyzeSRclusters_simulated(obj, center, xmax, ymax, xy_region)
% Plot observations, collapsed localizations and true emitters (true
% localizatons) for simulated 2D data.
%
% Inputs:
%    center         original cluster centers
%    x_max          maximum x extent of the simulation ROI
%    y_max          maximum y extent of the simulation ROI
%    xy_region      (OPTIONAL) rectangular region defining the maximum extent
% Outputs:
%    results
%       nearest_true_dists   distance of nearest true emitter to each collapsed
%                            localization
%    SRclusterFig   figure handle

   if isempty(obj.XY)
      warning('No points to cluster!');
      return
   end

   switch obj.Algorithm
   case 'hierarchalSimple'
      obj.C = SRcluster.hierarchalSimple(obj.XY, obj.E);
      nC = length(obj.C);
   otherwise
      c = Clustering();
      if strcmp(obj.Algorithm, 'Voronoi')
         c.Alpha      = obj.Alpha;
         c.Valgorithm = obj.Valgorithm;
         c.Plotting   = obj.Plotting;
         c.PtIDs      = obj.PtIDs;
      end
      [nC, obj.C, ~, ptsI] = ...
         c.cluster(obj.Algorithm, obj.XY, obj.E, obj.minPts);

      for i = 1 : length(ptsI)
         obj.C{nC + i} = ptsI(i);
      end
      nC = nC + length(ptsI);
   end

   if obj.PlotFigures
      SRclusterFig = figure();
   else
      SRclusterFig = figure('Visible', 'off');
   end
   hold on
   if obj.ProduceLegend
      % Fake points used to make the legend.
      plot(-1e+10, -1e+10, 'c.', 'MarkerSize', 10);
      plot(-1e+10, -1e+10, 'b.', 'MarkerSize', 10);
      plot(-1e+10, -1e+10, 'r.', 'MarkerSize', 10);
   end
   % Plot original data.
   plot(obj.XY_orig(:, 1), obj.XY_orig(:, 2), 'c.', 'MarkerSize', 10);
   % Plot processed data.
   plot(obj.XY(:, 1), obj.XY(:, 2), 'b.', 'MarkerSize', 10);
   % Plot original centers.
   plot(center(:, 1), center(:, 2), 'ro', 'MarkerSize', 10);
   if exist('xy_region', 'var')
      plot(xy_region(:, 1), xy_region(:, 2), 'k-', 'LineWidth', 3);

      axis([min([0; xy_region(:, 1)]), max([xmax; xy_region(:, 1)]), ...
            min([0; xy_region(:, 2)]), max([ymax; xy_region(:, 2)])]);
   else
      axis([0, xmax, 0, ymax]);
   end
   xlabel('x (nm)');
   ylabel('y (nm)');
   title(sprintf(''));
   if obj.ProduceLegend
      legend('observations', 'collapsed localizations', 'true emitters', ...
             'Location', 'Best');
   end
   hold off

   % Find the nearest true emitter to each collapsed localization.
   nearest_true_dists = [];
   for i = 1 : size(obj.XY, 1)
      d2 = (obj.XY(i, 1) - center(:, 1)).^2 + (obj.XY(i, 2) - center(:, 2)).^2;
      nearest_true_dists = [nearest_true_dists, sqrt(min(d2))];
   end
   results.nearest_true_dists = nearest_true_dists;

end

% =============================================================================

function [results, SRclusterFig] = ...
   analyzeSRclusters_simulated3(obj, center, xmax, ymax, zmax, xy_region)
% Plot observations, collapsed localizations and true emitters (true
% localizatons) for simulated 3D data.
%
% Inputs:
%    center         original cluster centers
%    x_max          maximum x extent of the simulation ROI
%    y_max          maximum y extent of the simulation ROI
%    z_max          maximum z extent of the simulation ROI
%    xy_region      (OPTIONAL) rectangular region defining the maximum extent
% Outputs:
%    results
%       nearest_true_dists   distance of nearest true emitter to each collapsed
%                            localization
%    SRclusterFig   figure handle

   if isempty(obj.XY)
      warning('No points to cluster!');
      return
   end

   switch obj.Algorithm
   case 'hierarchalSimple'
      obj.C = SRcluster.hierarchalSimple(obj.XY, obj.E);
      nC = length(obj.C);
   otherwise
      c = Clustering();
      if strcmp(obj.Algorithm, 'Voronoi')
         c.Alpha      = obj.Alpha;
         c.Valgorithm = obj.Valgorithm;
         c.Plotting   = obj.Plotting;
         c.PtIDs      = obj.PtIDs;
      end
      [nC, obj.C, ~, ptsI] = ...
         c.cluster(obj.Algorithm, obj.XY, obj.E, obj.minPts);

      for i = 1 : length(ptsI)
         obj.C{nC + i} = ptsI(i);
      end
      nC = nC + length(ptsI);
   end

   if obj.PlotFigures
      SRclusterFig = figure();
   else
      SRclusterFig = figure('Visible', 'off');
   end
   hold on
   if obj.ProduceLegend
      % Fake points used to make the legend.
      plot3(-1e+10, -1e+10, -1e+10, 'c.', 'MarkerSize', 10);
      plot3(-1e+10, -1e+10, -1e+10, 'b.', 'MarkerSize', 10);
      plot3(-1e+10, -1e+10, -1e+10, 'r.', 'MarkerSize', 10);
   end
   % Plot original data.
   plot3(obj.XY_orig(:, 1), obj.XY_orig(:, 2), obj.XY_orig(:, 3), ...
         'c.', 'MarkerSize', 10);
   % Plot processed data.
   plot3(obj.XY(:, 1), obj.XY(:, 2), obj.XY(:, 3), 'b.', 'MarkerSize', 10);
   % Plot original centers.
   plot3(center(:, 1), center(:, 2), center(:, 3), 'ro', 'MarkerSize', 10);
   if exist('xy_region', 'var')
      plot3(xy_region(:, 1), xy_region(:, 2), xy_region(:, 3), ...
            'k-', 'LineWidth', 3);

      axis([min([0; xy_region(:, 1)]), max([xmax; xy_region(:, 1)]), ...
            min([0; xy_region(:, 2)]), max([ymax; xy_region(:, 2)]), ...
            min([0; xy_region(:, 3)]), max([ymax; xy_region(:, 3)])]);
   else
      axis([0, xmax, 0, ymax, 0, zmax]);
   end
   xlabel('x (nm)');
   ylabel('y (nm)');
   ylabel('z (nm)');
   title(sprintf(''));
   if obj.ProduceLegend
      legend('observations', 'collapsed localizations', 'true emitters', ...
             'Location', 'Best');
   end
   hold off

   % Find the nearest true emitter to each collapsed localization.
   nearest_true_dists = [];
   for i = 1 : size(obj.XY, 1)
      d2 = (obj.XY(i, 1) - center(:, 1)).^2 + (obj.XY(i, 2) - center(:, 2)).^2;
      nearest_true_dists = [nearest_true_dists, sqrt(min(d2))];
   end
   results.nearest_true_dists = nearest_true_dists;

end

% =============================================================================

function Reject = reject(obj, SRtest)
% Show statistics of rejected fits for an SRtest object.
%
% Input:
%   SRtest   SRtest object
% Output:
%   Reject   structure containing various statistics (see below)

   Reject = [];

   SRT = struct(SRtest);
   Params     = SRT.Params;
   Results    = SRT.Results;
   Statistics = SRT.Statistics;

   n_files         = length(SRT.FileList);
   n_frames        = max(SRT.Params(:, 5)) + 1;
   AcceptedFits    = Statistics.AcceptedFits;
   PostConnectFits = Statistics.PostConnectFits;
   TotalFits       = Statistics.TotalFits;

   fprintf('TotalFits       = %7d   (%.3f)\n', TotalFits, 1);
   fprintf('AcceptedFits    = %7d   (%.3f)\n', AcceptedFits, ...
           AcceptedFits / TotalFits);
   if SRT.FrameConnect == 1
      fprintf('PostConnectFits = %7d   (%.3f)\n', PostConnectFits, ...
              PostConnectFits / TotalFits);
   end
   fprintf('\n');
   fprintf('MeanPhotons = %.3f\n', Statistics.MeanPhotons);
   fprintf('MeanBG      = %.3f\n', Statistics.MeanBg);
   fprintf('\n');

   X2_CDF = inline('gammainc(x/2,k/2)', 'k', 'x');
   k = SRT.FitBoxSize^2 - 4;
   X2 = -2*SRT.LL;
   pvalue = 1 - X2_CDF(k, X2);

   rejected = zeros(TotalFits, 4);
   rejected(:, 1) = SRT.MinPValue > pvalue;
   rejected(:, 2) = SRT.Params(:, 3)  <= SRT.MinPhotons;
   xMaxCRLBSigma = SRT.CRLB_STD(:, 1) >= SRT.MaxCRLBSigma;
   yMaxCRLBSigma = SRT.CRLB_STD(:, 2) >= SRT.MaxCRLBSigma;
   rejected(:, 3) = xMaxCRLBSigma | yMaxCRLBSigma;
   rejected(:, 4) = SRT.Params(:, 4)  >= SRT.MaxBg;

   total_rejects = sum(rejected, 2);
   %fprintf('Accepted = %d\n', sum(total_rejects == 0));
   fprintf('Rejected:\n');
   s = sum(rejected(:, 1));
   fprintf('   PValue       = %7d   (%.3f)       [>  %.3f]\n', ...
           s, s / TotalFits, SRT.MinPValue);
   s = sum(rejected(:, 2));
   fprintf('   MinPhotons   = %7d   (%.3f)       [<= %d]\n', ...
           s, s / TotalFits, SRT.MinPhotons);
   s = sum(rejected(:, 3));
   fprintf('   MaxCRLBSigma = %7d   (%.3f)       [>= %.3f]\n', ...
           s, s / TotalFits, SRT.MaxCRLBSigma);
   s = sum(xMaxCRLBSigma);
   fprintf('      xMaxCRLBSigma = %7d   (%.3f)\n', s, s / TotalFits);
   s = sum(yMaxCRLBSigma);
   fprintf('      yMaxCRLBSigma = %7d   (%.3f)\n', s, s / TotalFits);
   s = sum(rejected(:, 4));
   fprintf('   MaxBG        = %7d   (%.3f)       [>= %d]\n', ...
           s, s / TotalFits, SRT.MaxBg);
   s = sum(total_rejects > 1);
   fprintf('   --------------------------------\n');
   fprintf('   MULTIPLE     = %7d   (%.3f)\n', s, s / TotalFits);
   s = sum(total_rejects == 2);
   if s > 0
      fprintf('      2 criteria    = %7d   (%.3f)\n', s, s / TotalFits);
   end
   s = sum(total_rejects == 3);
   if s > 0
      fprintf('      3 criteria    = %7d   (%.3f)\n', s, s / TotalFits);
   end
   s = sum(total_rejects == 4);
   if s > 0
      fprintf('      4 criteria    = %7d   (%.3f)\n', s, s / TotalFits);
   end
   s = sum(total_rejects == 5);
   if s > 0
      fprintf('      5 criteria    = %7d   (%.3f)\n', s, s / TotalFits);
   end
   fprintf('   --------------------------------\n');
   s = AcceptedFits - PostConnectFits;
   fprintf('   Connect      = %7d   (%.3f)\n', s, s / TotalFits);

   if obj.Histograms
      figure();
      hold on
      hist(pvalue, 100);
      title(sprintf('Accept PValue <= MinPValue (%.3f)', SRT.MinPValue));
      ax = axis;   ymin = ax(3);   ymax = ax(4);
      plot([SRT.MinPValue, SRT.MinPValue], [ymin, ymax], 'r--');
      hold off

      figure();
      hold on
      hist(SRT.Params(:, 3), 100);
      title(sprintf('Accept Photons > MinPhotons (%d)', SRT.MinPhotons));
      ax = axis;   ymin = ax(3);   ymax = ax(4);
      plot([SRT.MinPhotons, SRT.MinPhotons], [ymin, ymax], 'r--');
      hold off

      figure();
      hold on
      hist(max(SRT.CRLB_STD(:, 1), SRT.CRLB_STD(:, 2)), 100);
      title(sprintf('Accept CRLBSigma < MaxCRLBSigma (%.3f)', ...
                    SRT.MaxCRLBSigma));
      ax = axis;   ymin = ax(3);   ymax = ax(4);
      plot([SRT.MaxCRLBSigma, SRT.MaxCRLBSigma], [ymin, ymax], 'r--');
      hold off

      figure();
      hold on
      hist(SRT.Params(:, 4), 100);
      title(sprintf('Accept Bg <= MaxBg (%d)', SRT.MaxBg));
      ax = axis;   ymin = ax(3);   ymax = ax(4);
      plot([SRT.MaxBg, SRT.MaxBg], [ymin, ymax], 'r--');
      hold off
   end

   SRcluster.plot_All_Localizations(pvalue, SRT);
%  SRcluster.plot_All_Localizations(pvalue, SRT, 2, 2);

   if obj.DetailedStatistics
      % Accepted localizations.
      mask = SRT.MinPValue <= pvalue;
      mask = mask & (SRT.Params(:, 3)   > SRT.MinPhotons);
      mask = mask & (SRT.CRLB_STD(:, 1) < SRT.MaxCRLBSigma);
      mask = mask & (SRT.CRLB_STD(:, 2) < SRT.MaxCRLBSigma);
      mask = mask & (SRT.Params(:, 4)   < SRT.MaxBg);

      %x = zeros(n_files, n_frames);
      %y = zeros(n_files, n_frames);
      n_files_seq  = 1 : n_files;
      n_frames_seq = 1 : n_frames;
      %for i = n_files_seq
      %   x(i, :) = n_frames_seq;
      %end
      %for j = n_frames_seq
      %   y(:, j) = n_files_seq;
      %end
      Reject.Statistics = zeros(n_files, n_frames, 4);
      for i = n_files_seq
         all_file_mask = Reject.all_FileNum == i;
         file_mask     = Reject.FileNum == i;
         for j = 0 : n_frames - 1
            all_frame_mask = all_file_mask & Reject.all_FrameNum == j;
            frame_mask     = file_mask     & Reject.FrameNum == j;
            sum_all      = sum(all_frame_mask);
            sum_rejected = sum(frame_mask);
            Reject.Statistics(i, j + 1, :) = [i, j, sum_all, sum_rejected];
         end
      end
      %figure();
      %plot3(x, y, Reject.Statistics(:, :, 3) - Reject.Statistics(:, :, 4), ...
      %      'ko');

      AbsoluteFrame = max(Params(:, 5)).*(Params(:, 6) - 1) + Params(:, 5);
      % Rejected localizations.
      RejectMask = ~mask;

      Reject.all_X        = SRT.Params(:, 1);
      Reject.all_Y        = SRT.Params(:, 2);
      Reject.all_FrameNum = SRT.Params(:, 5);
      Reject.all_FileNum  = SRT.Params(:, 6);

      Reject.Photons       = Params(RejectMask, 3);
      Reject.Pvalue        = pvalue(RejectMask);
      Reject.BG            = Params(RejectMask, 4);
      Reject.FrameNum      = Params(RejectMask, 5);
      Reject.FileNum       = Params(RejectMask, 6);
      Reject.AbsoluteFrame = AbsoluteFrame(RejectMask);
   end

end

% -----------------------------------------------------------------------------

function [XY_combined, sigma_combined, nodes_combined] = ...
   hierarchalSingleLabelP(obj, XY, sigma, Sigma_Reg, LoS)
% Collect Pvalues from the singleLabelTest vs. linkage distances
%
% Combine multiple clustered points into single labels when appropriate via a
% top-down descent through a hierarchal dendrogram relationship between points.
% n is the original number of points and m is the dimension.
% n' is the final number of points after combinations have occurred.
%
% Inputs:
%    XY          n x m matrix of coordinates
%    sigma       n x m matrix of position uncertainties (1 standard deviation)
%    Sigma_Reg   1 x m array of registration error (1 standard deviation)
%    LoS         level of significance
%
% Outputs:
%    XY_combined      n' x m final coordinate matrix
%    sigma_combined   n' x m final position uncertainty matrix
%    nodes_combined   cell array of indices of combined points per cluster

   XY_combined    = XY;
   sigma_combined = sigma;
   nodes_combined = {};

   n_pts   = size(XY, 1);
   n_nodes = n_pts - 1;   % does not include leaf nodes

   if n_pts <= 1
      return
   end

   Z = linkage(XY, 'single');
   %figure(); dendrogram(Z);

   % Find all the leaf nodes contained by each composite node by parsing the
   % Z matrix (see MATLAB linkage documentation).  Note that composite nodes
   % are indexed as node # - n_pts.  E.g., if n_pts = 10, then Z(2, :) is
   % composite node 2 and overall node 12.  Z(2, 1:2) are the nodes (in
   % overall node numbering) that are joined by node 12 which, for example,
   % might be 7 (a leaf node since it is <= 10) and 11 (a composite node whose
   % components are given in Z(1, 1:2)).
   cn = cell(1, n_nodes);
   ld = zeros(n_nodes, 2);
   for i = 1 : n_nodes
      l1 = Z(i, 1);   % child node 1
      l2 = Z(i, 2);   % child node 2
      if l1 <= n_pts
         leaf_nodes = l1;
      else
         leaf_nodes = cn{l1 - n_pts};
      end
      if l2 <= n_pts
         leaf_nodes = [leaf_nodes, l2];
      else
         leaf_nodes = [leaf_nodes, cn{l2 - n_pts}];
      end
      cn{i} = sort(leaf_nodes);
      ld(i, 1) = Z(i, 3);   % linkage distance
      ld(i, 2) = numel(cn{i});
   end

   PD = [];
   k = 0;
   % Start with the top-level node and descend down through the tree (the tree
   % is taken to have its root at the top and its leaves at the bottom).
   for i = n_nodes : -1 : 1
      LN = cn{i};
%     if ~isempty(LN)
      if ~isempty(LN) & ld(i, 1) < obj.MaxLD & ld(i, 2) < obj.MaxLN
         [Pvalue, XY_wm, sigma_wm] = ...
            SRcluster.singleLabelTest(XY(LN, :), sigma(LN, :), Sigma_Reg);
         PD = [PD; Pvalue, ld(i, :)];
         % Delete leaf nodes that are part of a composite single label.
         % Retain the first leaf node to hold the new information.
         if Pvalue > LoS
            cn{i} = LN(1);
            LN_deleted = LN(2 : end);
            % To be deleted at the end, but in order to keep the numbering
            % unchanged within the loop, set to a fantasy value for now.
            XY_combined(LN_deleted, :)    = NaN;
            sigma_combined(LN_deleted, :) = NaN;
            % Replace XY and sigma for this first leaf node with the weighted
            % means of the collapsed cluster.
            XY_combined(LN(1), :)    = XY_wm;
            sigma_combined(LN(1), :) = sigma_wm;
            k = k + 1;
            nodes_combined{k} = LN;
            for j = 1 : i - 1
               % cn{j} = setdiff(cn{j}, LN);
               % --->
               % cn{j} = cn{j}(~ismember(cn{j}, LN));
               % Optimized version of the line above: 
               cn{j} = cn{j}(~SRcluster.my_ismemberBuiltinTypes(cn{j}, LN));
            end
         end
      end
   end

   XY_combined(isnan(XY_combined(:, 1)), :) = [];
   sigma_combined(isnan(sigma_combined(:, 1)), :) = [];

   obj.PD = PD;

   x_max = ceil(max(PD(find(PD(:, 1) > LoS), 2)));

   if obj.PlotFigures
      obj.PDfig(1) = figure();
   else
      obj.PDfig(1) = figure('Visible', 'off');
   end
   hold on
   plot(PD(:, 2), PD(:, 1), 'k.', 'MarkerSize', 10);
   plot([0, max(PD(:, 2))], [LoS, LoS], 'r--');
   title(sprintf('Collapsing Pvalue > LoS = %.3f', LoS));
   xlabel('linkage distance');
   ylabel('Pvalue');
   axis([0, x_max, 0, 1]);
   hold off

   x_max = ceil(max(PD(find(PD(:, 1) > LoS), 3)));

   if obj.PlotFigures
      obj.PDfig(2) = figure();
   else
      obj.PDfig(2) = figure('Visible', 'off');
   end
   hold on
   plot(PD(:, 3), PD(:, 1), 'k.', 'MarkerSize', 10);
   plot([0, max(PD(:, 2))], [LoS, LoS], 'r--');
   title(sprintf('Collapsing Pvalue > LoS = %.3f', LoS));
   xlabel('# leaf nodes');
   ylabel('Pvalue');
   axis([0, x_max, 0, 1]);
   hold off

   k = find(PD(:, 1) > LoS);
   j = setdiff(1 : size(PD, 1), k);

   if obj.PlotFigures
      obj.PDfig(3) = figure();
   else
      obj.PDfig(3) = figure('Visible', 'off');
   end
   hold on
   cm = jet;
   n = size(cm, 1);
   colormap(cm);
   plot(PD(j, 3), PD(j, 2), 'k+');
   for i = 1 : numel(k)
      plot(PD(k(i), 3), PD(k(i), 2), '.', 'MarkerSize', 10, ...
           'MarkerEdgeColor', cm(round((n - 1) * PD(k(i), 1) + 1), :));
   end
   colorbar;
   title(sprintf('Collapsing Pvalue > LoS = %.3f', LoS));
   xlabel('# leaf nodes');
   ylabel('linkage distance');
   hold off

end

% =============================================================================
end % methods

methods(Static)
% =============================================================================

function min_d = nn_distances(xy)
% Minimum nearest neighbor distances from each point in xy to the other points
%
% Inputs:
%    xy      n x 2 (or n x 3) set of coordinates.
% Outputs:
%    min_d   minimum nearest neighbor distance for each point in xy

   %x_r = xy(:, 1);
   %y_r = xy(:, 2);
   %min_d = min(squareform(pdist([x_r, y_r])) + 1.0e+10 * eye(numel(x_r)));
   % Below is much faster, less memory intensive and works in 3D as well as 2D!
   [~, min_d] = knnsearch(xy, xy, 'K', 2);
   min_d = min_d(:, 2)';

end

% -----------------------------------------------------------------------------

function [XY_combined, sigma_combined, nodes_combined] = ...
   hierarchalSingleLabel(XY, sigma, Sigma_Reg, LoS)
% Combine multiple clustered points into single labels when appropriate via a
% top-down descent through a hierarchal dendrogram relationship between points.
% n is the original number of points and m is the dimension.
% n' is the final number of points after combinations have occurred.
%
% Inputs:
%    XY          n x m matrix of coordinates
%    sigma       n x m matrix of position uncertainties (1 standard deviation)
%    Sigma_Reg   1 x m array of registration error (1 standard deviation)
%    LoS         level of significance
%
% Outputs:
%    XY_combined      n' x m final coordinate matrix
%    sigma_combined   n' x m final position uncertainty matrix
%    nodes_combined   cell array of indices of combined points per cluster

   XY_combined    = XY;
   sigma_combined = sigma;
   nodes_combined = {};

   n_pts   = size(XY, 1);
   n_nodes = n_pts - 1;   % does not include leaf nodes

   if n_pts <= 1
      return
   end

   Z = linkage(XY, 'single');
   %figure(); dendrogram(Z);

   % Find all the leaf nodes contained by each composite node by parsing the
   % Z matrix (see MATLAB linkage documentation).  Note that composite nodes
   % are indexed as node # - n_pts.  E.g., if n_pts = 10, then Z(2, :) is
   % composite node 2 and overall node 12.  Z(2, 1:2) are the nodes (in
   % overall node numbering) that are joined by node 12 which, for example,
   % might be 7 (a leaf node since it is <= 10) and 11 (a composite node whose
   % components are given in Z(1, 1:2)).
   cn = cell(1, n_nodes);
   for i = 1 : n_nodes
      l1 = Z(i, 1);   % child node 1
      l2 = Z(i, 2);   % child node 2
      if l1 <= n_pts
         leaf_nodes = l1;
      else
         leaf_nodes = cn{l1 - n_pts};
      end
      if l2 <= n_pts
         leaf_nodes = [leaf_nodes, l2];
      else
         leaf_nodes = [leaf_nodes, cn{l2 - n_pts}];
      end
      cn{i} = sort(leaf_nodes);
   end

   k = 0;
   % Start with the top-level node and descend down through the tree (the tree
   % is taken to have its root at the top and its leaves at the bottom).
   for i = n_nodes : -1 : 1
      LN = cn{i};
      if ~isempty(LN)
         [Pvalue, XY_wm, sigma_wm] = ...
            SRcluster.singleLabelTest(XY(LN, :), sigma(LN, :), Sigma_Reg);
         % Delete leaf nodes that are part of a composite single label.
         % Retain the first leaf node to hold the new information.
         if Pvalue > LoS
            cn{i} = LN(1);
            LN_deleted = LN(2 : end);
            % To be deleted at the end, but in order to keep the numbering
            % unchanged within the loop, set to a fantasy value for now.
            XY_combined(LN_deleted, :)    = NaN;
            sigma_combined(LN_deleted, :) = NaN;
            % Replace XY and sigma for this first leaf node with the weighted
            % means of the collapsed cluster.
            XY_combined(LN(1), :)    = XY_wm;
            sigma_combined(LN(1), :) = sigma_wm;
            k = k + 1;
            nodes_combined{k} = LN;
            for j = 1 : i - 1
               % cn{j} = setdiff(cn{j}, LN);
               % --->
               % cn{j} = cn{j}(~ismember(cn{j}, LN));
               % Optimized version of the line above: 
               cn{j} = cn{j}(~SRcluster.my_ismemberBuiltinTypes(cn{j}, LN));
            end
         end
      end
   end

   XY_combined(isnan(XY_combined(:, 1)), :) = [];
   sigma_combined(isnan(sigma_combined(:, 1)), :) = [];

end

% -----------------------------------------------------------------------------

function [lia] = my_ismemberBuiltinTypes(a,b)
% Extracted and simplified from MATLAB's ismember.m for the above usage.
% General handling.
% Use FIND method for very small sizes of the input vector to avoid SORT.
% Handle empty arrays and scalars.  
numelA = numel(a);
numelB = numel(b);
if numelA == 0 || numelB <= 1
    if numelA > 0 && numelB == 1
        lia = (a == b);
    else
        lia = false(size(a));
    end
    return
end

scalarcut = 5;
if numelA <= scalarcut
    lia = false(size(a));
        for i=1:numelA
            lia(i) = any(a(i)==b(:));   % ANY returns logical.
        end
else
    % Use method which sorts list, then performs binary search.
    % Convert to full to work in C helper.
    
%       % Find out whether list is presorted before sort
%       % If the list is short enough, SORT will be faster than ISSORTED
%       % If the list is longer, ISSORTED can potentially save time
%       checksortcut = 1000;
%       if numelB > checksortcut
%           sortedlist = issorted(b(:));
%       else
%           sortedlist = 0;
%       end
%       if ~sortedlist
%           b = sort(b(:));
%       end
    
    % Use builtin helper function ISMEMBERHELPER:
    % [LIA,LOCB] = ISMEMBERHELPER(A,B) Returns logical array LIA indicating
    % which elements of A occur in B and a double array LOCB with the
    % locations of the elements of A occuring in B. If multiple instances
    % occur, the first occurence is returned. B must be already sorted.
    
            lia = builtin('_ismemberhelper',a,b);
end
end

% -----------------------------------------------------------------------------

function [X_Combined, S_Combined, nodes_Combined] = ...
   keith_HC(X_All, Sigma_All, Sigma_Reg, LoS)
% Keith Lidke's version of H-SET.

   X_Combined = X_All;
   S_Combined = Sigma_All;
   nodes_combined = {};
   NObs = size(X_All,1);
   k = 0;
   for cc =2:NObs  %Need to treat '1' case seperately (not yet done)
       T  = clusterdata(X_Combined,cc);
       
       for tt =1:max(T)
           ID = (T==tt);
           X = X_Combined(ID,:);
           S = S_Combined(ID,:);
           [Pvalue, X_Point, Sigma_Point] = ...
              SRcluster.singleLabelTest(X,S,Sigma_Reg);
           
           if Pvalue > LoS % Combine all points in cluster
               F = find(ID); %Find first point with this cluster ID
               ID(F(1))=0;        %set to zero so we don't erase in next step
               X_Combined(ID,:)=[];  %Remove others from array
               S_Combined(ID,:)=[];
               X_Combined(F(1),:) = X_Point;
               S_Combined(F(1),:) = Sigma_Point;
               k = k + 1;
               nodes_Combined{k} = F';
           end
           T  = clusterdata(X_Combined,cc);
       end
       
       if cc>size(X_Combined,1); break;end
   end

end

% =============================================================================

function C = hierarchalSimple(XY, E)
% Form clusters such that any point in a cluster is within E of some other
% point in the same cluster.
%
% Inputs:
%    XY       n x m matrix of coordinates
%    E        cutoff distance
%
% Outputs:
%    C        cell array of indices (wrt XY) per cluster

   n = size(XY, 1);
   if n <= 1
      % 1 point so 1 cluster!
      if n == 1
         C{1} = 1;
      % No points so no clusters!
      else
         C{1} = [];
      end
      return
   end

   Z = linkage(XY, 'single');
   %figure(); dendrogram(Z);
   T = cluster(Z, 'Cutoff', E, 'Criterion', 'distance', 'Depth', 2);
   nC = max(T);

   C = [];
   for i = 1 : nC
      c = find(T == i);
      C{i} = c';
   end

end

% =============================================================================

function [Pvalue, X_Point, Sigma_Point] = singleLabelTest(X, Sigma, Sigma_Reg)
%singleLabelTest Tests if a cluster of points came from point source 
%   This function calculates a p-value for a cluster of points. 
%   The meaning of of the p-value is the probabilty that a more extreme set
%   of N points came from a point source, where N is the number of observed
%   points. 
%
% INPUTS
%   X,Y  NxM Vectors of positions. N is number of particles, M is dimension
%   X_Sigma, Y_Sigma    NxM Vectors of position uncertainty (1 STD)
%   Sigma_Reg           1xM array of registration error (1 sigma)
% OUTPUTS
%   Pvalue:             Probability of more exterme cluster
%   X_Point:            Weighted mean location value of point source
%   Sigma_Point:        Weighted uncertainty of point source

N=size(X,1);
M=size(X,2);
DOF = M*N-M;

if DOF==0
    DOF=1;
end

%Make uncertainty larger due to registration error
Sigma = sqrt(Sigma.^2+repmat(Sigma_Reg.^2,[N,1]));

% MLE of center position:
X_Point = sum(X./Sigma.^2,1)./sum(1./Sigma.^2,1);
Sigma_Point = sqrt(1./sum(1./Sigma.^2,1));

%Likelihood at MLE
L=normpdf(X,repmat(X_Point,[N,1]),Sigma);

%Likelihood at Null
L0 = normpdf(X,X,Sigma);

%Calculate likelihood ratio:
R=-2*sum(sum(log(L./L0)));

X2_CDF=inline('gammainc(x/2,k/2)','k','x');
Pvalue=1-X2_CDF(DOF,R);

end

% =============================================================================

function plot_All_Localizations(pvalue, SRT, file, frame)
% Plot localizations, colored by acceptance or reason for rejection.
%
% Inputs:
%    SRT     SRtest object
%    file    (OPTIONAL) file number  to be restricted to
%    frame   (OPTIONAL) frame number to be restricted to

   Params = SRT.Params;
   if nargin == 2
      indx = 1 : length(Params(:, 1));
   else
      indx = find(Params(:, 6) == file & Params(:, 5) == frame);
   end
   k = 0;
   if length(indx) > 0
      k = k + 1;
      legends{k} = 'accepted';
   end

   indx_pvalue  = indx(find(SRT.MinPValue > pvalue(indx)));
   if length(indx_pvalue) > 0
      k = k + 1;
      legends{k} = 'PValue';
      indx = setdiff(indx, indx_pvalue);
   end
   indx_photons = indx(find(Params(indx, 3)   <= SRT.MinPhotons));
   if length(indx_photons) > 0
      k = k + 1;
      legends{k} = 'Photons';
      indx = setdiff(indx, indx_photons);
   end
   indx_xsigma  = indx(find(SRT.CRLB_STD(indx, 1) >= SRT.MaxCRLBSigma));
   if length(indx_xsigma) > 0
      k = k + 1;
      legends{k} = 'xSigma';
      indx = setdiff(indx, indx_xsigma);
   end
   indx_ysigma  = indx(find(SRT.CRLB_STD(indx, 2) >= SRT.MaxCRLBSigma));
   if length(indx_ysigma) > 0
      k = k + 1;
      legends{k} = 'ySigma';
      indx = setdiff(indx, indx_ysigma);
   end
   indx_bg      = indx(find(Params(indx, 4)   >= SRT.MaxBg));
   if length(indx_bg) > 0
      k = k + 1;
      legends{k} = 'BG';
      indx = setdiff(indx, indx_bg);
   end

   MS = {'MarkerSize', 20};

   figure();
   hold on
   plot(Params(indx, 1),         Params(indx, 2),         'k+');
   plot(Params(indx_pvalue, 1),  Params(indx_pvalue, 2),  'r.', MS{:});
   plot(Params(indx_photons, 1), Params(indx_photons, 2), 'b.', MS{:});
   plot(Params(indx_xsigma, 1),  Params(indx_xsigma, 2),  'g.', MS{:});
   plot(Params(indx_ysigma, 1),  Params(indx_ysigma, 2),  'g.', MS{:});
   plot(Params(indx_bg, 1),      Params(indx_bg, 2),      'm.', MS{:});
   if nargin > 2
      title(sprintf('file = %d, frame = %d', file, frame + 1));
   end
   xlabel('x');
   ylabel('y');
   legend(legends, 'Location', 'Best');

   if (0)
      sigma = SRT.MaxCRLBSigma; % * SRT.SRImageZoom;
      n = 4;
      delta = (0 : n) / n * 2 * pi; 
      x = sigma * cos(delta);
      y = sigma * sin(delta);
      for i = 1 : length(indx)
         plot(Params(indx(i), 1) + x, Params(indx(i), 2) + y, 'k-'); 
      end
   end
   hold off
   
end

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
