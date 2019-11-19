classdef Clustering < handle

% Clustering class written by Michael Wester, Keith Lidke, Carolyn Pehlke, Flor
%    Espinoza Hidalgo, Stanly Steinberg and others as noted internally
%    (8/9/2018) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2018 by Michael J. Wester and Keith A. Lidke
%
% DBSCAN_Daszykowski is the recommended DBSCAN algorithm as it is both fast and
% stable under coordinate reordering.
%
% DBSCAN_Daszykowski_noE computes its own value for E.
%
% DBSCAN_Kovesi is both slow and NOT stable under coordinate reordering.
%
% DBSCAN_Pehlke is somewhat slow and NOT stable under coordinate reordering.
%
% DBSCAN_Tran is supposed to be stable under coordinate reordering, but in
% practice, this does not seem to be the case!
%
% Getis is Carolyn Pehlke's Getis statistic based algorithm.
%
% Hierarchal is Matlab's hierarchal clustering algorithm with some small
% additions.
%
% Voronoi is Florian Levet et al's Voronoi based algorithm.
%
% Example main program:
%
%    XY = load('../../analysis/methods/data/pts.csv');  % nm
%    E = 30;   % nm
%    minPts = 3;
%
%    c = Clustering();
%    algorithm = 'Hierarchal';
%    [nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
%    fprintf('number of clusters = %d\n', nC);
%    results = c.clusterStats(XY, C, centers)
%    clusterFig = c.plotClusters(XY, C, centers, ptsI, algorithm);
%    %showm(clusterFig);
%
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

% =============================================================================
properties
% =============================================================================

   % Inputs to cluster.
   Algorithm = []; % input clustering algorithm(s) to perform
   XY = [];        % input point coordinates
   E = [];         % input cluster epsilon
   MinPts = [];    % input minimum points to make a cluster
   % Outputs from cluster.
   NC = [];        % output number of clusters
   C = [];         % output indices of the points in each cluster
   Centers = [];   % output cluster centers
   PtsI = [];      % output isolated points

   % Used by Getis and pair_{auto/cross}correlation.
   % If ROI is provided, it will be used, otherwise if ROI_size is defined,
   % then the (x_min, y_min) computed from the data will be used, otherwise if
   % neither is provided, then the xy_size will be calculated using the
   % (x_min, y_min, x_max, y_max) computed from the data as
   %    xy_size = max(x_max - x_min, y_max - y_min)
   ROI      = [];   % [x_min, y_min, x_size, y_size] (nm)
   ROI_size = [];   % Set x_size = y_size = ROI_size (nm)

   % Properties used by voronoi_Levet.
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

   PlotSkip = [];        % Plots to skip for plot*Hists ('f', 'p', 'c', 'C').

   % Properties below are used by cluster_stats.
   Font_props = {'FontSize', 15, 'FontWeight', 'bold'};
   Line_props = {'LineWidth', 3};
   Fig_ext = 'png';
   Results = '.';   % Directory to store results.
   Xlim = [];   % x-axis limits if defined
   Ylim = [];   % y-axis limits if defined

   Rate          =  20;   % Sampling rate for statistical functions.
   Dendro_cutoff =  50;   % Cluster cutoff for Dendrogram analysis (nm).
   Ripley_cutoff = 200;   % Ripley distance cutoff (nm).
   Confidence    = 2.576; % Bivariate Ripley confidence interval.
                          % 95% confidence interval: Confidence = 1.96
                          % 99% confidence interval: Confidence = 2.576
   Nsims         = 100;   % Bivariate Ripley simulations to run.

% =============================================================================
end % properties

methods
% =============================================================================

function [nC, C, centers, ptsI] = cluster(obj, algorithm, XY, E, minPts)
% Main interface to clustering algorithms.
%
% Inputs:
%    algorithm   one of the clustering algorithms below:
%                DBSCAN_Daszykowski or DBSCAN
%                DBSCAN_Daszykowski_noE (computes its own value for E)
%                DBSCAN_Kovesi
%                DBSCAN_Pehlke
%                DBSCAN_Tran
%                Getis      (Carolyn Pehlke's Getis statistic based algorithm)
%                Hierarchal (Matlab's hierarchal clustering algorithm)
%                Voronoi    (Florian Levet et al's Voronoi based algorithm)
%    XY          point coordinates [N x 2]
%    E           epsilon or cutoff distance
%    minPts      minimum number of points allowed in a cluster
%
% Outputs:
%    nC          number of clusters found
%    C           cell array of XY indices forming each cluster [nC x 1]
%    centers     coordinates of the center of each cluster [nC x n_dim]
%    ptsI        indices of points not found in any cluster

   % Save inputs in object.
   obj.Algorithm = algorithm;
   obj.XY        = XY;
   obj.E         = E;
   obj.MinPts    = minPts;

   if size(XY, 1) == 0
      warning('No points to cluster!');
      return
   end

   n_dim = size(XY, 2);

   switch algorithm
   case {'DBSCAN_Daszykowski', 'DBSCAN'}

      [ptsC, Ctype] = Clustering.dbscan_Daszykowski(XY, minPts - 1, E);
      nC = max([0, ptsC]);
      C = cell(1, nC);
      centers = zeros(n_dim, nC);
      for j = 1 : nC
         C{j} = find(ptsC == j);
         centers(:, j) = mean(XY(C{j}, :), 1);
      end
      ptsI = find(ptsC <= 0);
      %isolated = XY(ptsI, :);

   case 'DBSCAN_Daszykowski_noE'

      [ptsC, Ctype, E] = Clustering.dbscan_Daszykowski(XY, minPts - 1, []);
      obj.E = E;
      nC = max([0, ptsC]);
      C = cell(1, nC);
      centers = zeros(n_dim, nC);
      for j = 1 : nC
         C{j} = find(ptsC == j);
         centers(:, j) = mean(XY(C{j}, :), 1);
      end
      ptsI = find(ptsC <= 0);
      %isolated = XY(ptsI, :);

   case 'DBSCAN_Kovesi'

      [C, ptsC, centers] = Clustering.dbscan_Kovesi(XY', E, minPts + 3);
      nC = max(ptsC);
      ptsI = find(ptsC <= 0)';
      %isolated = XY(ptsI, :);

   case 'DBSCAN_Pehlke'

      [ptsC, nC, centers, N] = Clustering.dbscan_Pehlke(XY, E, minPts);
      centers = centers';
      C = cell(1, nC);
      for j = 1 : nC
         C{j} = find(ptsC == j)';
      end
      ptsI = find(ptsC <= 0)';
      %isolated = XY(ptsI, :);

   case 'DBSCAN_Tran'

      [ptsC, Clabscore] = Clustering.dbscan_Tran(XY, E, minPts);
      nC = max([0; ptsC]);
      C = cell(1, nC);
      centers = zeros(n_dim, nC);
      for j = 1 : nC
         C{j} = find(ptsC == j)';
         centers(:, j) = mean(XY(C{j}, :), 1);
      end
      ptsI = find(ptsC <= 0)';
      %isolated = XY(ptsI, :);

   case 'Getis'

      if n_dim == 3
         error('Getis only works in 2D!');
      end

      if ~isempty(obj.ROI)
         ROI = obj.ROI;
      else
         x_min = min(XY(:, 1));
         y_min = min(XY(:, 2));
         if ~isempty(obj.ROI_size)
            ROI = [x_min, y_min, obj.ROI_size, obj.ROI_size];
         else
            x_max = max(XY(:, 1));
            y_max = max(XY(:, 2));
            xy_size = max(x_max - x_min, y_max - y_min);
            ROI = [x_min, y_min, xy_size, xy_size];
         end
      end

      Gplotting = false;
      [inner_clusters, outer_clusters] = ...
         Clustering.getis(XY, ROI, minPts, Gplotting);
      Cxy = inner_clusters;
      nC = length(Cxy);
      C = cell(1, nC);
      centers = zeros(2, nC);
      for i = 1 : nC
         C{i} = arrayfun(@(x, y) find(XY(:, 1) == x & XY(:, 2) == y), ...
                         Cxy{i}(:, 1), Cxy{i}(:, 2))';
         centers(:, i) = mean(Cxy{i}, 1);
      end
      if nC > 0
         isolated = setdiff(XY, vertcat(Cxy{:}), 'rows');
      else
         isolated = XY;
      end
      ptsI = arrayfun(@(x, y) find(XY(:, 1) == x & XY(:, 2) == y), ...
                      isolated(:, 1), isolated(:, 2))';
      ptsI = sort(ptsI);
      %isolated = XY(ptsI, :);

   case {'Hierarchical', 'Hierarchal'}

      [C, ptsI] = Clustering.hierarchal(XY, E, minPts);
      nC = length(C);
      centers = zeros(n_dim, nC);
      for i = 1 : nC
         centers(:, i) = mean(XY(C{i}, :), 1);
      end
      %isolated = XY(ptsI, :);

   case 'Voronoi'

      %[area, rho, nC, C] = ...
      %   voronoi_Levet(XY, alpha, epsilon, minPts, algorithm);
      [~, ~, nC_a, C_a] = obj.voronoi_Levet(XY, obj.Alpha, E, minPts, ...
                                            obj.Valgorithm);
      nC = nC_a{obj.Valgorithm};
      C  = C_a{obj.Valgorithm};
      centers = zeros(n_dim, nC);
      for i = 1 : nC
         centers(:, i) = mean(XY(C{i}, :), 1);
      end
      ptsI = setdiff(1 : size(XY, 1), horzcat(C{:}));

   otherwise

      error('Unknown algorithm: %s\n', algorithm);

   end

   % Save outputs in object.
   obj.NC      = nC;
   obj.C       = C;
   obj.Centers = centers;
   obj.PtsI    = ptsI;

end

% -----------------------------------------------------------------------------

function cluster_stats(obj, algorithm, P, base_name, particle_types, ...
                            H_nm, V_nm, D_nm)
% Produce various cluster statistics (Hopkin's, Ripley's, etc.) and figures.
%
% Inputs:
%    algorithm        one of the statistical clustering algorithms below
%    P                cell array of point coordinates (n x 2 per cell element)
%                     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%                     where P{1} is (n1 x 2), P{2} is (n2 x 2)
%    base_name        name to identify saved plots, which will have various
%                     descriptive suffixes attached 
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};
%    H_nm             horizontal size of the ROI (nm)
%    V_nm             vertical   size of the ROI (nm)
%    D_nm             depth      size of the ROI (nm) for 3D stats [optional]

   if length(P) == 0
      error('No points to consider!');
   end

   if ~exist('D_nm', 'var')
      D_nm = -1;
   end

   switch algorithm
   case 'PairwiseDist'
      obj.pairwiseDist(P, base_name, particle_types);

   case 'PairwiseMutualDist'
      obj.pairwiseMutualDist(P, base_name, particle_types);

   case 'Hopkins'
      obj.hopkins(P, H_nm, V_nm, D_nm, base_name, particle_types);

   case 'Ripley'
      obj.ripley(P, H_nm, V_nm, D_nm, base_name, particle_types);

   case 'BivariateRipley'
      obj.bivripley(P, H_nm, V_nm, D_nm, base_name, particle_types);

   case 'Dendrogram' 
      obj.dendro_sep(P, base_name, particle_types);

   otherwise
      error('Unknown algorithm: %s\n', algorithm);
   end

end

% =============================================================================

function P = plotCombined(y, bin_width, base_name, x_label, legend_labels, ...
                          x_abbrev, property)
% Plot combined frequency, CDF, PDF and plotSpread plots of the arrays y.
%
% Input:
%    y           cell array of data arrays (need not be the same length)
%    bin_width   bin width for histogram plots
%    base_name   name to identify saved plots, which can have various
%                descriptive suffixes attached 
%    x_label     text for the x-label
%    legend_labels   {legend_label1, legend_label2}
%    x_abbrev    abbreviated identifier for plots used in constructing the
%                filename
%    property    various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Line_props      [{'LineWidth', 3}]
%       Fig_ext         ['png']  figure extension
%       PlotSkip        []       plots that can be skipped:
%                                   'f'   frequency
%                                   'p'   PDF
%                                   'c'   CDF
%                                   'C'   CDF (alternatve)
%                                   's'   PlotSpread
%       Results         ['.']    directory to store results
%       Xlim            []       x-axis limits if defined
% Output:
%    P           P-value for KS test on y1 and y2

   n = numel(y);
   base_text = regexprep(base_name, '_', '\\_');
   colors = ['b', 'r', 'g', 'k', 'c', 'm'];

   % frequency
   if ~any(property.PlotSkip == 'f')
      figure;
      %axes(property.Font_props{:})
      hold on
      for i = 1 : n
         h(i) = histogram(y{i});
         if ~isempty(bin_width)
            h(i).BinWidth = bin_width;
         end
      end
      if ~isempty(property.Xlim)
         xlim(property.Xlim);
      end
      title(base_text);
      xlabel(x_label);
      ylabel('frequency');
      if ~isempty(legend_labels)
         legend(legend_labels, 'Location', 'best');
      end
      hold off
      name = fullfile(property.Results, [base_name, x_abbrev, '_freq']);
      if ~isempty(property.Fig_ext)
         saveas(gcf, name, property.Fig_ext)
      end
      saveas(gcf, name, 'fig');
      close
   end

   % PDF
   if ~any(property.PlotSkip == 'p')
      figure;
      %axes(property.Font_props{:})
      hold on
      for i = 1 : n
         h(i) = histogram(y{i});
         h(i).Normalization = 'probability';
         if ~isempty(bin_width)
            h(i).BinWidth = bin_width;
         end
      end
      if ~isempty(property.Xlim)
         xlim(property.Xlim);
      end
      title(base_text);
      xlabel(x_label);
      ylabel('PDF');
      if ~isempty(legend_labels)
         legend(legend_labels, 'Location', 'best');
      end
      hold off
      name = fullfile(property.Results, [base_name, x_abbrev, '_PDF']);
      if ~isempty(property.Fig_ext)
         saveas(gcf, name, property.Fig_ext)
      end
      saveas(gcf, name, 'fig');
      close
   end

   % CDF
   if ~any(property.PlotSkip == 'c')
      figure;
      %axes(property.Font_props{:})
      hold on
      for i = 1 : n
         h(i) = histogram(y{i});
         h(i).Normalization = 'cdf';
         if ~isempty(bin_width)
            h(i).BinWidth = bin_width;
         end
      end
      BinLimits(1) = min([h.BinLimits]);
      BinLimits(2) = max([h.BinLimits]);
      for i = 1 : n
         h(i).BinLimits = BinLimits;
      end
      if ~isempty(property.Xlim)
         xlim(property.Xlim);
      end
      title(base_text);
      xlabel(x_label);
      ylabel('CDF');
      if ~isempty(legend_labels)
         legend(legend_labels, 'Location', 'best');
      end
      hold off
      %X1 = h1.BinEdges;
      %Y1 = h1.Values;
      name = fullfile(property.Results, [base_name, x_abbrev, '_CDF']);
      if ~isempty(property.Fig_ext)
         saveas(gcf, name, property.Fig_ext)
      end
      saveas(gcf, name, 'fig');
      close
   end

   % CDF (alternative)
   if ~any(property.PlotSkip == 'C')
      isnull = all(cellfun(@isempty, y));
      figure;
      %axes(property.Font_props{:})
      hold on
      if ~isempty(property.Xlim)
         x_max = property.Xlim(2);
      elseif ~isnull
         x_max = max(cellfun(@max, y));
      else
         x_max = 1;
      end
      if ~isnull
         T = cell(1, n);
         X = cell(1, n);
         Y = cell(1, n);
         for i = 1 : n
            T{i} = tabulate(y{i});
            X{i} = [0; T{i}(:, 1); x_max];
            Y{i} = [0; cumsum(T{i}(:, 2)) / numel(y{i}); 1];
            plot(X{i}, Y{i}, [colors(i), '-'], property.Line_props{:});
         end
      end
      if ~isempty(property.Xlim)
         xlim(property.Xlim);
      end
      title(base_text);
      xlabel(x_label);
      ylabel('CDF');
      if ~isempty(legend_labels)
         legend(legend_labels, 'Location', 'best');
      end
      hold off
      name = fullfile(property.Results, [base_name, x_abbrev, '_CDF2']);
      if ~isempty(property.Fig_ext)
         saveas(gcf, name, property.Fig_ext)
      end
      saveas(gcf, name, 'fig');
      close
   end

   % PlotSpread
   if ~any(property.PlotSkip == 's')
      figure;
      %axes(property.Font_props{:})
      if ~isempty(property.Xlim)
         x_max = property.Xlim(2);
      elseif ~all(cellfun(@isempty, y))
         x_max = max(cellfun(@max, y));
      else
         x_max = 1;
      end
      plotSpread(y, 'showMM', 1, 'xNames', legend_labels);
      hold on
      if ~isempty(property.Xlim)
         ylim(property.Xlim);
      end
      title(base_text);
      ylabel(x_label);
      hold off
      name = fullfile(property.Results, [base_name, x_abbrev, '_PS']);
      if ~isempty(property.Fig_ext)
         saveas(gcf, name, property.Fig_ext)
      end
      saveas(gcf, name, 'fig');
      close
   end

   P = zeros(n);
   for i = 1 : n
      P(i, i) = 1;
      for j = i+1 : n
         [~, P(i, j)] = kstest2(y{i}, y{j});
         P(j, i) = P(i, j);
      end
   end

end

% -----------------------------------------------------------------------------

function pairwiseDist(obj, P, base_name, particle_types)
% Plot pairwise distances and CDFs.
%
% Inputs:
%    P                cell array of point coordinates (n x 2 per cell element)
%                     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%                     where P{1} is (n1 x 2), P{2} is (n2 x 2)
%    base_name        name to identify saved plots, which will have various
%                     descriptive suffixes attached 
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};

   bins = 100;
   sims = 20;

   for m = 1 : length(particle_types)
      if ~isempty(obj.Fig_ext)
         figure('Visible', 'off');
      else
         figure;
      end
      axes(obj.Font_props{:})
      hold on

      X = P{m};
      D = pdist(X);
      M = mean(D);
      S = std(D);

      % Establish the positions of the bin centers (x) based on a normal random
      % distribution that inherits its characteristics from the coordinates in
      % X given the specified number of bins (also, see comment below).
      Xr = M + S*randn(size(X));
      Dr = pdist(Xr);
      [~, x] = hist(Dr, bins);

      y = hist(D, x);
      plot(x, y, 'k-', obj.Line_props{:});

      % Create a random distribution with the same characteristics as X:
      % based on the same number of points, mean, standard deviation.  Average
      % over the specified number of simulations (sims).
      YR = 0;
      for i = 1 : sims
         Xr = M + S*randn(size(X));
         Dr = pdist(Xr);
         yr = hist(Dr, x);
         YR = YR + yr;
      end
      yr = YR / sims;
      plot(x, yr, 'r--', obj.Line_props{:});

      legend('data', 'random', 'Location', 'NorthEast');
      title(['Pairwise Distance PDF for ', base_name, '\_', ...
             particle_types{m}]);
      xlabel('distance (nm)');
      ylabel('frequency');
      hold off
      name = fullfile(obj.Results, ...
                      [base_name, '_', particle_types{m}, '_pairwisePDF']);
      if ~isempty(obj.Fig_ext)
         print(['-d', obj.Fig_ext], name);
      else
         saveas(gcf, name);
         delete(gcf);
      end

      if ~isempty(obj.Fig_ext)
         figure('Visible', 'off');
      else
         figure;
      end
      axes(obj.Font_props{:})
      hold on

      [f, xx] = ecdf(sort(D));
      plot(xx, f, 'k-', obj.Line_props{:});

      %f = cumsum(1 : numel(D));
      %f = f / f(end);
      %xx = sort(D);
      %plot(xx, f, 'b-', obj.Line_props{:});

      [fr, xxr] = ecdf(sort(Dr));
      plot(xxr, fr, 'r--', obj.Line_props{:});

      %fr = cumsum(1 : numel(Dr));
      %fr = fr / fr(end);
      %xxr = sort(Dr);
      %plot(xxr, fr, 'm--', obj.Line_props{:});

      legend('data', 'random', 'Location', 'SouthEast');
      title(['Pairwise Distance CDF for ', base_name, '\_', ...
             particle_types{m}]);
      xlabel('distance (nm)');
      ylabel('frequency');
      hold off
      name = fullfile(obj.Results, ...
                      [base_name, '_', particle_types{m}, '_pairwiseCDF']);
      if ~isempty(obj.Fig_ext)
         print(['-d', obj.Fig_ext], name);
      else
         saveas(gcf, name);
         delete(gcf);
      end
   end

end

% -----------------------------------------------------------------------------

function pairwiseMutualDist(obj, P, base_name, particle_types)
% Plot pairwise distances and CDFs between two populations of particles.
%
% Inputs:
%    P                cell array of point coordinates (n x 2 per cell element)
%                     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%                     where P{1} is (n1 x 2), P{2} is (n2 x 2)
%    base_name        name to identify saved plots, which will have various
%                     descriptive suffixes attached 
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};

   bins = 100;
   sims = 20;

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on

   if length(particle_types) ~= 2
      error('Mutual distances can only be computed between 2 particle types!');
   end

   X = P{1};
   Y = P{2};
   D = pdist2(X, Y);
   D = D(:);
%  M = mean(D);
%  S = std(D);

   % Establish the positions of the bin centers (x) based on a normal random
   % distribution that inherits its characteristics from the coordinates in
   % X given the specified number of bins (also, see comment below).
%  Xr = M + S*randn(size(X));
%  Yr = M + S*randn(size(Y));
%  Dr = pdist2(Xr, Yr);
%  [~, x] = hist(Dr, bins);

%  y = hist(D, x);
   [y, x] = hist(D, bins);
   plot(x, y, 'k-', obj.Line_props{:});

   % Create a random distribution with the same characteristics as X:
   % based on the same number of points, mean, standard deviation.  Average
   % over the specified number of simulations (sims).
%  YR = 0;
%  for i = 1 : sims
%     Xr = M + S*randn(size(X));
%     Dr = pdist(Xr);
%     yr = hist(Dr, x);
%     YR = YR + yr;
%  end
%  yr = YR / sims;
%  plot(x, yr, 'r--', obj.Line_props{:});

%  legend('data', 'random', 'Location', 'NorthEast');
   title(['Pairwise Distance PDF for ', base_name, '\_', ...
          particle_types{1}, ',', particle_types{2}]);
   xlabel('distance (nm)');
   ylabel('frequency');
   hold off
   name = fullfile(obj.Results, [base_name, '_', particle_types{1}, ',', ...
                                 particle_types{2}, '_pairwisePDF2']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on

   [f, xx] = ecdf(sort(D));
   plot(xx, f, 'k-', obj.Line_props{:});

   %f = cumsum(1 : numel(D));
   %f = f / f(end);
   %xx = sort(D);
   %plot(xx, f, 'b-', obj.Line_props{:});

%  [fr, xxr] = ecdf(sort(Dr));
%  plot(xxr, fr, 'r--', obj.Line_props{:});

   %fr = cumsum(1 : numel(Dr));
   %fr = fr / fr(end);
   %xxr = sort(Dr);
   %plot(xxr, fr, 'm--', obj.Line_props{:});

%  legend('data', 'random', 'Location', 'SouthEast');
   title(['Pairwise Distance CDF for ', base_name, '\_', ...
          particle_types{1}, ',', particle_types{2}]);
   xlabel('distance (nm)');
   ylabel('frequency');
   hold off
   name = fullfile(obj.Results, [base_name, '_', particle_types{1}, ',', ...
                                 particle_types{2}, '_pairwiseCDF2']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

end

% -----------------------------------------------------------------------------

function H = hopkins_ROIcombined(obj, base_name, n_ROIs, ROIs)
% Use Hopkins' statistics to test the clustering of a series of ROIs.
%
% Originally written by Michael Wester and Stanly Steinberg in 2008; extended
% to 3D and combined ROIs in 2017--2018.
%
% Inputs:
%    base_name   name to identify saved plots, which can have various
%                descriptive suffixes attached 
%    n_ROIs      number of ROIs to combine
%    ROIs        n_ROIs cell array containing the following fields:
%       XY1      nX sets of (x, y {, z}) coordinates
%       ROI      ROI limits in the form
%                   [x_min, x_max, y_min, y_max {, z_min, z_max}]
% Output:
%    H           mean Hopkin's statistic (over ntests) for each ROI

base_text = regexprep(base_name, '_', '\\_');

test = 5;        % The number of test points.
ntests = 1000;   % The number of Hopkins statistics to use.

% Dimension (2D or 3D)
dim = 2;
if numel(size(ROIs{1}.XY1)) == 3
   dim = 3;
end

% Compute ntests Hopkins' statistics (using test probes) for each ROI.
H  = zeros(1, n_ROIs);
HS = zeros(1, ntests);
for i = 1 : n_ROIs
   X = ROIs{i}.XY1;
   ROI = ROIs{i}.ROI;
   x_min = ROI(1);
   x_max = ROI(2);
   y_min = ROI(3);
   y_max = ROI(4);
   X(:, 1) = X(:, 1) - x_min;
   X(:, 2) = X(:, 2) - y_min;
   if dim == 2
      for j = 1 : ntests
         HS(j) = Clustering.hopkinstat(X, x_max - x_min, y_max - y_min, test);
      end
   else
      z_min = ROI(5);
      z_max = ROI(6);
      X(:, 3) = X(:, 3) - z_min;
      for j = 1 : ntests
         HS(j) = Clustering.hopkinstat3(X, x_max - x_min, y_max - y_min, ...
                                           z_max - z_min, test);
      end
   end
   % H(i) will be the mean of the ntests Hopkins' statistics for ROI i.
   H(i) = mean(HS);
end

% Plot a histogram of the findings.
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
histogram(H, 25);
xlim([0, 1]);
if ~isempty(obj.Xlim)
    xlim(obj.Xlim);
end
title([base_text, ' [all ROIs]']);
xlabel('H (Hopkins statistic)');
ylabel('frequency');
hold off
name = fullfile(obj.Results, [base_name, '_hopkins_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot a PDF of the findings.
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
% For some reason, Normalization = PDF does not work correctly, but probability
% does, so ...
%histogram(H, 25, 'Normalization', 'PDF');
histogram(H, 25, 'Normalization', 'probability');
xlim([0, 1]);
ylim([0, 1]);
if ~isempty(obj.Xlim)
    xlim(obj.Xlim);
end
if ~isempty(obj.Ylim)
    ylim(obj.Ylim);
end
title([base_text, ' [all ROIs]']);
xlabel('H (Hopkins statistic)');
ylabel('probability');
hold off
name = fullfile(obj.Results, [base_name, '_hopkins_PDF_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

end

% -----------------------------------------------------------------------------

function hopkins(obj, P, H_nm, V_nm, D_nm, base_name, particle_types)
% Use the Hopkins' statistic to test the clustering of the points in P.
%
% Written by Michael Wester and Stanly Steinberg in 2008.
%
% Inputs:
%    P                cell array of point coordinates (n x 2 per cell element)
%                     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%                     where P{1} is (n1 x 2), P{2} is (n2 x 2)
%    H_nm             horizontal size of the ROI (nm)
%    V_nm             vertical   size of the ROI (nm)
%    D_nm             depth      size of the ROI (nm) for 3D stats [optional]
%    base_name        name to identify saved plots, which will have various
%                     descriptive suffixes attached 
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};

test = 5;        % The number of test points.
ntests = 1000;   % The number of Hopkins statistics to use.
bins = 100;      % The number of bins for the probability graphs.
fitting = false; % Gaussian fitting of the data.

% Compute the Hopkins' statistic for each negative.
fprintf('Compute the Hopkins statistic for %s.\n', base_name);
for m = 1 : length(particle_types)
   X = [];
   if D_nm <= 0   % 2D
      for i = 1 : ntests
         X = [X, Clustering.hopkinstat(P{m}, H_nm, V_nm, test)];
      end
   else   % 3D
      for i = 1 : ntests
         X = [X, Clustering.hopkinstat3(P{m}, H_nm, V_nm, D_nm, test)];
      end
   end
   H{m} = X;
end

% The analytic distribution for the Hopkins test:
xa = linspace(0, 1, obj.Rate);
ya = ((xa.^(test-1)).*((1-xa).^(test-1)))/(gamma(test)^2/gamma(2*test));

% Plot the summary statistics.
% The PDF for each negative in the experiment.
fprintf('   Plot the Hopkins statistics.\n');
for m = 1 : length(particle_types)
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on
   [X, V] = Clustering.histogram(H{m}, [0, 1], bins);
   Pr = V*bins/ntests;
   bar(X, Pr, 1)
   hold on
   if fitting
      % Gaussian fit to the data.
      f = fit(X', Pr', 'gauss2');
      [fit_max, i_max] = max(f(X));
      x_max = X(i_max);
      h = plot(f, 'g');
      set(h, obj.Line_props{1}, obj.Line_props{2});
   end
   plot(xa, ya, 'r', obj.Line_props{:})
   axis([0 1 0 ceil(max(max(Pr), max(ya)))]);
   xlabel('H -- The Hopkins Statistic');
   ylabel('PDF');
   if fitting
      title({['Hopkins PDF for ', base_name, '\_', particle_types{m}], ...
             sprintf('(x, fit)_{max} = (%5.3f, %.3f)', x_max, fit_max)});
   else
      title(['Hopkins PDF for ', base_name, '\_', particle_types{m}]);
   end
   if fitting
      legend('data', 'Gaussian fit', 'random', 'Location', 'NorthWest');
   else
      legend('data', 'random', 'Location', 'NorthWest');
   end
   name = fullfile(obj.Results, ...
                   [base_name, '_', particle_types{m}, '_hopkinspdf']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end
end

hold off

end

% -----------------------------------------------------------------------------

function ripley_ROIcombined(obj, base_name, n_ROIs, ROIs)
% Use Ripley's statistics to test the clustering of a series of ROIs all of the
% same size.
%
% Originally written by Michael Wester and Stanly Steinberg in 2008; extended
% to 3D and combined ROIs in 2017--2018.
%
% Inputs:
%    base_name   name to identify saved plots, which can have various
%                descriptive suffixes attached 
%    n_ROIs      number of ROIs to combine
%    ROIs        n_ROIs cell array containing the following fields:
%       XY1      nX sets of (x, y {, z}) coordinates
%       ROI      ROI limits in the form
%                   [x_min, x_max, y_min, y_max {, z_min, z_max}]

% Dimension (2D or 3D)
dim = 2;
if numel(size(ROIs{1}.XY1)) == 3
   dim = 3;
end

% Divide the interval [0, cutoff] into Rate number of pieces of length dr.
% The look at the distance d between two points and add 1 to R(d/dr).
% This provides a PDF for the distances.
dr=obj.Ripley_cutoff/obj.Rate;
K = zeros(1, obj.Rate + 1);
for m = 1 : n_ROIs
   x_min = ROIs{m}.ROI(1);
   x_max = ROIs{m}.ROI(2);
   y_min = ROIs{m}.ROI(3);
   y_max = ROIs{m}.ROI(4);
   if dim == 3
      z_min = ROIs{m}.ROI(5);
      z_max = ROIs{m}.ROI(6);
   end

   X = ROIs{m}.XY1;
   nX = length(X);
   R = zeros(1, obj.Rate);
   if dim == 2
      for i = 1 : nX
         for j = 1 : i - 1
            p = ceil(sqrt((X(j,1) - X(i,1))^2 + (X(j,2) - X(i,2))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   else
      for i = 1 : nX
         for j = 1 : i - 1
            p = ceil(sqrt((X(j,1) - X(i,1))^2 + (X(j,2) - X(i,2))^2 + ...
                          (X(j,3) - X(i,3))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   end
   % Add a zero to the beginning of R for nicer plotting.
   R = [0, R];
   % Compensate for the j=1:i-1 savings.
   R = 2*R;
   % Convert the PDF to a CDF.
   for i = 2 : length(R)
      R(i) = R(i) + R(i-1);
   end
   % Average over the number of particles.
   R = R/nX;
   % Compute the intensity and normalize R.
   if dim == 2
      lambda = nX/((x_max - x_min)*(y_max - y_min));
   else
      lambda = nX/((x_max - x_min)*(y_max - y_min)*(z_max - z_min));
   end
   R = R/lambda;
   % Accumulate over each ROI.
   K = K + R;
end
% Normalize over all ROIs.
K = K / n_ROIs;

% Plot Ripley
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
% The analytic expected value.
%plot([0:obj.Rate]*dr, pi*([0:obj.Rate]*dr).^2, 'k', obj.Line_props{:})
A = pi*([0:obj.Rate]*dr).^2;
plot([0:obj.Rate]*dr, A, 'k', obj.Line_props{:})
plot([0:obj.Rate]*dr, K, 'r', obj.Line_props{:})
legend({'random', 'probes'}, 'Location', 'NorthWest');
axis([0 obj.Ripley_cutoff 0 max(K)]);
xlabel('r (nm)');
ylabel('K(r)');
title(['Ripley K [RC] for ', base_name]);
hold off
name = fullfile(obj.Results, [base_name, '_ripleyK_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot sqrt of Ripley
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
grid on
M = sqrt(A/pi) - [0 : obj.Rate]*dr;
N = sqrt(K/pi) - [0 : obj.Rate]*dr;
plot([0:obj.Rate]*dr, M, 'k', obj.Line_props{:})
plot([0:obj.Rate]*dr, N, 'r', obj.Line_props{:})
legend({'random', 'probes'}, 'Location', 'Best');
xlabel('r (nm)');
ylabel('L(r) - r ');
title(['Ripley L(r) - r [RC] for ', base_name]);
hold off
name = fullfile(obj.Results, [base_name, '_ripleyL_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot Ripley/(pi*r^2) (in nm)
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
grid on
% Avoid dividing by 0.
M = A(2:length(A));
M = M./(pi*(dr*[1:obj.Rate]).^2);
N = K(2:length(K));
N = N./(pi*(dr*[1:obj.Rate]).^2);
plot([0:obj.Rate]*dr, [1, M], 'k', obj.Line_props{:});
plot([0:obj.Rate]*dr, [1, N], 'r', obj.Line_props{:});
legend({'random', 'probes'}, 'Location', 'Best');
axis([0 obj.Ripley_cutoff 0 max(N)]);
xlabel('r (nm)');
ylabel('K(r)/(\pi r^2)');
title(['Normalized Ripley [RC] for ', base_name]);
hold off
name = fullfile(obj.Results, [base_name, '_ripleyNormalized_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

end

% -----------------------------------------------------------------------------

function ripley(obj, P, H_nm, V_nm, D_nm, base_name, particle_types)
% Use Ripley's statistics to test the clustering of the points in P.
%
% Written by Michael Wester and Stanly Steinberg in 2008.
%
% Inputs:
%    P                cell array of point coordinates (n x 2 per cell element)
%                     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%                     where P{1} is (n1 x 2), P{2} is (n2 x 2)
%    H_nm             horizontal size of the ROI (nm)
%    V_nm             vertical   size of the ROI (nm)
%    D_nm             depth      size of the ROI (nm) for 3D stats [optional]
%    base_name        name to identify saved plots, which will have various
%                     descriptive suffixes attached 
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};

M = length(particle_types);
for i = 1 : M
   %probes{i} = [particle_types{i}, ' nm'];
   probes{i} = particle_types{i};
end

colors = {'r', 'g', 'b'};

% Dimension (2D or 3D)
dim = 2;
if D_nm > 0
   dim = 3;
end

% Compute the Ripley statistcs.
fprintf('Compute Ripley statistic for %s.\n', base_name);

% Divide the interval [0, cutoff] into Rate number of pieces of length dr.
% The look at the distance d between two points and add 1 to R(d/dr).
% This provides a PDF for the distances.
dr=obj.Ripley_cutoff/obj.Rate;
for m = 1 : length(particle_types)
   X = P{m};
   nX = length(X);
   R = zeros(1, obj.Rate);
   if dim == 2
      %for i = 1 : length(X)
      for i = 1 : nX
         for j = 1 : i - 1
            p = ceil(sqrt((X(j,1) - X(i,1))^2 + (X(j,2) - X(i,2))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   else
      for i = 1 : nX
         for j = 1 : i - 1
            p = ceil(sqrt((X(j,1) - X(i,1))^2 + (X(j,2) - X(i,2))^2 + ...
                          (X(j,3) - X(i,3))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   end
   % Add a zero to the beginning of R for nicer plotting.
   R = [0, R];
   % Compensate for the j=1:i-1 savings.
   R = 2*R;
   % Convert the PDF to a CDF.
   for i = 2 : length(R)
      R(i) = R(i) + R(i-1);
   end
   % Average over the number of particles.
   %R = R/length(X);
   R = R/nX;
   % Compute the intensity and normalize R.
   if dim == 2
      %lambda = length(X)/(H_nm*V_nm);
      lambda = nX/(H_nm*V_nm);
   else
      lambda = nX/(H_nm*V_nm*D_nm);
   end
   R = R/lambda;
   K(m) = {R};
end

% Plot Ripley
fprintf('   Plot Ripley statistic.\n');
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
% The analytic expected value.
plot([0:obj.Rate]*dr, pi*([0:obj.Rate]*dr).^2, 'k', obj.Line_props{:})
mx = 0;
for m = 1 : length(particle_types)
   N = K{m};
   mx = max(mx, max(N));
   plot([0:obj.Rate]*dr, N, colors{m}, obj.Line_props{:})
end
legend({'random', probes{:}}, 'Location', 'NorthWest');
axis([0 obj.Ripley_cutoff 0 mx]);
xlabel('r (nm)');
ylabel('K(r)');
title(['Ripley for ', base_name]);
hold off
name = fullfile(obj.Results, [base_name, '_ripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot sqrt of Ripley
fprintf('   Plot the L(r) statistic.\n');
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
mx = 0;
for m = 1 : length(particle_types)
   N = sqrt(K{m}/pi) - [0 : obj.Rate]*dr;
%  mx = max(mx, max(N));
   plot([0:obj.Rate]*dr, N, colors{m}, obj.Line_props{:})
end
legend(probes, 'Location', 'NorthWest');
xlabel('r (nm)');
ylabel('L(r) - r ');
title(['Ripley L(r) - r for ', base_name]);
hold off
name = fullfile(obj.Results, [base_name, '_l_r_ripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot Ripley/(pi*r^2) (in microns)
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
mx = 0;
fprintf('   Plot the K(r)/(pi r^2) statistic.\n');
for m = 1 : length(particle_types)
   % Avoid dividing by 0.
   N = K{m}(2:length(K{m}));
   N = N./(pi*(dr*[1:obj.Rate]).^2);
   mx = max(mx, max(N));
   plot([0:obj.Rate]*dr, [0, N], colors{m}, obj.Line_props{:});
end
legend(probes, 'Location', 'NorthEast');
axis([0 obj.Ripley_cutoff 0 mx]);
xlabel('r (nm)');
ylabel('K(r)/(\pi r^2)');
title(['Normalized Ripley for ', base_name]);
hold off
name = fullfile(obj.Results, [base_name, '_norm_ripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% m       is the number of different types of particles.
% n       is the number of images in the experiment.
% P{m}    is a cell array of matrices of coordinates.
% Rate    is the sampling rate.

end

% -----------------------------------------------------------------------------

function bivripley(obj, P, H_nm, V_nm, D_nm, base_name, particle_types)
% Use the bivariate Ripley's statistic to test the clustering of the points in
% P.
%
% Written by Michael Wester and Stanly Steinberg in 2008.
%
% Inputs:
%    P                cell array of point coordinates (n x 2 per cell element)
%                     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%                     where P{1} is (n1 x 2), P{2} is (n2 x 2)
%    H_nm             horizontal size of the ROI (nm)
%    V_nm             vertical   size of the ROI (nm)
%    D_nm             depth      size of the ROI (nm) for 3D stats [optional]
%    base_name        name to identify saved plots, which will have various
%                     descriptive suffixes attached 
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};

% Dimension (2D or 3D)
dim = 2;
if D_nm > 0
   dim = 3;
end

% Compute the bivariate Ripley statistics.
% Add Q for second species
fprintf('Compute the bivariate Ripley statistic for %s.\n', base_name);

% Divide the interval [0, cutoff] into Rate number of pieces of length dr.
% The look at the distance d between two points and add 1 to R(d/dr).
% This provides a PDF for the distances.
dr = obj.Ripley_cutoff/obj.Rate;
X = P{1};
Y = P{2};
nX = length(X);
nY = length(Y);
R = zeros(1, obj.Rate);
if dim == 2
   %for i = 1 : length(X)
   %   for j = 1 : length(Y)
   for i = 1 : nX
      for j = 1 : nY
         p = ceil(sqrt((X(i,1) - Y(j,1))^2 + (X(i,2) - Y(j,2))^2)/dr);
         if p > 0 & p <= obj.Rate
            R(p) = R(p) + 1;
         end
      end
   end
else
   for i = 1 : nX
      for j = 1 : nY
         p = ceil(sqrt((X(i,1) - Y(j,1))^2 + (X(i,2) - Y(j,2))^2 + ...
                       (X(i,3) - Y(j,3))^2)/dr);
         if p > 0 & p <= obj.Rate
            R(p) = R(p) + 1;
         end
      end
   end
end

% Add a zero to the beginning of R for nicer plotting.
R = [0, R];
% Convert the PDF to a CDF.
for i = 2 : length(R)
   R(i) = R(i) + R(i-1);
end
   
% Compute the intensity and normalize R.
% This needs a fix. xxx
% xxx should use the area of the crop_area (err crop region).
% A = H_nm*V_nm
if dim == 2
   %lambda1 = length(X)/(H_nm*V_nm);
   %lambda2 = length(Y)/(H_nm*V_nm);
   lambda1 = nX/(H_nm*V_nm);
   lambda2 = nY/(H_nm*V_nm);
   R = R/(lambda1*lambda2*(H_nm*V_nm));
else
   lambda1 = nX/(H_nm*V_nm*D_nm);
   lambda2 = nY/(H_nm*V_nm*D_nm);
   R = R/(lambda1*lambda2*(H_nm*V_nm*D_nm));
end
K(1, 2) = {R};

% Plot sqrt of Ripley
fprintf('   Plotting the L(r) statistic for %s.\n', base_name);
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
N = sqrt(K{1, 2}/pi) - [0:obj.Rate]*dr;
plot([0:obj.Rate]*dr, N, '-r', obj.Line_props{:})
xlabel('Distance (nm)');
ylabel('L(r)-r ');
NN = N;

for s = 1 : obj.Nsims
   % Create random point vectors for simulations
   % and calculate Ripley statistic
   if dim == 2
      %Xr = [rand(length(X),1)*H_nm rand(length(X),1)*V_nm];
      %Yr = [rand(length(Y),1)*H_nm rand(length(Y),1)*V_nm];
      Xr = [rand(nX,1)*H_nm rand(nX,1)*V_nm];
      Yr = [rand(nY,1)*H_nm rand(nY,1)*V_nm];
      Rr = zeros(1, obj.Rate);
      %for i = 1 : length(Xr)
      %   for j = 1 : length(Yr)
      for i = 1 : nX
         for j = 1 : nY
            p=ceil(sqrt((Xr(i,1) - Yr(j,1))^2 + (Xr(i,2) - Yr(j,2))^2)/dr);
            if p > 0 & p <= obj.Rate
               Rr(p) = Rr(p) + 1;
            end
         end
      end
   else
      Xr = [rand(nX,1)*H_nm rand(nX,1)*V_nm rand(nX,1)*D_nm];
      Yr = [rand(nY,1)*H_nm rand(nY,1)*V_nm rand(nY,1)*D_nm];
      Rr = zeros(1, obj.Rate);
      for i = 1 : nX
         for j = 1 : nY
            p=ceil(sqrt((Xr(i,1) - Yr(j,1))^2 + (Xr(i,2) - Yr(j,2) + ...
                         Xr(i,3) - Yr(j,3))^2)/dr);
            if p > 0 & p <= obj.Rate
               Rr(p) = Rr(p) + 1;
            end
         end
      end
   end
   % Add a zero to the beginning of R for nicer plotting.
   Rr = [0, Rr];
   % Convert the PDF to a CDF.
   for i = 2 : length(Rr)
      Rr(i) = Rr(i) + Rr(i-1);
   end
   % Compute the intensity and normalize R.
% Fix area here also xxx.
   if dim == 2
      %lambda1r = length(Xr)/(H_nm*V_nm);
      %lambda2r = length(Yr)/(H_nm*V_nm);
      lambda1r = nX/(H_nm*V_nm);
      lambda2r = nY/(H_nm*V_nm);
      Rr = Rr/(lambda1*lambda2*(H_nm*V_nm));
   else
      lambda1r = nX/(H_nm*V_nm*D_nm);
      lambda2r = nY/(H_nm*V_nm*D_nm);
      Rr = Rr/(lambda1*lambda2*(H_nm*V_nm*D_nm));
   end
   N = sqrt(Rr/pi) - [0:obj.Rate]*dr;
   % Plot all random simulations
%  plot([0:obj.Rate]*dr, N, obj.Line_props{:})
   hold on
   % Put Ripley data into a matrix
   E(s, 1:length(Rr))=N;
end

for ss = 1 : length(Rr)
   % Calculate standard deviation of all simulations at particular dr
   S(ss, 1) = std(E(:, ss));
   % Calculate mean of all simulations at particular dr 
   A(ss) = mean(E(:, ss));
   High(ss) = A(ss)+obj.Confidence*S(ss);
   Low(ss)  = A(ss)-obj.Confidence*S(ss);
end

% Plot data
% grid off
% set(gca, 'Box', 'on', 'LineWidth', 3, 'FontSize', 18, 'FontWeight', 'bold')
plot([0:obj.Rate]*dr, High, '--k', obj.Line_props{:});
plot([0:obj.Rate]*dr, Low,  '--k', obj.Line_props{:});
plot([0:obj.Rate]*dr, H_nm, '*k', obj.Line_props{:});

% Fit the figure
b_n = base_name;
p_t = particle_types{1};
for i = 2 : length(particle_types)
   p_t = [p_t ',' particle_types{i}];
end

%axis([0 200 min(Low)-5 max(NN)+5])
axis([0 obj.Ripley_cutoff min(Low)-5 max(NN)+5])
title(['Bivariate Ripley Analysis for ', b_n, ' (', p_t, ')'])
legend('data', 'confidence', 'Location', 'Best')
name = fullfile(obj.Results, [b_n, '_bivripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% m       is the number of different types of particles.
% n       is the number of images in the experiment.
% P{m}    is a cell array of matrices of coordinates.
% Rate    is the sampling rate.

end

% -----------------------------------------------------------------------------

function dendro_sep(obj, P, base_name, particle_types)

% This function plots dendrograms and does clustering analyis.
% Flor Espinoza modified on July 14, 2010.
% The function dendrogram in matlat constructs clusters trees
% if we have lx<=30, the nodes of the dendrogram  equals lx, i.e, Node1= lx(1)
% (node one equals particle 1) but if lx>30 we can choose the number of nodes

%clear all;
% lx: integer, total number of particles
% NN=30: number of nodes, we can define how many nodes we want, i.e. NN=lx
% PN: array, particles in nodes
% nppn: vector, number of particles per node
% MP: array, particles in each cluster
% XM: vector, number of particles per cluster ; XM(1:lx) are singletons
% N1: scalar, total numner of clusters
% NC: vector, cluster numbering for clusters with 2 or more particles
% ncc: vector, cummulative number of clusters in cluster
% CNP: array, cummulative particles in clusters
% TNPC: vector, cummulative num of particles in clusters
% DMNC=d_I: scalar, (intrinsic distance) distance for maximum num of cluster
% R: scalar, cutt-off distance

%#########################################################################
% This results are print in the "analyis_dendro.txt" file!!
%#########################################################################
% results at the intrinsic distance, d_I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MNC: max number of cluster
% DMNC: distance at the MNC, also called intrinsic distance, d_I
% TPC: total number of particles in clusters
% PPC: total percentage of particles in clusters
% MCS(1), MCS(2), MCS(3): three max cluster sizes
% PSC,P24C,PB5C: percent singlets, clusters with 2-4, or more than 5 particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results at cut-off distance, R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCR:  max number of cluster
% dR:  distance (equal o very close to R)
% TPCR: total number of particles in clusters
% PPCR: total percentage of particles in clusters
% MCSR(1),MCSR(2),MCSR(3): three max cluster sizes
% PSCR,P24CR,PB5CR: percent singlets, clusters with 2-4, or more than 5
%    particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z(1,3); min distance between two particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = length(particle_types);
R = obj.Dendro_cutoff;

fprintf('Plotting the dendrogram for %s.\n', base_name);

for m = 1 : M
   X = P{m};
   x=X(1);
   y=X(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   lx=length(X);
   NN=30;
% distance vector i.e, if lx=10; Y(1:9)=dP1:(dP2:dP10),Y(10:17)=dP2:(dP3:dP10)
   Y = pdist(X,'euclidean');
% Z,(lx-1)*3 matrix, Z(:,1:2) contains particles or clusters, Z(:,3)
% is the distance from particles to cluster or from clusters to clusters. The
% first cluster form in Z is numbering lx+1, the second lx+2 and so on
   Z = linkage(Y,'single');

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})

   [H,T] = dendrogram(Z,NN,'colorthreshold','default');
   set(H,'LineWidth',2) % T is a vector that contains the nodes of the particles

   title(['Dendrogram for ', base_name, '\_', particle_types{m}]);
   ylabel('d(nm)')
   hold off

   name = fullfile(obj.Results, ...
                   [base_name, '_',  particle_types{m}, '_dendrogram']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% counting how many particles each node has
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if lx <=30 % each node has one particle
      for i=1:lx
         PN{i}=i;
      end
   else
      for i=1:NN
         k=find(T==i);
         PN{i}=k';
         nppn(i)=length(PN{i});
      end
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particles per cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i=1:lx
       MP{i}=i;
       XM(i)=1;
   end
   for i=1:lx-1
      MP{lx+i}=[MP{Z(i,1)} MP{Z(i,2)}];% identifying particles per cluster
      XM(lx+i)=length(MP{lx+i});
      MCN{i}= T(MP{lx+i})'; % identifying clusters with nodes
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accumulated number of clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   N1=max(max(Z(:,1:2)))+1;
   NC=lx+1:N1;
   DM{1}=NC(1);
   ncc(1)=1; % at distance d(1) we have one cluster
   lk=DM{1};

   for i=2:lx-1
      k=[lk lx+i];
      lk1=[];
      for j=1:length(k)
         if ((k(j)~= Z(i,1))  && (k(j)~= Z(i,2)))
            lk1=[lk1 k(j)];
         end
      end
      lk=[lk1];
      DM{i}=lk1;
      ncc(i)=length(lk1);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acummulated particles in clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i=1:lx-1
      ind=DM{i};
      le=length(ind);
      v=[];
      for j=1:le
         v=[v MP{ind(j)}];
      end
      CNP{i}=v;
      TNPC(i)=length(CNP{i});
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ploting accumulated distance and number of particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on

   
   plot(Z(:,3), TNPC',obj.Line_props{:})
   xlabel('d(nm)'); ylabel('number of particles')

   title(['Total number of particles in clusters for ', base_name,'\_', ...
           particle_types{m} ]);
   hold off

   name = fullfile(obj.Results, ...
                  [base_name, '_', particle_types{m}, '_num_particles']);

   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the maximun number of clusters at intrinsic distance d_I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   MNC=max(ncc);
   ind1=find(ncc==MNC);
   DMNC=round(Z(ind1(1),3)); % min distance where max number of clusters occurs
   PC=XM(DM{ind1(1)}); % identifying clusters and number of particle
   [indd,indpp] = sort(PC,2,'descend'); %finding the biggest clusters
   MCS=indd(1:3); % three max cluster sizes
   TPC=sum(PC); % total number of particles in clusters
   PSC=(lx-TPC)/lx; % percentage of singlestons
   n1=find(2<=PC & PC<=4);
   P24C=(sum(PC(n1)))/lx;
   PB5C=1-(PSC+P24C);
   y1=[PSC P24C PB5C];
   PPC=(TPC/lx); % percentage of  two or more particles in clusters
   M1=[MNC DMNC lx TPC PPC MCS]; %max NC, dist d, num  part, num part in clus,
   % of part in clus, 3 big clus

% DM{ind1(1)}
% PC
% mDm=max(DM{ind1(1)})
% MP{mDm}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the maximun number of clusters at cut_off distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   indR=find(Z(:,3)<=R); LL=length(indR);
   dR=round(Z(indR(LL),3));
   MCR=(ncc(indR(LL)));
   PCR=XM(DM{LL}); % identifying clusters and number of particles
   [inddR,indppR] = sort(PCR,2,'descend'); %finding the biggest clusters
   MCSR=inddR(1:3); % three max cluster sizes
   TPCR=sum(PCR); % total number of particles in clusters
   PSCR=(lx-TPCR)/lx;  % percentage of singlestons
   n2=find(2<=PCR & PCR<=4);
   P24CR=(sum(PCR(n2)))/lx;
   PB5CR=1-(PSCR+P24CR);
   y2=[PSCR P24CR PB5CR];
   PPCR=(TPCR/lx);% percentage of two or more particles in clusters
   M2=[MCR dR lx TPCR PPCR MCSR]; %%max NC, at cut-off distance, # part, # part in clus, 
   % of part in clus,%3 big clus
   
   pt=[particle_types{m} 'nm']; DMNCu =[num2str(DMNC) 'nm'];
   dRu =[num2str(dR) 'nm']; md=[num2str(round(Z(1,3))) 'nm'];

   fid1 = fopen([obj.Results, '/analysis_dendro.txt'], 'a');
   fprintf(fid1,'%-10s  %4s %4d %4d %6s %4d %5.2f %4d %4d %4d %5.2f %5.2f %5.2f %4d %6s %4d %5.2f %4d %4d %4d %5.2f %5.2f %5.2f %5s\n', ...
   base_name, pt, lx, MNC, DMNCu, TPC, PPC, MCS(1), MCS(2),  ...
   MCS(3), PSC, P24C, PB5C, MCR, dRu, TPCR, PPCR, MCSR(1), MCSR(2), MCSR(3), ...
   PSCR, P24CR, PB5CR, md);
   fclose(fid1);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on

   d_I=ones(size(ncc))*DMNC; % intrinsic distance
   plot(Z(:,3),ncc,obj.Line_props{:})
   %hold on
   plot(d_I,ncc,'.',obj.Line_props{:})
   xlabel('distance(nm)'); ylabel('number of clusters');
   %hold off

   title(['Total number of clusters for ', base_name, '\_', ...
          particle_types{m}]);
   hold off

   name = fullfile(obj.Results, ...
                   [base_name, '_', particle_types{m}, '_num_clusters']);

   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   clear X TNPC Z lx ind1 ncc MNC DMNC lx TPC PPC MCS MCR dR lx TPCR PPCR ...
         MCSR PSC P24C PB5C PSCR P24CR PB5CR
end % m

end

% =============================================================================

function [area, rho, nC, C] = voronoi_Levet(obj, xy, alpha, epsilon, ...
                                            minPts, algorithm, ShrinkFactor)
% Florian Levet, Eric Hosy, Adel Kechkar, Corey Butler, Anne Beghin, Daniel
% Choquet and Jean-Baptiste Sibarita, ``SR-Tesseler: a method to segment and
% quantify localization-based super-resolution microscopy data'', _Nature
% Methods_, 2015 (DOI:10.1038/NMETH.3579).
%
% Inputs:
%    xy          n x 2 array of point coordinates
%    alpha       ratio of local density / overall density for a point's Voronoi
%                region to be considered sufficiently dense for clustering
%                purposes
%    epsilon     distance constraint for clustering; no constraint if <= 0
%    minPts      minimum number of points needed to form a cluster
%    algorithm   cell array of algorithms to apply:
%                1   [0] calculations consider Voronoi regions only
%                2   [1] calculations consider Voronoi regions and their
%                    adjacent neighbors
%                3   [1M] consider the median density of each cell and its
%                    neighbors
% Outputs:
%    area   cell array of Voronoi cell areas per point
%    rho    cell array of Voronoi cell densities per point
%    nC     cell array of number of clusters discovered
%    C      cell array of indices of points in each cluster

   if any(algorithm < 1 | algorithm > 3)
      error('Invalid algorithm selected: %d', algorithm);
   end

   if ~exist('ShrinkFactor', 'var')
      ShrinkFactor = 0.5;
   end

   XY = double(xy);
   dim = size(XY, 2);

   X = XY(:, 1);
   Y = XY(:, 2);
   n_XY = length(X);
   if dim == 2
      area_XY = (max(X) - min(X)) * (max(Y) - min(Y));
   else
      Z = XY(:, 3);
      area_XY = (max(X) - min(X)) * (max(Y) - min(Y)) * (max(Z) - min(Z));
   end
   rho_XY = n_XY / area_XY;

   [v, c] = voronoin(XY);
   n_v = size(v, 1);
   n_c = length(c);

   % Find what Voronoi cells (c) contain each vertex (v).
   v2c = cell(1, n_v);  
   for i = 1 : n_c
      c_i = c{i};
      for j = 1 : length(c_i)
         v2c{c_i(j)} = [v2c{c_i(j)}, i];
      end
   end

   % For each cell, collect together the indices for itself and all its
   % neighbors.
   self_nbrs = cell(1, n_c);
   for i = 1 : n_c
      c_i = c{i};
      for j = 1 : length(c_i)
         k = v2c{c_i(j)};
         self_nbrs{i} = [self_nbrs{i}, k];
      end
      self_nbrs{i} = unique(self_nbrs{i});
   end

   % Level 0: consider each cell individually.
   area_0 = zeros(1, n_XY);
   rho_0  = zeros(1, n_XY);
   for i = 1 : n_c
      c_i = c{i};
      if all(c_i ~= 1)
         if dim == 2
            area_0(i) = polyarea(v(c_i, 1), v(c_i, 2));
         else
            [~, area_0(i)] = boundary(v(c_i, 1), v(c_i, 2), v(c_i, 3), ...
                                      ShrinkFactor);
         end
         rho_0(i) = 1 / area_0(i);
      else
         area_0(i) = Inf;
         rho_0(i) = 0;
      end
   end

   if any(algorithm == 1)
      [~, i_rho_0] = find(rho_0 >= alpha*rho_XY);
      [nC_0, C_0] = Clustering.cluster_voronoi(i_rho_0, self_nbrs, epsilon, ...
                                               minPts, XY);

      if obj.Plotting
         if dim == 2
            Clustering.plot_voronoi(X, Y, v, c, rho_0 ./ rho_XY, ...
                                    'rho_0 / rho_a', i_rho_0, obj.PtIDs);
         else
            Clustering.plot_voronoi3(X, Y, Z, v, c, rho_0 ./ rho_XY, ...
                                     'rho_0 / rho_a', i_rho_0, obj.PtIDs);
         end
      end

      area{1} = area_0;
      rho{1}  = rho_0;
      nC{1}   = nC_0;
      C{1}    = C_0;
   end

   if any(algorithm == 2)
      % Level 1: consider each cell and its neighbors.
      area_1 = zeros(1, n_XY);
      rho_1  = zeros(1, n_XY);
      for i = 1 : n_c
         c_i = c{i};
         if all(c_i ~= 1)
            % Cells with vertex 1, which is the point at Infinity, are
            % excluded.
            n_self_nbrs = length(self_nbrs{i});
            for j = 1 : n_self_nbrs
               area_1(i) = area_1(i) + area_0(self_nbrs{i}(j));
            end
            rho_1(i) = n_self_nbrs / area_1(i);
         else
            area_1(i) = Inf;
            rho_1(i) = 0;
         end
      end

      [~, i_rho_1] = find(rho_1 >= alpha*rho_XY);
      [nC_1, C_1] = Clustering.cluster_voronoi(i_rho_1, self_nbrs, epsilon, ...
                                               minPts, XY);

      if obj.Plotting
         if dim == 2
            Clustering.plot_voronoi(X, Y, v, c, rho_1 ./ rho_XY, ...
                                    'rho_1 / rho_a', i_rho_1, obj.PtIDs);
         else
            Clustering.plot_voronoi3(X, Y, Z, v, c, rho_1 ./ rho_XY, ...
                                     'rho_1 / rho_a', i_rho_1, obj.PtIDs);
         end
      end

      area{2} = area_1;
      rho{2}  = rho_1;
      nC{2}   = nC_1;
      C{2}    = C_1;

      %area_1A = zeros(1, n_XY);
      %rho_1A  = zeros(1, n_XY);
      %for i = 1 : n_c
      %   c_i = c{i};
      %   if all(c_i ~= 1)
      %      n_self_nbrs = 0;
      %      for j = 1 : n_c
      %         if ~isempty(intersect(c{i}, c{j}))
      %            n_self_nbrs = n_self_nbrs + 1;
      %            area_1A(i) = area_1A(i) + area_0(j);
      %         end
      %      end
      %      rho_1A(i) = n_self_nbrs / area_1A(i);
      %   else
      %      area_1A(i) = Inf;
      %      rho_1A(i) = 0;
      %   end
      %end
      %plot_voronoi(X, Y, v, c, rho_1A);
   end

   if any(algorithm == 3)
      % Level 1M: consider the median density of each cell and its neighbors.
      area_1M = area_0;
      rho_1M  = zeros(1, n_XY);
      for i = 1 : n_c
         rho_1M(i) = median(rho_0(self_nbrs{i}));
      end

      [~, i_rho_1M] = find(rho_1M >= alpha*rho_XY);
      [nC_1M, C_1M] = Clustering.cluster_voronoi(i_rho_1M, self_nbrs, ...
                                                 epsilon, minPts, XY);

      if obj.Plotting
         if dim == 2
            Clustering.plot_voronoi(X, Y, v, c, rho_1M ./ rho_XY, ...
                                    'rho_{1M} / rho_a', i_rho_1M, obj.PtIDs);
         else
            Clustering.plot_voronoi3(X, Y, Z, v, c, rho_1M ./ rho_XY, ...
                                     'rho_{1M} / rho_a', i_rho_1M, obj.PtIDs);
         end
      end

      area{3} = area_1M;
      rho{3}  = rho_1M;
      nC{3}   = nC_1M;
      C{3}    = C_1M;
   end

end

% -----------------------------------------------------------------------------

function nn_dists = nn_ROIcombined(obj, base_name, n_ROIs, ROIs)
% Plot the mean particle nearest neighbor distances for a series of ROIs.
%
% Inputs:
%    base_name   name to identify saved plots, which can have various
%                descriptive suffixes attached 
%    n_ROIs      number of ROIs to combine
%    ROIs        n_ROIs cell array containing the following fields:
%       XY1      nX sets of (x, y {, z}) coordinates
%       ROI      ROI limits in the form
%                   [x_min, x_max, y_min, y_max {, z_min, z_max}]
% Output:
%    nn          nearest nighbor distances over all particles for each ROI

   nn_dists = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      nn_dists{i} = Clustering.nn_distances(ROIs{i}.XY1);
   end

   % Histogram
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on
   histogram(arrayfun(@(i) mean(nn_dists{i}), 1 : n_ROIs), 25);
   if ~isempty(obj.Xlim)
       xlim(obj.Xlim);
   end
   title(base_name);
   xlabel('nearest neighbor distance (nm)');
   ylabel('frequency');
   hold off
   name = fullfile(obj.Results, [base_name, '_nn_RC']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   % PDF
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on
   histogram(arrayfun(@(i) mean(nn_dists{i}), 1 : n_ROIs), 25, ...
             'Normalization', 'probability');
   ylim([0, 1]);
   if ~isempty(obj.Xlim)
       xlim(obj.Xlim);
   end
   if ~isempty(obj.Ylim)
       ylim(obj.Ylim);
   end
   title([base_name, ' [all ROIs]']);
   xlabel('nearest neighbor distance (nm)');
   ylabel('probability');
   hold off
   name = fullfile(obj.Results, [base_name, '_nn_PDF_RC']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

end

% =============================================================================
end % methods

methods(Static)
% =============================================================================

function [class,type,Eps]=dbscan_Daszykowski(x,k,Eps)
% Written by Michal Daszykowski
% -------------------------------------------------------------------------
% Function: [class,type]=dbscan(x,k,Eps)
% -------------------------------------------------------------------------
% Aim: 
% Clustering the data with Density-Based Scan Algorithm with Noise (DBSCAN)
% -------------------------------------------------------------------------
% Input: 
% x - data set (m,n); m-objects, n-variables
% k - number of objects in a neighborhood of an object 
% (minimal number of objects considered as a cluster)
% Eps - neighborhood radius, if not known avoid this parameter or put []
% -------------------------------------------------------------------------
% Output: 
% class - vector specifying assignment of the i-th object to certain 
% cluster (m,1)
% type - vector specifying type of the i-th object 
% (core: 1, border: 0, outlier: -1)
% -------------------------------------------------------------------------
% Example of use:
% x=[randn(30,2)*.4;randn(40,2)*.5+ones(40,1)*[4 4]];
% [class,type]=dbscan(x,5,[]);
% -------------------------------------------------------------------------
% References:
% [1] M. Ester, H. Kriegel, J. Sander, X. Xu, A density-based algorithm for 
% discovering clusters in large spatial databases with noise, proc. 
% 2nd Int. Conf. on Knowledge Discovery and Data Mining, Portland, OR, 1996, 
% p. 226, available from: 
% www.dbs.informatik.uni-muenchen.de/cgi-bin/papers?query=--CO
% [2] M. Daszykowski, B. Walczak, D. L. Massart, Looking for 
% Natural Patterns in Data. Part 1: Density Based Approach, 
% Chemom. Intell. Lab. Syst. 56 (2001) 83-92 
% -------------------------------------------------------------------------
% Written by Michal Daszykowski
% Department of Chemometrics, Institute of Chemistry, 
% The University of Silesia
% December 2004
% http://www.chemometria.us.edu.pl

[m,n]=size(x);

if nargin<3 | isempty(Eps)
   [Eps]=Clustering.epsilon(x,k);
end

x=[[1:m]' x];
[m,n]=size(x);
type=zeros(1,m);
no=1;
touched=zeros(m,1);

for i=1:m
    if touched(i)==0;
       ob=x(i,:);
       D=Clustering.dist_Daszykowski(ob(2:n),x(:,2:n));
       ind=find(D<=Eps);
    
       if length(ind)>1 & length(ind)<k+1       
          type(i)=0;
          class(i)=0;
       end
       if length(ind)==1
          type(i)=-1;
          class(i)=-1;  
          touched(i)=1;
       end

       if length(ind)>=k+1; 
          type(i)=1;
          class(ind)=ones(length(ind),1)*max(no);
          
          while ~isempty(ind)
                ob=x(ind(1),:);
                touched(ind(1))=1;
                ind(1)=[];
                D=Clustering.dist_Daszykowski(ob(2:n),x(:,2:n));
                i1=find(D<=Eps);
     
                if length(i1)>1
                   class(i1)=no;
                   if length(i1)>=k+1;
                      type(ob(1))=1;
                   else
                      type(ob(1))=0;
                   end

                   for i=1:length(i1)
                       if touched(i1(i))==0
                          touched(i1(i))=1;
                          ind=[ind i1(i)];   
                          class(i1(i))=no;
                       end                    
                   end
                end
          end
          no=no+1; 
       end
   end
end

i1=find(class==0);
class(i1)=-1;
type(i1)=-1;

end

%...........................................
function [Eps]=epsilon(x,k)

% Function: [Eps]=epsilon(x,k)
%
% Aim: 
% Analytical way of estimating neighborhood radius for DBSCAN
%
% Input: 
% x - data matrix (m,n); m-objects, n-variables
% k - number of objects in a neighborhood of an object
% (minimal number of objects considered as a cluster)



[m,n]=size(x);

Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);

end

%............................................
function [D]=dist_Daszykowski(i,x)

% function: [D]=dist(i,x)
%
% Aim: 
% Calculates the Euclidean distances between the i-th object and all objects in x	 
%								    
% Input: 
% i - an object (1,n)
% x - data matrix (m,n); m-objects, n-variables	    
%                                                                 
% Output: 
% D - Euclidean distance (m,1)



[m,n]=size(x);
D=sqrt(sum((((ones(m,1)*i)-x).^2)'));

if n==1
   D=abs((ones(m,1)*i-x))';
end

end

% =============================================================================

function [C, ptsC, centres] = dbscan_Kovesi(P, E, minPts)
% DBSCAN clustering algorithm
%
% Usage:  [C, ptsC, centres] = dbscan(P, E, minPts)
%
% Arguments:
%         P - dim x Npts array of points.
%         E - Distance threshold.
%    minPts - Minimum number of points required to form a cluster.
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.
%      ptsC - Array of length Npts listing the cluster number associated with
%             each point.  If a point is denoted as noise (not enough nearby
%             elements to form a cluster) its cluster number is 0.
%   centres - dim x Nc array of the average centre of each cluster.

% Reference:
% Martin Ester, Hans-Peter Kriegel, Jörg Sander, Xiaowei Xu (1996). "A
% density-based algorithm for discovering clusters in large spatial databases
% with noise".  Proceedings of the Second International Conference on Knowledge
% Discovery and Data Mining (KDD-96). AAAI Press. pp. 226-231.  
% Also see: http://en.wikipedia.org/wiki/DBSCAN

% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% PK January 2013
    
    [dim, Npts] = size(P);
    
    ptsC  = zeros(Npts,1);
    C     = {};
    Nc    = 0;               % Cluster counter.
    Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.
    
    for n = 1:Npts
       if ~Pvisit(n)                            % If this point not visited yet
           Pvisit(n) = 1;                       % mark as visited
           neighbourPts = Clustering.regionQuery(P, n, E); % and find its neighbours

           if length(neighbourPts) < minPts-1  % Not enough points to form a cluster
%          if length(neighbourPts) < minPts    % Not enough points to form a cluster
               ptsC(n) = 0;                    % Mark point n as noise.
           
           else                % Form a cluster...
               Nc = Nc + 1;    % Increment number of clusters and process
                               % neighbourhood.
           
               C{Nc} = [n];    % Initialise cluster Nc with point n
               ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.
               
               ind = 1;        % Initialise index into neighbourPts array.
               
               % For each point P' in neighbourPts ...
               while ind <= length(neighbourPts)
                   
                   nb = neighbourPts(ind);
                   
                   if ~Pvisit(nb)        % If this neighbour has not been visited
                       Pvisit(nb) = 1;   % mark it as visited.
                       
                       % Find the neighbours of this neighbour and if it has
                       % enough neighbours add them to the neighbourPts list
                       neighbourPtsP = Clustering.regionQuery(P, nb, E);
                       if length(neighbourPtsP) >= minPts
                           neighbourPts = [neighbourPts  neighbourPtsP];
                       end
                   end            
                   
                   % If this neighbour nb not yet a member of any cluster add it
                   % to this cluster.
                   if ~ptsC(nb)  
                       C{Nc} = [C{Nc} nb];
                       ptsC(nb) = Nc;
                   end
                   
                   ind = ind + 1;  % Increment neighbour point index and process
                                   % next neighbour
               end
           end
       end
    end
    
    % Find centres of each cluster
    centres = zeros(dim,length(C));
    for n = 1:length(C)
        for k = 1:length(C{n})
            centres(:,n) = centres(:,n) + P(:,C{n}(k));
        end
        centres(:,n) = centres(:,n)/length(C{n});
    end

end % of dbscan    
    
%------------------------------------------------------------------------

function neighbours = regionQuery(P, n, E)
% Find indices of all points within distance E of point with index n
% This function could make use of a precomputed distance table to avoid
% repeated distance calculations, however this would require N^2 storage.
% Not a big problem either way if the number of points being clustered is
% small.   For large datasets this function will need to be optimised.

% Arguments:
%              P - the dim x Npts array of data points
%              n - Index of point of interest
%              E - Distance threshold
    
    E2 = E^2;   
    [dim, Npts] = size(P);
    neighbours = [];
    
    for i = 1:Npts
        if i ~= n
            % Test if distance^2 < E^2 
            v = P(:,i)-P(:,n);
            dist2 = v'*v;
            if dist2 < E2 
               neighbours = [neighbours i];     
            end
        end
    end
    
end % of regionQuery

% =============================================================================

function [C Nc Centers N] = dbscan_Pehlke(X, epsilon, minpts)
% DBSCAN algorithm included by Carolyn Pehlke in SuperCluster.

   V=zeros(size(X,1),1);   %Visited
   C=zeros(size(X,1),1);   %Cluster ID, 0 means noise

   cid=0;

   for ii=1:size(X,1)
       if V(ii); continue; end
       V(ii)=1;
       [nnpts n]=Clustering.getnnpts(X,ii,epsilon);
       if n<minpts;  C(ii)=0; continue; end
       cid=cid+1;
       C(nnpts)=ones(n,1)*cid; %claim all into cluster
       [V C]=Clustering.exapandcluster(X,nnpts,cid,epsilon,minpts,V,C);  %follow points
    end

    for nn=1:cid
        Centers(nn,:)=mean(X(C==nn,:),1);
        N(nn)=sum(C==nn);
    end

    Nc=cid;
end

function [nnpts n]=getnnpts(X,ii,epsilon)
   d=Clustering.dist_Pehlke(X(ii,:),X);
   nnpts=find(d<epsilon);
   n=length(nnpts);
end

function [V C]=exapandcluster(X,nnpts,cid,epsilon,minpts,V,C)

   %nnpts is the set of un-evaluated points

   while ~isempty(nnpts) %still points left to explore neighbors
      ii=nnpts(1);    %take index from top of list
      nnpts(1)=[];    %remove ii from list
      V(ii)=1;        %mark ii as visited

      [nnpts_p n_p]=Clustering.getnnpts(X,ii,epsilon);   %find points close to point ii
      for jj=1:n_p    %test each point for more neighbors
         if V(nnpts_p(jj))==0    %if point isn't claimed
            V(nnpts_p(jj))=1;   %mark as vistited
            C(nnpts_p(jj))=cid; %claim to current cluster
            nnpts=cat(2,nnpts,nnpts_p(jj)); %add point to list to find it's neighbors  
         end
      end
   end
end

% distance function - get distance from point in question to all other
% points
function distlist = dist_Pehlke(point,mat)
   x1 = point(1);
   y1 = point(2);
   for rr = 1:size(mat,1)
      x2 = mat(rr,1);
      y2 = mat(rr,2);
      distlist(rr) = sqrt((x1-x2)^2+(y1-y2)^2);
   end
end

% =============================================================================

function [labs labscore] = dbscan_Tran(a,Eps,MinPts)
%  DBSCAN clustering
% [labs labscore] = dbscan(A,Eps,MinPts)
%
% DBSCAN clustering of data matrix in A. labels is a vector with 
% cluster labels for each vector.
% 
% In case of publication of any application of this method,
% please, cite the original work:
% Thanh N. Tran*, Klaudia Drab, Michal Daszykowski, "Revised DBSCAN algorithm 
% to cluster data with dense adjacent clusters", Chemometrics and Intelligent 
% Laboratory Systems, 120:9296.
% DOI: 10.1016/j.chemolab.2012.11.006 


UNCLASSIFIED = 0;
BORDER = -2;
    % Square Eps in order not to square all distances of points
    %eps = eps^2;
    
     m = size(a,1);
     labs = zeros(m,1);
     ClusterId = 1;
     for i=1:m
        if labs(i) == UNCLASSIFIED                 
            
           % Expand cluster ClusterId
           % Get a set of points of distance < eps
           [ExpandClusterReturn labs]= Clustering.expandCluster(a,labs,i,ClusterId,Eps, MinPts);
           if ExpandClusterReturn
               ClusterId = ClusterId +1;
           end
        end
     end

     % Step 3:
     labscore = labs; core_index = find(labscore > 0);
     border_points = find(labs==BORDER);
     % For xborder in border_list but has no ClusterId
     for i=1:length(border_points)
            % xborder in border_list
            currentB = border_points(i);
            d=Clustering.distance(+a(currentB,:),+a(core_index,:),1);
            % the closest core-points 
            [tmp nearest_core]=min(d);
            nearest_core_index=core_index(nearest_core);
            %Assign xborder to ClusterId of the closest core-points 
            labs(currentB)=labs(nearest_core_index);
     end
end

function [ExpandClusterReturn labs]= expandCluster(a,labs,i,ClusterId, Eps, MinPts)
UNCLASSIFIED = 0;
NOISE = -1;
BORDER = -2;
           % calculate distances
           d=Clustering.distance(+a(i,:),+a(:,:),1);
           
           % seeds = Retrieve_Neighbors(xi, Eps) 
           seeds = find(d < Eps);
           
           % If |seeds | < MinPts
           if size(seeds,2) < MinPts, 
               labs(i) = NOISE;                 % Assign xi as noise
               ExpandClusterReturn = 0;         % Return without expansion success 
           else
               	% STEP 1: xi is a identified as a starting core-point for ClusterId
                %labs(i) = ClusterId;  %        % Thanh changed in May 2012
                % Exclude xi from the Seeds
                %seeds = setdiff(seeds, i);     % Thanh changed in May 2012
                % STEP 2:  Identify chains
                % For all xj in seeds
                while ~isempty(seeds) % Not an empty seeds
                    % current point is the first point in seeds

                    currentP = seeds(1);
                    d=Clustering.distance(+a(currentP,:),+a(:,:),1);
                    % NEps(xj)= Retrieve_Neighbors(xj,Eps)
                    result = find(d <= Eps);
                    %If | NEps(xj) | >= MinPts                              // xj is a core point
                    if length(result) >= MinPts 
                        % Assign xj  to ClusterId
                        labs(currentP) = ClusterId;
                        % Add all UNCLASSIFIED in NEps(xj) to seeds
                        result_unclassified = result(find(labs(result)==UNCLASSIFIED));
                        result_noise = result(find(labs(result)==NOISE));
                        % Temperary complete the intermediate chain ... the chain can be extended in the comming steps.
                        % The border points have not yet been assigned to any cluster
                        labs([result_unclassified result_noise]) = BORDER;   % first assign to a border-point... later will be reassigned to e.g. core-point
                        seeds = union(seeds,result_unclassified);
                    end
                  % Exclude the current point in seeds and go back to the loop
                  seeds = seeds(2:size(seeds,2));
                end % end while
                % Return with expansion success 
                ExpandClusterReturn = 1;    % return true 
           end
end

function d = distance(X,Y, dosqrt)
% function d = distance(X,Y, dosqrt);
% dosqrt: default 0: No square root to reduce computation time (only for comparation perpose)
% Euclidean distance matrix between row vectors in X and Y

if nargin < 3, dosqrt = 0; end

U=~isnan(Y); Y(~U)=0;
V=~isnan(X); X(~V)=0;
if ~dosqrt d=abs(X.^2*U'+V*Y'.^2-2*X*Y');
else
    d=sqrt(X.^2*U'+V*Y'.^2-2*X*Y');
end

end

% =============================================================================

function [inner_clusters, outer_clusters] = getis(xy, ROI, cutoff, plotting)
% Getis-based clustering developed by Carolyn Pehlke.

   regs(1).points = xy;
   regs(1).ROI = ROI;
   regs(1).area = ROI(3) * ROI(4);
   regs(1).name = 'ROI';

   ar = 10;   % area of single subregion

   regs = Clustering.plotClusts(regs, ar, cutoff, plotting);

   inner_clusters = regs.clusts{1};
   outer_clusters = regs.clusts{2};
   n_inner_clusters = length(inner_clusters);
   n_outer_clusters = length(outer_clusters);
   for i = 1 : n_inner_clusters
      inner_clusters{i} = inner_clusters{i}(:, [2, 1]);
   end
   for i = 1 : n_outer_clusters
      outer_clusters{i} = outer_clusters{i}(:, [2, 1]);
   end
end

% callback for cluster quantification
%   function plotClusts(src,eventdata)
    function regs = plotClusts(regs, ar, cutoff, plotting)
%       set(handles.workTxt,'Visible','on') % set the "working" banner so people don't click like crazy while matlab is thinking
%       drawnow % and do it right now
%       regs = get(handles.pickButt,'UserData'); % get the ROI data
        if isempty(regs) % if you haven't picked ROIs, do nothing
            return;
        end
        num = length(regs); % get the number of ROIs
%       ar = str2num(get(handles.arBox,'String')); % area of single subregion
%       cutoff = str2num(get(handles.minCtBox,'String')); % get the minimum cluster size from GUI
%       set([handles.clustTab,handles.getisTab],'Enable','on') % turn on the cluster and cluster stats tabs
%       clustShow % flip to the cluster view panel
%       set(handles.clustTab,'Value',1) % and set the tab value appropriately
        % lay groundwork for figuring out proper axes later on
        tmpMax = 0;
        tmpMin = 10000;
        tmpCon = 0;
        conMin = 10000;
        tmpDist = 0;
        distMin = 10000;
        passval = zeros(size(regs));
        for tt = 1:num
            pts = double(regs(tt).points);
            ntests = 1000;
            %test = 5;
            test = min(5, size(pts, 1));
            if size(pts, 1) == 1
               regs.clusts = cell(1, 2);
               passval(tt) = 1;
               continue;   % No clustering possible.
            end
            bins = 100;
            xx = [];
            xrange = ceil(max(regs(tt).points(:,1)) - min(regs(tt).points(:,1))); % get the proper x and y range for the ROI
            yrange = ceil(max(regs(tt).points(:,2)) - min(regs(tt).points(:,2)));
            for kk = 1 : ntests
                xx = [xx, Clustering.sR_hopkinstat(regs(tt).points, xrange, yrange, test)];
            end
            H{tt} = xx;
            if mean(H{tt}) > .3 && mean(H{tt}) < .75
%               disp(['Data is random, unable to find clusters. Skipping region ',num2str(tt)])
                regs.clusts = cell(1, 2);
                passval(tt) = 1;
            end
        end

        indx = find(~passval);
        % declare some variables ahead of time
        Ggroup = cell(1,sum(~passval));
        Gmat = cell(1,num);
        clusters = cell(1,sum(~passval));
        hull = cell(1,sum(~passval));
        sig = cell(1,sum(~passval));

        leg = [{'Inner'},{'Outer'}];
        % for each ROI
        if any(~passval)
        for aa = 1:length(indx)
            try
            rect = regs(indx(aa)).ROI; % get the size of the ROI
            pts = double(regs(indx(aa)).points);
            %%%%% This could be faster if we do it outside the loop and
            %%%%% assume that people are being sensible and making their
            %%%%% ROIs all the same size. But are they?
            poss = 0:5:regs(indx(aa)).ROI(3)/2; % get a range of distances half the size of the ROI
            [K dr] = Clustering.sR_Ripley(regs(indx(aa)).ROI(3)/2,regs,length(poss)-1); % find the Ripley's K function over that range
            % estimate maximum G cutoff based on peak of Ripley's L
            % function
            L = sqrt(K{aa}/pi) - poss; % get L(r)-r
            L = diff(L);
            I = find(L<=0,1,'first'); % and find the first critical point
            maxcut = ceil(ceil(poss(I+1)/2)/10)*10; % get the max cutoff, rounded to the nearest 10nm
            gcut = 0:ar:maxcut; % make the range of cutoff values
            cutval = ceil(gcut/ar); % getis cutoffs in units of subregions
            cutval = cutval(2:end);
            % get the average distance (the square root of the inverse
            % density)
            minper = round(sqrt((rect(3)*rect(4))/size(pts,1))/ar);
            % get the two max lengths - the smaller of the two for the
            % inner clusters, the larger of the two for the outer clusters
            threshval = [min(minper,max(cutval)),max(max(cutval),minper)];
            % get the proper R for making the alpha hull.
            if 2*minper*ar < 70
                R = 9*minper*ar;
            else
                R = 2*minper*ar;
            end
            % divide ROI into equal-area subregions or pixels
            % assign points to proper subregion to make a histogram image
            % of the ROI
            xmin = rect(1); xmax = rect(1) + rect(3); ymin = rect(2); ymax = rect(2) + rect(4);
            imszX = ceil((xmax-xmin)/ar)+1;  % image size x
            imszY = ceil((ymax-ymin)/ar)+1; % image size y
            im = zeros(imszX,imszY); % make empty image
            locations = cell(imszX,imszY); % make a cell for saving the SR localizations that are in each pixel
            x = round((pts(:,1)-xmin)/ar)+1; % get the x and y points ready for binning
            y = round((pts(:,2)-ymin)/ar)+1;
            for ii = 1:size(x,1) % bin localizations into pixels
                im(y(ii),x(ii))= im(y(ii),x(ii))+1;
                locations{y(ii),x(ii)} = [locations{y(ii),x(ii)};[pts(ii,2),pts(ii,1)]]; % and save them in the localization variable
            end
            [row col] = find(im); % find all positive pixels in im
            [Ggroup{aa} Gmat{aa}] = Clustering.calcGetis(im,cutval,row,col); % run getis
            maxes = [];
            getg = Gmat{aa}(:,:,1);
            maxes = imregionalmax(getg,4);
% % %             bins = unique(sort(getg)); % get all G values greater than 0, sort and make bins
% % %             n = hist(getg,bins); % make histogram to get counts
% % %             tmpn = sum(n)*.75; % get level that is 75% of the total
            % find bin that is closest to that level, that is threshold
% % %             for vv = 1:length(bins)
% % %                 if (sum(n(vv:end)) - tmpn) < 0
% % %                    tval = bins(vv-1);
% % %                    break
% % %                 end
% % %             end
            % the G values greater than the threshold are the seed points
            % for clusters
% % %             maxes =  getg(:,:) >= tval;
            [r c] = find(maxes);
            % these are to be sure that outer clusters only combine
            % existing inner clusters rather than adding single points to
            % them
            newRow = [];
            newCol = [];
            dt = diff(Ggroup{aa},1,2); % get derivative G vs d curves

            for zz = 1:2 % we want to make two passes of this to get inner and outer clusters
             if zz == 2 % if this is the second (outer) pass, limit points to just points contained in existing clusters
                 row = newRow;
                 col = newCol;
             end
            % find the intersection of seedpoints out of total points
            [C ia ib] = intersect([row col],[r c],'rows');
            whichclust = zeros(size(row)); % set cluster definition variable
            for ww = 1:size(ia,1) % for all seed points
                ind = find(~dt(ia(ww),:),1,'first'); % find first zero in derivative
                minthreshval = [min(ind,minper),max(ind,minper)]; % use either ind or the smaller of minper/ind as search radius
                locmax = [row(ia(ww)),col(ia(ww))]; % find current seed point
                tmpdist = Clustering.dist_Pehlke(locmax,[row col]); % find distances from seed point to other points in region
                if ~isempty(ind) % if there is zero in first derivative
                    thresh = minthreshval(zz); % use smaller of ind,minper for first round and ind for second round
                else % if ind is empty (ie if no zero in derivative)
                    thresh = threshval(zz); % use either minper or cutval depending on if inner or outer clusters
                end
                newind = find(tmpdist < ceil(thresh)); % find all distances less than the search radius
                va = [];
                    if ~isempty(newind) % if any are found
                        for pp = 1:length(newind) % get value of each pixel within thresh distance
                            va(pp) = im(row(newind(pp)),col(newind(pp)));
                        end
                        if sum(va) > cutoff % make sure we're not adding things smaller than the allowed cluster size
                        test = whichclust(newind(:)) > 0; % check to see if they are already part of a cluster
                        if any(test) % if yes,
                            vals = unique(whichclust(newind(test))); % figure out which clusters are represented
                            for tt = 1:length(vals)
                                whichclust(whichclust == vals(tt)) = ww; % and combine them all with a new cluster number
                            end
                        end
                        whichclust(newind) = ww; % set the cluster values of these pixels to the current cluster number
                        end
                    end
            end
            clustInds = unique(whichclust); % find out the ID numbers of each cluster
            rem = [];
            clusts = cell(size(clustInds'));
            for ss = 1:length(clustInds) % for each cluster
                if clustInds(ss) % if the value is not zero (noise)
                    clusts{ss} = [row(whichclust == clustInds(ss)),col(whichclust == clustInds(ss))]; % then fill in pixel locations for each cluster
                    newRow = [newRow;row(whichclust == clustInds(ss))]; % and add those pixels to the list of possible pixels for the next round
                    newCol = [newCol;col(whichclust == clustInds(ss))];
                else
                    rem = ss; % if value is 0, remove
                end
            end
            clusts(rem) = []; % get rid of noise pixels
            cluster = cell(size(clusts)); % make a new cell to hold the localizations within each cluster
            pixa = zeros(size(clusts)); % declare some variables ahead of time
            hll = cell(size(cluster));
            sg = cell(size(cluster));
            aspe = zeros(size(clusts));
            conv = zeros(size(clusts));
            rem = [];
            for qq = 1:length(cluster) % for each cluster
               for rr = 1:size(clusts{qq},1) % get the localizations that are in each pixel and add those to the cluster
                   cluster{qq} = [cluster{qq};locations{clusts{qq}(rr,1),clusts{qq}(rr,2)}];
               end
              if size(cluster{qq},1) > cutoff % if the cluster is larger than the minimum cutoff value
                   if size(cluster{qq},1) > 3 % has to be more than 3 points in order to do the convex hull or delaunay triangulation
                       tri = delaunay(cluster{qq}(:,1),cluster{qq}(:,2)); % get triangulation of points
                       [~,rcc] = circumcenters(TriRep(tri,cluster{qq}));
                       tri = tri(rcc < R,:);
                       rcc = rcc(rcc < R);
                       alpha = freeBoundary(TriRep(tri,cluster{qq})); % make alpha hull
                       tmpHull = cluster{qq}(alpha(:,:),:); % store the alpha hull
                       chull = convhull(cluster{qq}); % get the convex hull of the same points
                       [C ia ic] = unique(tmpHull,'rows');
                       tmpHull = tmpHull(sort(ia),:); % organize the alpha hull points into the right order for plotting
                       tmpHull = [tmpHull;tmpHull(1,:)];
                   else
                       tmpHull = cluster{qq}; % if there are fewer than 3 points (but more than cutoff) use the points as the hulls
                       chull = cluster{qq};
                   end
                hll{qq} = tmpHull; % store the alpha hull for saving
                pixa(qq) = polyarea(hll{qq}(:,1),hll{qq}(:,2)); % get the area of the alpha hull
                for ii = 1:size(tmpHull,1)-1 % get the perimeter of the alpha hull
                    slen(ii) = sqrt((tmpHull(ii,1)-tmpHull(zz+1,1))^2 + (tmpHull(ii,2)-tmpHull(ii+1,2))^2);
                end
                perm = sum(slen); % get total perimeter
                slen = []; % clear it to avoid errors later
                sg{qq} = 2*sqrt(pixa(qq)/pi); % get equivalent diameter and store
                conv(qq) = pixa(qq)/polyarea(cluster{qq}(chull,1),cluster{qq}(chull,2)); % divide alpha hull area by convex hull area to get convexity
                aspe(qq) = (4*pi*pixa(qq))/perm^2; % get compactness
             else
                rem = [rem,qq]; % get rid of clusters smaller than allowed cluster size
             end
             end
            rem = unique(rem); % remove cells representing excluded too-small clusters
            pixa(rem(:)) = [];
            aspe(rem(:)) = [];
            cluster(rem(:)) = [];
            hll(rem(:)) = [];
            sg(rem(:)) = [];
            sg = cell2mat(sg);
            conv(rem(:)) = [];
            % find nearest edge-to-edge neighbor cluster
            nD = [];
            for bb = 1:length(hll)
                 for dd = 1:length(hll)
                    for ee = 1:size(hll{bb},1)-1
                        tmpdst = Clustering.dist_Pehlke(hll{bb}(ee,:),hll{dd}(:,:));
                        tdst(ee) = min(tmpdst);
                    end
                    if any(tdst)
                        td(dd) = min(tdst(tdst>0));
                    else
                        td(dd) = 0;
                    end
                    tdst = [];
                end
                nD{bb} = min(td(td>0));
                td = [];
            end
            % store variables from this round in the proper place for
            % saving
            convty{aa}{zz} = conv;
             nD = cell2mat(nD);
            pixarea{aa}{zz} = pixa;
            compact{aa}{zz} = aspe;
            clusters{aa}{zz} = cluster;
            hull{aa}{zz} = hll;
            sig{aa}{zz} = sg;
            nearDist{aa}{zz} = nD;
            clear pixa clusts aspe cluster hll sg nD
            end
             % set min and max for plotting variables
             maxDia = max([max([sig{aa}{:}]),tmpMax]);
             tmpMax = maxDia;
             minDia = min([min([sig{aa}{:}]),tmpMin]);
             tmpMin = minDia;
             maxCon = max([max([convty{aa}{:}]),tmpCon]);
             tmpCon = maxCon;
             minCon = min([min([convty{aa}{:}]),conMin]);
             conMin = minCon;
             maxDist = max([max([nearDist{aa}{:}]),tmpDist]);
             tmpDist = maxDist;
             minDist = min([min([nearDist{aa}{:}]),distMin]);
             distMin = minDist;
             % save results to regs structure
             regs(indx(aa)).G = Ggroup{aa};
             regs(indx(aa)).compact = compact{aa};
             regs(indx(aa)).eqDiameter = sig{aa};
             regs(indx(aa)).convex = convty{aa};
             regs(indx(aa)).maxG = Gmat{aa};
             regs(indx(aa)).neardist = nearDist{aa};
             regs(indx(aa)).estMnDiam = [mean(sig{aa}{1}),mean(sig{aa}{2})];
             regs(indx(aa)).clusts = clusters{aa};
             regs(indx(aa)).ahull = hull{aa};
             regs(indx(aa)).Clarea = pixarea{aa};
            catch err
                disp('No clusters found in this region')
                convty{aa} = [];
                pixarea{aa} = [];
                compact{aa} = [];
                clusters{aa} = [];
                hull{aa} = [];
                sig{aa} = [];
                nearDist{aa} = [];
                regs.clusts = cell(1, 2);
                try
                    continue
                catch err
%                   error('Unable to find clusters with current parameters. Try lowering the minimum cluster size.')
                    fprintf('Unable to find clusters with current parameters. Try lowering the minimum cluster size.')
                end

            end
        end
    if plotting
    try
        % plot all results
%       subplot(handles.clustAx)
        figure;
        plotnum = 1;
        for dd = 1:length(indx)
            % show plot of clusters
            % first plot all points in gray and set limits, set title, etc.
            hh = subplot(num,2,plotnum);
            plot(regs(indx(dd)).points(:,1),regs(indx(dd)).points(:,2),'Color',[.5 .5 .5],'LineStyle','none','Marker','.');
            set(hh,'YDir','normal');
            set(hh,'XLim',[regs(indx(dd)).ROI(1),regs(indx(dd)).ROI(1)+regs(indx(dd)).ROI(3)])
            set(hh,'YLim',[regs(indx(dd)).ROI(2),regs(indx(dd)).ROI(2)+regs(indx(dd)).ROI(4)])
            title([regs(indx(dd)).name,', data size: ',num2str(length(regs(indx(dd)).points))])
            % then plot inner clusters with random colored points
            hold on
            if ~isempty(clusters{dd})
            for qq = 1:length(clusters{dd}{1})
                if ~isempty(clusters{dd}{1}{qq})
                    plot(clusters{dd}{1}{qq}(:,2),clusters{dd}{1}{qq}(:,1),'Color',rand([1 3]),'LineStyle','none','Marker','.');
                end
            end
            hold off
            hold on % and plot inner cluster outlines in magenta
            for mm = 1:length(hull{dd}{1})
                cvh = hull{dd}{1}{mm};
                if ~isempty(cvh)
                plot(cvh(:,2),cvh(:,1),'m')
                end
                hold on
            end
            hold off
            % Then plot outer cluster outlines in black
            hold on
            for nn = 1:length(hull{dd}{2})
                plot(hull{dd}{2}{nn}(:,2),hull{dd}{2}{nn}(:,1),'color','k','linewidth',1.5)
            end
            hold off
            plotnum = plotnum + 1;
            % histogram of cluster equivalent diameters, both inner and
            % outer
            jj(dd) = subplot(num,2,plotnum);
            bins = 0:10:maxDia;
            for ll = 1:length(sig{dd})
                barvals(ll,:) = hist(jj(dd),sig{dd}{ll},bins);
            end
            bar(bins,barvals');
            axis tight
            title([regs(indx(dd)).name,', data size: ',num2str(length(regs(indx(dd)).points))])
            xlabel({'Equivalent Diameter (nm)';['Inner mean: ',num2str(round(mean(sig{dd}{1})))];['Outer mean: ',num2str(round(mean(sig{dd}{2})))]}) % label with both means
            legend(leg,'location','best')
            plotnum = plotnum + 1;
            end
        end
        barvals = [];
        for mm = 1:length(jj)
            tmp = get(jj(mm),'YLim'); % get limits for each plot, then find the max
            maxy(mm) = tmp(2);
        end
        maxy = max(maxy);
        set([jj(:)],'Ylim',[0 maxy]) % set the limits for all plots to the max so that scale can be appreciated
        % plot convexity in the same way
        plotnum = 1;
        conBins = 0:.1:1;
        if minDist ~= maxDist
            distBins = linspace(minDist,maxDist,10);
        else
            distBins = 0:10:maxDist;
        end
        asbins = 0:.1:1;
%       subplot(handles.getisAx)
        figure;
        for tt = 1:length(indx)
            if ~isempty(clusters{tt})
            gg(tt) = subplot(num,3,plotnum);
            for ll = 1:length(convty{tt})
                barvals(ll,:) = hist(convty{tt}{ll},conBins);
            end
            bar(conBins,barvals')
            xlabel({'Convexity';' (A.U.) '})
            ylabel('Count')
            axis tight
            title([regs(indx(tt)).name,', data size: ',num2str(length(regs(indx(tt)).points))])
            legend(leg,'location','best')
            plotnum = plotnum+1;
            % then plot compactness
            ii(tt) = subplot(num,3,plotnum);
            for ll = 1:length(compact{tt})
                asvals(ll,:) = hist(compact{tt}{ll},asbins);
            end
            b = bar(asbins,asvals','BarWidth',1);
            set(b,'BarWidth',1)
            xlabel('Compactness')
            ylabel('Count')
            axis tight
            title([regs(indx(tt)).name,', data size: ',num2str(length(regs(indx(tt)).points))])
            legend(leg,'location','best')
            plotnum = plotnum+1;
            % then plot nearest neighbor distances
            kk(tt) = subplot(num,3,plotnum);
            for ll = 1:length(nearDist{tt})
                distvals(ll,:) = hist(nearDist{tt}{ll},distBins);
            end
            bar(distBins,distvals')
            xlabel('Distance to Nearest Cluster (nm)')
            ylabel('Count')
            axis tight
            title([regs(indx(tt)).name,', data size: ',num2str(length(regs(indx(tt)).points))])
            legend(leg,'location','best')
            plotnum = plotnum+1;
            end
        end
        if exist('gg','var')
        for mm = 1:length(gg) % get maximum limits for all plots
            tmp = get(gg(mm),'YLim');
            maxy(mm) = tmp(2);
            tmp2 = get(ii(mm),'YLim');
            maxy2(mm) = tmp2(2);
            tmp3 = get(ii(mm),'XLim');
            maxx(mm) = tmp3(2);
            tmp4 = get(kk(mm),'YLim');
            maxy3(mm) = tmp4(2);
        end
        end
        maxy = max(maxy); % and then set appropriately
        set([gg(:)],'Ylim',[0 maxy])
        maxy2 = max(maxy2);
        set([ii(:)],'Ylim',[0 maxy2])
        maxx = max(maxx);
        set([ii(:)],'Xlim',[0 maxx])
        maxy3 = max(maxy3);
        set([kk(:)],'YLim',[0 maxy3])

    catch err % an error will occur if you try to plot again when the axes has already been used. throw an error message but don't abort program
%       if strcmp(err.identifier,'MATLAB:subplot:InvalidAxesHandle')
%           error(err.identifier,'You already did that.')
%       else
%           error('Unable to find clusters with current parameters. Try lowering the minimum cluster size.')
            fprintf('Unable to find clusters with current parameters. Try lowering the minimum cluster size.')
%           return
%       end
    end
    end
        end
%       set(handles.pickButt,'UserData',regs)  % store regs variable
%       set(handles.workTxt,'Visible','off')   % turn off the "working" message because it is done now
    end

% distance function - get distance from point in question to all other
% points
%   function distlist = dist_Pehlke(point,mat)
%       x1 = point(1);
%       y1 = point(2);
%       for rr = 1:size(mat,1)
%           x2 = mat(rr,1);
%           y2 = mat(rr,2);
%           distlist(rr) = sqrt((x1-x2)^2+(y1-y2)^2);
%       end
%   end

% chi-square test function
%   function [H p] = doChi2(grp1,grp2)
%      chi2 = sum(((grp1 - grp2).^2)./(grp1 + grp2));
%      p = 1 - chi2cdf(chi2,1);
%      H = p <= .05;
%   end

    function [Ggroup Gmat] = calcGetis(mBW,cutval,gx,gy)
% function for performing getis analysis
        v = mBW(:); % get locations of all pixels in image
        v = v(v>0); % only consider non-zero pixels
        v = double(v);
        for qq = 1:length(cutval) % loop through cutoff distances
            w = Clustering.weight_matrix([gy gx],2,cutval(qq)); % make weight matrix for each cutoff
            [G(:,qq),Gstar(:,qq)] = Clustering.getis_statistic(w,v); % do getis statistic
            tmp = zeros(size(mBW)); % make blank image size of input image
            for uu = 1:length(gx)
                tmp(gx(uu),gy(uu)) = G(uu,qq);  % get the G value for each pixel at this cutoff distance
                Gmat(:,:,qq) = tmp; % store in 3D matrix of G values
            end
        end
       v = [];
       Ggroup = G;
    end

    function [K dr] = sR_Ripley(cutoff,regs,rate)
% ripley's function
    % Written by Michael Wester and Stanly Steinberg in 2008.
    % step size
    dr = cutoff/rate;
        % for each region
        for aa = 1:length(regs)
            area = regs(aa).area;
            Xx = regs(aa).points;
            R = zeros(1,rate);
            for ii = 1:length(Xx)
                for jj = 1:ii-1
                    p = ceil(sqrt((Xx(jj,1)-Xx(ii,1))^2 + (Xx(jj,2)-Xx(ii,2))^2)/dr);
                    if p > 0 && p <= rate
                        R(p) = R(p) + 1;
                    end
                end
            end
            R = [0,R];
        R = 2*R;

        % convert pdf to cdf
        for bb = 2:length(R)
            R(bb) = R(bb) + R(bb-1);
        end

        % average over number of particles
        R = R/length(Xx);
        % compute the intensity and normalize R
        lambda = length(Xx)/area;
        R = R/lambda;
        K{aa} = R;
        end
    end

    function H = sR_hopkinstat(P,A,B,mm)
% Hopkin's function
% Written by Michael Wester and Stanly Steinberg in 2008.
        % Compute a Hopkin's statistic for the particles P.
        % The particle region is [0,A] x [0,B] .
        % m is the number of tests points.
        % n is the number of particles.
        %n = length(P);
        n = numel(P) / 2;
        m = min(mm, n);
        % Create the indices of m test particles.
        index = [0,m+1];
        while ( (length(index) < m) || (index(1) < 1) || (index(end) > n) )
           index = unique(round(1/2+n*rand(m,1)));
        end

        tmp1 = P(:,1);
        tmp1 = (tmp1 - min(tmp1))/(max(tmp1)-min(tmp1));
        P(:,1) = tmp1*A;

        tmp2 = P(:,2);
        tmp2 = (tmp2 - min(tmp2))/(max(tmp2)-min(tmp2));
        P(:,2) = tmp2*B;


        T = P(index,:);
        % Create m test points
        S = rand(m,2)*diag([A,B]);
        % Compute the minimum distance.
        U = ones(1,m)*sqrt(A^2+B^2);
        for k = 1:m
            for i = 1:n
                    dist = sqrt((S(k,1)-P(i,1))^2 +(S(k,2)-P(i,2))^2);
                    if dist > 0
                        U(k) = min(U(k), dist);
                    end
            end
        end
        % Compute the minimum distance.
        W = ones(1,m)*sqrt(A^2+B^2);
        for k = 1:m
            for i = 1:n
                    dist = sqrt((T(k,1)-P(i,1))^2 +(T(k,2)-P(i,2))^2);
                    if dist > 0
                        W(k) = min(W(k), dist);
                    end
            end
        end
        H = (sum(U.^2)/(sum(U.^2)+sum(W.^2)));
    end

    function [G, Gstar] = getis_statistic(w, v)
% function for calculating Getis statistic
% written by Michael Wester 2013
% J. K. Ord and Arthur Getis, ``Local Spatial Autocorrelation Statistics:
% Distributional Issues and an Application'', _Geographical Analysis_, Volume
% 27, Number 4, October 1995, 286--306.
%
% For n data points,
%    w is a n x n (sparse) weight matrix
%    x is a n x d matrix of coordinates, one set per row
%    v is a n x 1 vector of values

           n = length(v);

           Gstar_n = w * v;
           Gstar_d = sum(v);
           if abs(Gstar_d) <= 5*eps
              fprintf('WARNING (getis): Gstar denominator nearly zero!   %f\n', ...
                      Gstar_d);
           end
           Gstar = Gstar_n ./ Gstar_d;

           G = zeros(n, 1);
           for i = 1 : n
              G(i) = (Gstar_n(i) - w(i, i)*v(i)) ./ (Gstar_d - v(i));
           end
    end

    function w = weight_matrix(x, method, d)
% function for creating Getis weight matrix
% written by Michael Wester 2013

       n = length(x);

       % Compute the distances between the points.
       D = squareform(pdist(x));

       switch method
          case 1,   % pure distances
             w = D;
          case 2,   % 1 if distance <= d, otherwise 0
             w = D <= d;
          case 3,   % distance decay
             w = min(1 ./D, 1);
          case 4,   % exponential distance decay
             w = exp(-d/100 * D);
          otherwise
             fprintf('ERROR (weight_matrix): Unknown method!   %d\n', method);
       end
    end

% =============================================================================

function [C, ptsI] = hierarchal(XY, E, minPts)
% Form clusters such that any point in a cluster is within E of some other
% point in the same cluster.
%
% Inputs:
%    XY       matrix of coordinates [N x 2]
%    E        cutoff distance
%    minPts   minimum number of points required for a cluster
%
% Outputs:
%    C        cell array of indices (wrt XY) per cluster [N x 1]
%    ptsI     indices (wrt XY) of isolated points

   if size(XY, 1) == 1
      C{1} = [1];
      ptsI = [];
      return;
   end

   Z = linkage(XY, 'single');
   %figure; dendrogram(Z);
   T = cluster(Z, 'Cutoff', E, 'Criterion', 'distance', 'Depth', 2);
   nC = max(T);

   % Remove clusters of size < minPts, taking the points to be isolated.
   C = [];
   j = 0;
   ptsI = [];
   for i = 1 : nC
      c = find(T == i);
      n = length(c);
      if n >= minPts
         j = j + 1;
         C{j} = c';
      else
         ptsI = [ptsI, c'];
      end
   end
   nC = j;

   ptsI = sort(ptsI);

%  % Check if the isolated points really are isolated or if they can be added
%  % to an existing cluster.  Do a single check here, although this might be
%  % iterated in the general case (or the algorithm rewritten entirely).
%  E2 = E^2;
%  not_isolated = [];
%  for i = 1 : length(ptsI)
%     I = ptsI(i);
%     p = XY(I, :);
%     MIN_d2 = 1.0e10;  
%     for j = 1 : nC
%        c = XY(C{j}, :);
%        min_d2 = min((c(:, 1) - p(1)).^2 + (c(:, 2) - p(2)).^2);
%        if min_d2 < MIN_d2
%           MIN_d2 = min_d2;
%           indx = j;
%        end
%     end
%     % Point is not isolated after all!  Add it into the nearest cluster.
%     if MIN_d2 < E2
%        C{indx} = [C{indx}, I];
%        not_isolated = [not_isolated, I];
%     end
%  end
%  % Removed non-isolated points (as determined above) from the list of
%  % isolated points.
%  ptsI = setdiff(ptsI, not_isolated);

end

% =============================================================================

function plot_voronoi(X, Y, v, c, rho, str, dense, ptIDs)
% Plot the Voronoi diagram corresponding to (X, Y), coloring the cells
% according to the density rho.  Also,
%    v       vertices [N x 2]
%    c       vertex indices cooresponding to each cell [Nc x 1]
%    rho     density of each cell [Nc x 1]
%    dense   indices of those cells >= prescribed density criterion
%    ptIDs   write out point IDs if true

   n_XY = length(X);
   n_v  = size(v, 1);
   n_c  = length(c);

   delta = 0.001 * (max(X) - min(X));
   figure();
   hold on
   voronoi(X, Y);
   limits = axis;
   colormap(jet);
   caxis([min(rho), max(rho)]);
   for i = 1 : n_c
      c_i = c{i};
      if all(c_i ~= 1)
         fill(v(c_i, 1), v(c_i, 2), rho(i));
         %fill(v(c_i, 1), v(c_i, 2), floor((i - 1)/(n_c - 1) * 255 + 1));
      end
   end
   plot(X, Y, 'k.', 'MarkerSize', 10);
   plot(X(dense), Y(dense), 'r.', 'MarkerSize', 10);
   if exist('ptIDs', 'var') & ptIDs
      for i = 1 : n_XY
         text(X(i) + delta, Y(i), sprintf('%d', i), 'Color', 'k');
      end
   end
%  plot(v(:, 1), v(:, 2), 'g.', 'MarkerSize', 10);
%  for i = 2 : n_v
%     text(v(i, 1) + delta, v(i, 2), sprintf('%d', i), 'Color', 'g');
%  end
   axis(limits);
   cb = colorbar;
   cb.Label.String = 'relative density';
   title(str);
   hold off

end

% -----------------------------------------------------------------------------

function plot_voronoi3(X, Y, Z, v, c, rho, str, dense, ptIDs)
% Plot the Voronoi diagram corresponding to (X, Y, Z), coloring the cells
% according to the density rho.  Also,
%    v       vertices [N x 2]
%    c       vertex indices cooresponding to each cell [Nc x 1]
%    rho     density of each cell [Nc x 1]
%    dense   indices of those cells >= prescribed density criterion
%    ptIDs   write out point IDs if true

   n_XY = length(X);
   n_v  = size(v, 1);
   n_c  = length(c);

   delta = 0.001 * (max(X) - min(X));
   figure();
   hold on
   dt = delaunayTriangulation(X, Y, Z);
   tetramesh(dt);
   limits = axis;
   colormap(jet);
   caxis([min(rho), max(rho)]);
   for i = 1 : n_c
      c_i = c{i};
      if all(c_i ~= 1)
         fill3(v(c_i, 1), v(c_i, 2), v(c_i, 3), rho(i));
         %fill(v(c_i, 1), v(c_i, 2), floor((i - 1)/(n_c - 1) * 255 + 1));
      end
   end
   plot3(X, Y, Z, 'k.', 'MarkerSize', 10);
   plot3(X(dense), Y(dense), Z(dense), 'r.', 'MarkerSize', 10);
   if exist('ptIDs', 'var') & ptIDs
      for i = 1 : n_XY
         text(X(i) + delta, Y(i), Z(i), sprintf('%d', i), 'Color', 'k');
      end
   end
%  plot(v(:, 1), v(:, 2), 'g.', 'MarkerSize', 10);
%  for i = 2 : n_v
%     text(v(i, 1) + delta, v(i, 2), sprintf('%d', i), 'Color', 'g');
%  end
   axis(limits);
   cb = colorbar;
   cb.Label.String = 'relative density';
   title(str);
   hold off

end

% -----------------------------------------------------------------------------

function [nC, C] = cluster_voronoi(i_rho, self_nbrs, epsilon, minPts, XY)
% Taking the density indices i_rho that identify points to be clustered,
% generate the clusters C (their number given by nC).
%
% Inputs:
%    i_rho       indices of all the cells exceeding the density criterion
%    self_nbrs   for each cell, the indices for itself and all its neighbors
%    epsilon     epsilon or cutoff distance
%    minPts      minimum number of points allowed in a cluster
%    XY          point coordinates [N x 2]
% Outputs:
%    nC          number of clusters found
%    C           cell array of XY indices forming each cluster [nC x 1]

   nC = 0;
   C = {};
   while ~isempty(i_rho)
      % There is at least one more point left to be clustered, so create a new
      % cluster and stuff the point in it while deleting it off the list of
      % points remaining to be clustered.
      C_nC = i_rho(1); 
      i_rho(1) = [];
      % Find the point's neighbors and see if any of them are on the list of
      % points remaining to be clustered.
      lo = 1;
      hi = 1;
      nbrs = self_nbrs{C_nC(lo:hi)};
      l = intersect(nbrs, i_rho);
      while ~isempty(l)
         % If there is an overlap between the neighbors of the points just
         % added to the current cluster and those remaining to be clustered,
         % stuff the overlap in the current cluster, delete them off the list
         % of points remaining to be clustered, and then compute the overlap
         % between the neighbors of these new cluster points and those
         % remaining.
         C_nC = [C_nC, l];
         i_rho = setdiff(i_rho, l);
         lo = hi + 1;
         hi = hi + length(l);
         % Find the unique neighbors of the new points just added to the
         % cluster (and themselves).
         nbrs = unique([ self_nbrs{C_nC(lo:hi)} ]);
         l = intersect(nbrs, i_rho);
      end
      % Eliminate clusters smaller than minPts.
      if length(C_nC) >= minPts
         nC = nC + 1;  
         C{nC} = sort(C_nC);
      end
   end

   % If epsilon > 0, apply further restrictions on the clusters found above.
   % For each Voronoi cluster, apply a secondary clustering algorithm that
   % separates points based on epsilon.  This will, in general, separate the
   % Voronoi clusters into smaller (or same size) clusters and isolated points.
   % Collect together the newly separated clusters and return these.
   if epsilon > 0
      c = Clustering();
      Algorithm = 'Hierarchal';
      nB = 0;
      B  = {};
      for i = 1 : nC
         xy = XY(C{i}, :);
         [nCC, CC, ~, ~] = c.cluster(Algorithm, xy, epsilon, minPts);
         for j = 1 : nCC
            nB = nB + 1;
            % CC{j} are the point indices of a new cluster, so map these back
            % into the point indices of the original cluster C{i}.
            B{nB} = C{i}(CC{j});
         end
      end

      nC = nB;
      C  = B;
   end

end

% =============================================================================

function results = clusterStats(XY, C, centers, shrinkFactor)
% Compute statistics on computed clusters.
%
% Inputs:
%    XY             coordinates of all the points that were processed
%    C              cell array of the indices of the points in each cluster
%    centers        array of the coordinates of the center of each cluster
%
% Outputs (contained in results):
%    nC                number of clusters
%    n_points          total number of points
%    n_clustered       number of points in clusters
%    n_isolated        number of points not in clusters
%    n_pts             number of points per cluster
%    numclust(1,2,3)   number of singlet, double, multiple clusters, where
%                      singlet clusters include isolated points (see
%                      SRcluster.m for an equivalent definition)
%    singlet_faction   numclust(1) / sum(numclust)
%    sigma_actual      actual (computed) sigma of each cluster, that is, the
%                      standard deviation of the intracluster distances
%    indices_hull      cell array of boundary hull indices relative to XY per
%                      cluster
%    areas             area of each cluster
%    equiv_radii       equivalent radius of each cluster
%    n_pts_per_area    number of points per area for clusters containing 3 or
%                      more points
%    perimeters        perimeter of each cluster
%    compactness       4 pi area / perimeter^2 of each cluster
%    min_c2c_dists     minimum center-to-center distances for each cluster
%                      with respect to all the others
%    min_e2e_dists     minimum edge-to-edge distances for each cluster convex
%                      hull with respect to all the others
%    min_c2c_dist      min(min_c2c_dists)
%    min_e2e_dist      min(min_e2e_dists)

   if exist('shrinkFactor', 'var')
      ShrinkFactor = shrinkFactor;
   else
      ShrinkFactor = 0.5;
   end

   dim = size(XY, 2);

   nC = length(C);
   results.nC = nC;

   xy_hull       = cell(1, nC);
   min_e2e_dists = zeros(1, nC);
   results.n_points     = size(XY, 1);
   results.n_clustered  = 0;
   results.n_isolated   = 0;
   results.n_pts        = zeros(1, nC);
   results.numclust     = zeros(1, 3);
   results.sigma_actual = zeros(1, nC);
   results.indices_hull = cell(1, nC);
   results.areas        = zeros(1, nC);
   results.equiv_radii  = zeros(1, nC);
   results.n_pts_per_area = [];
   results.perimeters   = [];
   results.compactness  = [];
   for i = 1 : nC
      xy = double(XY(C{i}, :));
      n_pts = size(xy, 1);
      results.n_pts(i) = n_pts;
      results.n_clustered = results.n_clustered + n_pts;
      if n_pts == 1
         results.numclust(1) = results.numclust(1) + 1;
         results.sigma_actual(i) = 0;
         results.indices_hull{i} = 1;
         results.areas(i)        = 0;
         results.equiv_radii(i)  = 0;
         results.compactness(i)  = 1;
      elseif n_pts == 2
         results.numclust(2) = results.numclust(2) + 1;
         results.sigma_actual(i) = std(pdist(xy));
         results.indices_hull{i} = [1, 2];
         results.areas(i)        = 0;
         results.equiv_radii(i)  = pdist(xy) / 2;
         results.compactness(i)  = 0;
      else % n_pts >= 3
         results.numclust(3) = results.numclust(3) + 1;
         results.sigma_actual(i) = std(pdist(xy));
         try
            %[k, A] = convhull(xy(:, 1), xy(:, 2));
            if dim == 2
               [k, A] = boundary(xy(:, 1), xy(:, 2), ShrinkFactor);
            else
               [k, A] = boundary(xy(:, 1), xy(:, 2), xy(:, 3), ShrinkFactor);
            end
            if isempty(k)
               k = 1;
            end
         catch
            fprintf('boundary collinear (n_points = %d)\n', n_pts);
            %xy
            k = 1;
            A = 0;
         end
         results.indices_hull{i} = k;
         results.areas(i)        = A;
         if dim == 2
            % A = pi r^2
            results.equiv_radii(i)  = sqrt(A / pi);
         else
            % V (A) = 4/3 pi r^3
            results.equiv_radii(i)  = (3/4*A / pi)^(1/3);
         end
         if A > 0
            results.n_pts_per_area  = [results.n_pts_per_area, n_pts / A];
            if dim == 2
               perim = 0;
               for j = 1 : length(k) - 1
                  perim = perim + pdist([xy(k(j), :); xy(k(j + 1), :)]);
               end
               results.perimeters      = [results.perimeters, perim];
               results.compactness     = [results.compactness, ...
                                          4*pi*A / perim^2];
            end
         end
      end
      xy_hull{i} = xy(results.indices_hull{i}, :);
   end
   results.n_isolated = results.n_points - results.n_clustered;
   results.numclust(1) = results.numclust(1) + results.n_isolated;
   results.singlet_fraction = results.numclust(1) / sum(results.numclust);

   % The below can be an expensive operation timewise.
   if nC > 1 
      %results.min_c2c_dists = ...
      %   min(squareform(pdist(centers')) + 1.0e+10 * eye(nC));
      results.min_c2c_dists = Clustering.nn_distances(centers');
      results.min_c2c_dist = min(results.min_c2c_dists);

      results.min_e2e_dists = zeros(1, nC);
      min_e2e_dist = 1.0e+10;
      for i = 1 : nC
         min_e2e_dists(i) = 1.0e+10;
         for j = 1 : i - 1
            e2e = Clustering.edge2edge(xy_hull{i}, xy_hull{j});
            min_e2e_dist = min(min_e2e_dist, e2e);
            min_e2e_dists(i) = min(min_e2e_dists(i), e2e);
         end
         for j = i + 1 : nC
            e2e = Clustering.edge2edge(xy_hull{i}, xy_hull{j});
            min_e2e_dist = min(min_e2e_dist, e2e);
            min_e2e_dists(i) = min(min_e2e_dists(i), e2e);
         end
         results.min_e2e_dists(i) = min_e2e_dists(i);
      end
      results.min_e2e_dist = min_e2e_dist;
   else
      results.min_c2c_dists = [];
      results.min_c2c_dist  = [];
      results.min_e2e_dists = [];
      results.min_e2e_dist  = [];
   end

end

% -----------------------------------------------------------------------------

function clusterFig = plotClusters(XY, C, centers, ptsI, txt, shrinkFactor, ...
                                   options)
% Plot and label the 2D clusters.
%
% Inputs:
%    XY             point coordinates [N x 2]
%    C              cell array of XY indices forming each cluster [nC x 1]
%    centers        coordinates of the center of each cluster [nC x n_dim]
%    ptsI           indices of points not found in any cluster
%    txt            descriptive text to add to the plot's title
%    shrinkFactor   boundary shrink factor (0 = convex hull, 1 = as concave as
%                   possible, 0.5 = MATLAB default)
%    options        'L': label each cluster
%                   'O': outline each cluster
%                   'P': print the size of each cluster
%                   '1': use only one color for all clusters
% Output:
%    clusterFig   figure handle

   if exist('shrinkFactor', 'var')
      ShrinkFactor = shrinkFactor;
   else
      ShrinkFactor = 0.5;
   end

   if exist('options', 'var')
      labeling = strfind(options, 'L');
      outlines = strfind(options, 'O');
      printing = strfind(options, 'P');
      onecolor = strfind(options, '1');
   else
      labeling = false;
      outlines = true;
      printing = false;
      onecolor = false;
   end

   colors = 'rgbcmy';
   n_colors = length(colors);

   nC = length(C);
   n_isolated = length(ptsI);
   n_points = size(XY, 1);
   n_clustered = 0;

   clusterFig = figure('Visible', 'off');
   hold on
   fprintf('\n');
   for i = 1 : nC
      j = Clustering.nMODm(i, n_colors);
      if onecolor
         color = 'm';
      else
         color = colors(j);
      end
      n_pts = length(C{i});
      n_clustered = n_clustered + n_pts;
      if printing
         fprintf('cluster %d has %2d points (%s)\n', i, n_pts, color);
      end
      plot(XY(C{i}, 1), XY(C{i}, 2), [color, '.'], 'MarkerSize', 12);
      if labeling
         text(centers(1, i), centers(2, i), sprintf('%d', i));
      end
      if outlines
         if length(C{i}) <= 2
            plot(XY(C{i}, 1), XY(C{i}, 2), [color, '-'], 'LineWidth', 2);
         else %if length(C{i}) >= 3
            xy = double(XY(C{i}, :));
            %k = convhull(XY(C{i}, 1), XY(C{i}, 2));
            k = boundary(xy(:, 1), xy(:, 2), ShrinkFactor);
            k = C{i}(k);
            plot(XY(k, 1), XY(k, 2), [color, '-'], 'LineWidth', 2);
         end
      end
   end
   if n_isolated > 0
      plot(XY(ptsI, 1), XY(ptsI, 2), 'k.', 'MarkerSize', 12);
   end
   title(sprintf('%s (clusters = %d, clustered %% = %.3f)', ...
                regexprep(txt, '_', '\\_'), nC, n_clustered/n_points * 100));
   hold off

end

% -----------------------------------------------------------------------------

function clusterFig = plotClusters3(XY, C, centers, ptsI, txt, shrinkFactor, ...
                                    options)
% Plot and label the 3D clusters.
%
% Inputs:
%    XY             point coordinates [N x 2]
%    C              cell array of XY indices forming each cluster [nC x 1]
%    centers        coordinates of the center of each cluster [nC x n_dim]
%    ptsI           indices of points not found in any cluster
%    txt            descriptive text to add to the plot's title
%    shrinkFactor   boundary shrink factor (0 = convex hull, 1 = as concave as
%                   possible, 0.5 = MATLAB default)
%    options        'L': label each cluster
%                   'O': outline each cluster
%                   'P': print the size of each cluster
%                   '1': use only one color for all clusters
% Output:
%    clusterFig   figure handle

   if exist('shrinkFactor', 'var')
      ShrinkFactor = shrinkFactor;
   else
      ShrinkFactor = 0.5;
   end

   if exist('options', 'var')
      labeling = strfind(options, 'L');
      outlines = strfind(options, 'O');
      printing = strfind(options, 'P');
      onecolor = strfind(options, '1');
   else
      labeling = false;
      outlines = true;
      printing = false;
      onecolor = false;
   end

   colors = 'rgbcmy';
   n_colors = length(colors);

   nC = length(C);
   n_isolated = length(ptsI);
   n_points = size(XY, 1);
   n_clustered = 0;

   clusterFig = figure('Visible', 'off');
   hold on
   fprintf('\n');
   for i = 1 : nC
      j = Clustering.nMODm(i, n_colors);
      if onecolor
         color = 'm';
      else
         color = colors(j);
      end
      n_pts = length(C{i});
      n_clustered = n_clustered + n_pts;
      if printing
         fprintf('cluster %d has %2d points (%s)\n', i, n_pts, color);
      end
      plot3(XY(C{i}, 1), XY(C{i}, 2), XY(C{i}, 3), [color, '.'], ...
            'MarkerSize', 12);
      if labeling
         text(centers(1, i), centers(2, i), centers(3, i), sprintf('%d', i));
      end
      if outlines
         if length(C{i}) <= 3
            plot3(XY(C{i}, 1), XY(C{i}, 2), XY(C{i}, 3), [color, '-'], ...
                  'LineWidth', 2);
         else %if length(C{i}) >= 4
            xy = double(XY(C{i}, :));
            %k = convhull(XY(C{i}, 1), XY(C{i}, 2), XY(C{i}, 3));
            k = boundary(xy(:, 1), xy(:, 2), xy(:, 3), ShrinkFactor);
            k = C{i}(k);
            plot3(XY(k, 1), XY(k, 2), XY(k, 3), [color, '-'], 'LineWidth', 2);
            %trisurf(k, XY(:, 1), XY(:, 2), XY(:, 3), 'FaceColor', color);
         end
      end
   end
   if n_isolated > 0
      plot3(XY(ptsI, 1), XY(ptsI, 2), XY(ptsI, 3), 'k.', 'MarkerSize', 12);
   end
   title(sprintf('%s (clusters = %d, clustered %% = %.3f)', ...
                regexprep(txt, '_', '\\_'), nC, n_clustered/n_points * 100));
   hold off

end

% -----------------------------------------------------------------------------

function clusterFig = ...
   plotLabelClusters(xy, c, Centers, PtsI, txt, shrinkFactor)
% Plot the 2D clusters for sets of clusters representing different labels,
% using a single color per label.
%
% Inputs: (nL = number of labels)
%                  (nL x 1) cell arrays of
%    xy                point coordinates [N x 2]
%    c                 cell array of XY indices forming each cluster [nC x 1]
%    Centers           coordinates of the center of each cluster [nC x n_dim]
%    PtsI              indices of points not found in any cluster
%    txt            descriptive text to add to the plot's title
%    shrinkFactor   boundary shrink factor (0 = convex hull, 1 = as concave as
%                   possible, 0.5 = MATLAB default)
% Output:
%    clusterFig   figure handle

   if exist('shrinkFactor', 'var')
      ShrinkFactor = shrinkFactor;
   else
      ShrinkFactor = 0.5;
   end

   outlines = true;
   colors = 'gmkcryb';

   clusterFig = figure('Visible', 'off');
   hold on
   for j = 1 : numel(xy)
      XY = xy{j};
      C = c{j};
      centers = Centers{j};
      ptsI = PtsI{j};

      nC = length(C);
      n_isolated = length(ptsI);
      n_points = size(XY, 1);
      n_clustered = 0;
      color = colors(j);

      for i = 1 : nC
         n_pts = length(C{i});
         n_clustered = n_clustered + n_pts;
         plot(XY(C{i}, 1), XY(C{i}, 2), [color, '.'], 'MarkerSize', 12);
         if outlines
            if length(C{i}) <= 2
               plot(XY(C{i}, 1), XY(C{i}, 2), [color, '-'], 'LineWidth', 2);
            else
               k = boundary(XY(C{i}, 1), XY(C{i}, 2), ShrinkFactor);
               k = C{i}(k);
               plot(XY(k, 1), XY(k, 2), [color, '-'], 'LineWidth', 2);
            end
         end
      end
      if n_isolated > 0
         plot(XY(ptsI, 1), XY(ptsI, 2), 'k.', 'MarkerSize', 12);
      end
   end
   title(sprintf('%s', regexprep(txt, '_', '\\_')));
   hold off

end

% -----------------------------------------------------------------------------

function r = nMODm(n, m)
% Modulus such that r is in [1, m] rather than [0, m - 1].

   r = mod(n, m);
   if r <= 0
      r = r + m;
   end

end

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
   min_d = [];
   if size(xy, 1) > 1
      [~, min_d] = knnsearch(xy, xy, 'K', 2);
      min_d = min_d(:, 2)';
   end

end

% -----------------------------------------------------------------------------

function min_d = edge2edge(hull1, hull2)
% Minimum edge-to-edge distance between hull 1 and hull 2.
%
% Inputs:
%    hull1, hull2   (x, y {, z}) coordinates of the two hulls defining the
%                   boundaries of a cluster of points [N x dim]
% Outputs:
%    min_d          minimum edge-to-edge distance between hull 1 and hull 2

   dim = size(hull1, 2);
   min_d2 = 1.0e+10;
   n1 = size(hull1, 1);
   %n2 = size(hull2, 1);
   for i = 1 : n1
      %d2 = (hull1(i, 1) - hull2(:, 1)).^2 + (hull1(i, 2) - hull2(:, 2)).^2;
      d2 = 0;
      for j = 1 : dim
         d2 = d2 + (hull1(i, j) - hull2(:, j)).^2;
      end
      min_d2 = min(min_d2, min(d2));
   end
   min_d = sqrt(min_d2);

end

% =============================================================================

function H = hopkinstat(P,A,B,mm)
% Compute a Hopkin's statistic for the 2D particles P.
% Written by Michael Wester and Stanly Steinberg in 2008.
% The particle region is [0,A] x [0,B] .
% m is the number of tests points.
% n is the number of particles.
%n = length(P);
n = numel(P) / 2;
m = min(mm, n);
% Create the indices of m test particles.
index = [0,m+1];
while ( (length(index) < m) | (index(1) < 1) | (index(end) > n) )
   index = unique(round(1/2+n*rand(m,1)));
end
T = P(index,:);
% Create m test points
S = rand(m,2)*diag([A,B]);
% Compute the minimum distance.
U = ones(1,m)*sqrt(A^2+B^2);
% W = ones(1,m)*sqrt(A^2+B^2);
W = U;
for k = 1:m
   for i = 1:n
      dist = sqrt((S(k,1)-P(i,1))^2 + (S(k,2)-P(i,2))^2);
      if dist > 0 
         U(k) = min(U(k), dist);
      end

      dist = sqrt((T(k,1)-P(i,1))^2 + (T(k,2)-P(i,2))^2);
      if dist > 0 
         W(k) = min(W(k), dist);
      end
   end
end
H = (sum(U.^2)/(sum(U.^2)+sum(W.^2)));

end

% -----------------------------------------------------------------------------

function H = hopkinstat3(P,A,B,C,mm)
% Compute a Hopkin's statistic for the 3D particles P.
% Written by Michael Wester and Stanly Steinberg in 2008.
% The particle region is [0,A] x [0,B] x [0,C].
% m is the number of tests points.
% n is the number of particles.
%n = length(P);
n = numel(P) / 3;
m = min(mm, n);
% Create the indices of m test particles.
index = [0,m+1];
while ( (length(index) < m) | (index(1) < 1) | (index(end) > n) )
   index = unique(round(1/2+n*rand(m,1)));
end
T = P(index,:);
% Create m test points
S = rand(m,3)*diag([A,B,C]);
% Compute the minimum distance.
U = ones(1,m)*sqrt(A^2+B^2+C^2);
% W = ones(1,m)*sqrt(A^2+B^2+C^2);
W = U;
for k = 1:m
   for i = 1:n
      dist = sqrt((S(k,1)-P(i,1))^2 + (S(k,2)-P(i,2))^2 + (S(k,3)-P(i,3))^2);
      if dist > 0 
         U(k) = min(U(k), dist);
      end

      dist = sqrt((T(k,1)-P(i,1))^2 + (T(k,2)-P(i,2))^2 + (T(k,3)-P(i,3))^2);
      if dist > 0 
         W(k) = min(W(k), dist);
      end
   end
end
H = (sum(U.^3)/(sum(U.^3)+sum(W.^3)));

end

% -----------------------------------------------------------------------------

function [X,V]=histogram(A,ab,n)
% Written by Michael Wester and Stanly Steinberg in 2008.
% Sort A into bins whose edges are given by ab(1)+i*(ab(2)-ab(1))/n, i=0:n.
L = ab(2)-ab(1);
X = linspace(L/(2*n),1-L/(2*n),n);
V=zeros(1,n);
for i=1:length(A)
	k=floor(n*(A(i)-ab(1))/L)+1;
	V(k)=V(k) + 1;
end

end

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
