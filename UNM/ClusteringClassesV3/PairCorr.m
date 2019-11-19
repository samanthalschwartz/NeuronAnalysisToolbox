classdef PairCorr < handle

% PairCorr class written by Michael Wester, Keith Lidke, Carolyn Pehlke and
%    others as noted internally (2/26/2018) <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2018 by Michael J. Wester and Keith A. Lidke
%
% Example main program:
%
%    pc = PairCorr();
%    pc.Results = 'results';
%
%    [x, y] = textread('../data/9021_5.txt',  '%*u %u %u %*u', ...
%                                             'headerlines', 1);
%    XY_5 =  [x, y];
%    [x, y] = textread('../data/9021_10.txt', '%*u %u %u %*u', ...
%                                             'headerlines', 1);
%    XY_10 = [x, y];
%    pc.ROI = [0, 0, 7400, 6000];
%    results_pcc  = pc.pair_correlation(2.7559, '9021', XY_5, XY_10)
%    results_pac1 = pc.pair_correlation(2.7559, '9021_5',  XY_5)
%    results_pac2 = pc.pair_correlation(2.7559, '9021_10', XY_10)
%
%    results_Vpcc = ...
%       pc.pair_correlation_Veatch(XY_5, XY_10, 2.7559, '9021', 'cross');
%    results_Vpac1 = ...
%       pc.pair_correlation_Veatch(XY_5,  [], 2.7559, '9021_5', 'auto');
%    results_Vpac2 = ...
%       pc.pair_correlation_Veatch(XY_10, [], 2.7559, '9021_10','auto');

% =============================================================================
properties
% =============================================================================

   % If ROI is provided, it will be used, otherwise if ROI_size is defined,
   % then the (x_min, y_min) computed from the data will be used, otherwise if
   % neither is provided, then the xy_size will be calculated using the
   % (x_min, y_min, x_max, y_max) computed from the data as
   %    xy_size = min(x_max - x_min, y_max - y_min)
   ROI      = [];   % [x_min, y_min, x_size, y_size]
   ROI_size = [];   % Set x_size = y_size = ROI_size

   Font_props = {'FontSize', 15, 'FontWeight', 'bold'};
   %Line_props = {'LineWidth', 3};
   % If Fig_ext is empty, plot to the screen and save as .fig; otherwise, do
   % not plot to the screen, but save as .fig AND .Fig_ext .
   Fig_ext = 'pdf';
   Lines   = true;  % If true, plot lines rather than points for g(r) vs. r
   Results = '.';   % Directory to store results.

   Rmax_axis = -1;   % Sets plotting limit if > 0
   % Default fit model for pair_correlation_Veatch.  This can also be supplied
   % as the optional last (6th) argument to this function.
   Veatch_fit = 'exponential_and_gaussian';
   %Veatch_fit = 'exponential_and_cosine';
   %Veatch_fit = 'exponential';
   % Note that pair_correlation and pair_correlation_ROIcombined use a 2D
   % Gaussian fit model.

   % Properties for {PC,PCC,PAC}_SR:
   Verbose = false;   % verbose output and extra saved .mat files
   % Perform auto-correlations on each ROI as well as cross-correlations.
   Auto    = false;   % only applicable to PC_SR
   % Use the original Veatch code for pairwise correlation as well as code
   % adapted from SuperCluster (originally written by Carolyn Pehlke and Keith
   % Lidke).
   Veatch  = false;
   % If true, correlate on ROIs rather than the entire image.  This is usually
   % preferred except for small datasets.
   ROIs    = true;
   Pixel2nm = 104;   % conversion from pixels to nm
   % Histogram bin size for pairwise correlation---this is the number of pixels
   % per bin over which correlation statistics are collected.
   Hist_bin_size = 104 / 20;

   % H-SET pass 2 parameters (in which clustering is performed):
   % If true, perform H-SET on observations.  H-SET (Hierarchical Single
   % Emitter Hypothesis Test), is a top-down hierarchical clustering algorithm
   % that collapses clusters of observations of blinking fluorophores into
   % single estimates of the true locations (localizations) of the
   % fluorophores.
   HSET = false;
   Algorithm = 'DBSCAN_Daszykowski';
   E = 30;       % epsilon for clustering (nm)
   MinPts = 3;   % minimum number of points for clustering
   % The shrink factor is in the interval [0, 1].
   % ShrinkFactor = 0 produces convex hulls around clustered points.
   % ShrinkFactor = 1 produces compact boundaries around clustered points,
   %                  which may be very concave.
   % The MATLAB default value for the shrink factor is 0.5
   ShrinkFactor = 0.5;

% =============================================================================
end % properties

methods
% =============================================================================
function PC_SR(obj, SRDR1, SRDR2, basefile)
% Perform pair correlation before and after H-SET collapse on the given input
% SRD.Results datasets.
%
% Cross-correlation and optionally auto-correlation are performed.
%
% Use the mouse to select ROIs (regions of interest):
%    left click  chooses the center of a fixed size (x_size x y_size) region
%    right click chooses an adjustable rectangular size region
%    key press:
%       backspace or delete   deletes the previous region
%       anything else         terminates selection
%
% Inputs:
%    SRDR1      SRD.Results dataset label 1 
%    SRDR2      SRD.Results dataset label 2
%    basefile   base filename to which various text will be appended to
%               identify the results

   % Results directory.
   if ~isdir(obj.Results)
      mkdir(obj.Results);
   end
   base = fullfile(obj.Results, basefile);

   % ==========================================================================

   if obj.HSET
      SRc = SRcluster();
      SRc.Algorithm = obj.Algorithm;
      SRc.E = obj.E;
      SRc.minPts = obj.MinPts;
      SRc.PlotFigures = false;
      SRc.ShrinkFactor = obj.ShrinkFactor;
   end

   % --------------------------------------------------------------------------

   obj.Lines = true;

   % ==========================================================================

   RT = ROITools();
   if obj.ROIs
      % Choose some ROIs to analyze the data over.
      if ~isempty(obj.ROI)
         RT.ROI_sizes = [obj.ROI(3), obj.ROI(4)];
      elseif ~isempty(obj.ROI_size)
         RT.ROI_sizes = [obj.ROI_size, obj.ROI_size];
      else
         error('ROI sizes not defined!');
      end
      RT.pixel2nm = obj.Pixel2nm;

      file_display = regexprep(basefile, '_', '\\_');
      [n_ROIs, RoI, Sigma_Reg] = RT.getROI({SRDR1, SRDR2}, file_display);
      saveas(gcf, sprintf('%s_ROIs.fig', base));
      if ~isempty(obj.Fig_ext)
         print(sprintf('%s_ROIs.%s', base, obj.Fig_ext), ['-d', obj.Fig_ext]);
      end
      close
   else
      % Analyze the data over the entire region.
      n_ROIs = 1;

      Sigma_Reg{1} = [10, 10];
      Sigma_Reg{2} = [10, 10];

      [XY1, XY_STD1, Sigma_Reg1] = RT.import_XY(SRDR1, obj.Pixel2nm);
      [XY2, XY_STD2, Sigma_Reg2] = RT.import_XY(SRDR2, obj.Pixel2nm);
      RoI{1}.X{1} = XY1(:, 1);
      RoI{1}.Y{1} = XY1(:, 2);
      RoI{1}.X{2} = XY2(:, 1);
      RoI{1}.Y{2} = XY2(:, 2);
      if ~isempty(XY_STD1)
         RoI{1}.X_STD{1} = XY_STD1(:, 1);
         RoI{1}.Y_STD{1} = XY_STD1(:, 2);
      else
         RoI{1}.X_STD{1} = zeros(size(XY1, 1), 1);
         RoI{1}.Y_STD{1} = zeros(size(XY1, 1), 1);
      end
      if ~isempty(XY_STD2)
         RoI{1}.X_STD{2} = XY_STD2(:, 1);
         RoI{1}.Y_STD{2} = XY_STD2(:, 2);
      else
         RoI{1}.X_STD{2} = zeros(size(XY2, 1), 1);
         RoI{1}.Y_STD{2} = zeros(size(XY2, 1), 1);
      end

      RoI{1}.ROI = [min([RoI{1}.X{1}; RoI{1}.X{2}]), ...
                    max([RoI{1}.X{1}; RoI{1}.X{2}]), ...
                    min([RoI{1}.Y{1}; RoI{1}.Y{2}]), ...
                    max([RoI{1}.Y{1}; RoI{1}.Y{2}])];
   end
   save(sprintf('%s_ROIs.mat', base), 'n_ROIs', 'RoI', 'Sigma_Reg');

   % Perform H-SET (optionally) and pair correlation over each ROI in turn,
   % combining the results at the end.
   ROIs_pre  = cell(n_ROIs, 1);
   ROIs_post = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      fprintf('ROI %d\n\n', i);

      obj.ROI = [RoI{i}.ROI(1), RoI{i}.ROI(3), ...
                 RoI{i}.ROI(2) - RoI{i}.ROI(1), RoI{i}.ROI(4) - RoI{i}.ROI(3)];

      XY1    = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      XY2    = [ RoI{i}.X{2}, RoI{i}.Y{2} ];
      sigma1 = [ RoI{i}.X_STD{1}, RoI{i}.Y_STD{1} ];
      sigma2 = [ RoI{i}.X_STD{2}, RoI{i}.Y_STD{2} ];

      base_file_pre  = sprintf('%s_ROI%d_pre',  basefile, i);
      % --- BEGIN Auto --------------------------------------------------------
      if obj.Auto
         base_file_pre1  = sprintf('%s_L1_ROI%d_pre',  basefile, i);
         base_file_pre2  = sprintf('%s_L2_ROI%d_pre',  basefile, i);
      end
      % --- END Auto ----------------------------------------------------------
      if obj.HSET
         base_file_post = sprintf('%s_ROI%d_post', basefile, i);
         % --- BEGIN Auto -----------------------------------------------------
         if obj.Auto
            base_file_post1 = sprintf('%s_L1_ROI%d_post', basefile, i);
            base_file_post2 = sprintf('%s_L2_ROI%d_post', basefile, i);
         end
         % --- END Auto -------------------------------------------------------
         base1 = sprintf('%s_L1_ROI%d', base, i);
         base2 = sprintf('%s_L2_ROI%d', base, i);
      end

      ROIs_pre{i}.ROI = RoI{i}.ROI;
      ROIs_pre{i}.XY1 = XY1;
      ROIs_pre{i}.XY2 = XY2;

      results_pccA = obj.pair_correlation(obj.Hist_bin_size, base_file_pre, ...
                                          XY1, XY2)
      % --- BEGIN Auto --------------------------------------------------------
      if obj.Auto
         results_pac1A = obj.pair_correlation(obj.Hist_bin_size, ...
                                              base_file_pre1, XY1)
         results_pac2A = obj.pair_correlation(obj.Hist_bin_size, ...
                                              base_file_pre2, XY2)
      end
      % --- END Auto ----------------------------------------------------------
      if obj.Veatch
         results_VpccA = obj.pair_correlation_Veatch(XY1, XY2, ...
                            obj.Hist_bin_size, base_file_pre, 'cross');
         % --- BEGIN Auto -----------------------------------------------------
         if obj.Auto
            results_Vpac1A = obj.pair_correlation_Veatch(XY1, [], ...
                                obj.Hist_bin_size, base_file_pre1, 'auto');
            results_Vpac2A = obj.pair_correlation_Veatch(XY2, [], ...
                                obj.Hist_bin_size, base_file_pre2, 'auto');
         end
         % --- END Auto -------------------------------------------------------
      end

      if obj.HSET
         [XY1_SR, sigma1_SR, combined1] = ...
            SRc.clusterSR(XY1, sigma1, Sigma_Reg{1});
         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            SRclusterFig1 = SRc.plotSRclusters();
            SRclusterFig1.Visible = 'on';
            saveas(SRclusterFig1, sprintf('%s_%s.fig', base1, 'SRcluster'));
            close(SRclusterFig1);
            SRcollapseFig1 = SRc.plotSRcollapse();
            SRcollapseFig1.Visible = 'on';
            saveas(SRcollapseFig1, sprintf('%s_%s.fig', base1, 'SRcollapse'));
            close(SRcollapseFig1);
            [results1, analysisFigs1] = SRc.analyzeSRclusters();
            for j = 1 : length(analysisFigs1)
               analysisFigs1{j}.Visible = 'on';
               saveas(analysisFigs1{j}, ...
                      sprintf('%s_%s%1d.fig', base1, 'SR', j));
               close(analysisFigs1{j});
            end
         end
         % --- END Verbose ----------------------------------------------------

         [XY2_SR, sigma2_SR, combined2] = ...
            SRc.clusterSR(XY2, sigma2, Sigma_Reg{2});
         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            SRclusterFig2 = SRc.plotSRclusters();
            SRclusterFig2.Visible = 'on';
            saveas(SRclusterFig2, sprintf('%s_%s.fig', base2, 'SRcluster'));
            close(SRclusterFig2);
            SRcollapseFig2 = SRc.plotSRcollapse();
            SRcollapseFig2.Visible = 'on';
            saveas(SRcollapseFig2, sprintf('%s_%s.fig', base2, 'SRcollapse'));
            close(SRcollapseFig2);
            [results2, analysisFigs2] = SRc.analyzeSRclusters();
            for j = 1 : length(analysisFigs2)
               analysisFigs2{j}.Visible = 'on';
               saveas(analysisFigs2{j}, ...
                      sprintf('%s_%s%1d.fig', base2, 'SR', j));
               close(analysisFigs2{j});
            end
         end
         % --- END Verbose ----------------------------------------------------

         ROIs_post{i}.ROI = RoI{i}.ROI;
         ROIs_post{i}.XY1 = XY1_SR;
         ROIs_post{i}.XY2 = XY2_SR;

         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            save(sprintf('%s_ROI%d_XYsigma.mat', base, i), ...
                 'XY1', 'sigma1', 'XY2', 'sigma2',    ...
                 'XY1_SR', 'sigma1_SR', 'XY2_SR', 'sigma2_SR');
         end
         % --- END Verbose ----------------------------------------------------

         results_pccB = obj.pair_correlation(obj.Hist_bin_size, ...
                                             base_file_post, XY1_SR, XY2_SR)
         % --- BEGIN Auto -----------------------------------------------------
         if obj.Auto
            results_pac1B = obj.pair_correlation(obj.Hist_bin_size, ...
                                                 base_file_post1, XY1_SR)
            results_pac2B = obj.pair_correlation(obj.Hist_bin_size, ...
                                                 base_file_post2, XY2_SR)
         end
         % --- END Auto -------------------------------------------------------
         if obj.Veatch
            results_VpccB = obj.pair_correlation_Veatch(XY1_SR, XY2_SR, ...
                               obj.Hist_bin_size, base_file_post, 'cross');
            % --- BEGIN Auto --------------------------------------------------
            if obj.Auto
               results_Vpac1B = obj.pair_correlation_Veatch(XY1_SR, [], ...
                                   obj.Hist_bin_size, base_file_post1, 'auto');
               results_Vpac2B = obj.pair_correlation_Veatch(XY2_SR, [], ...
                                   obj.Hist_bin_size, base_file_post2, 'auto');
            end
            % --- END Auto ----------------------------------------------------
         end

         figure();
         hold on
         plot(XY1_SR(:, 1), XY1_SR(:, 2), 'k.', 'MarkerSize', 10);
         plot(XY2_SR(:, 1), XY2_SR(:, 2), 'b.', 'MarkerSize', 10);
         title(sprintf('Collapsed Two Color Plot ROI %d', i));
         xlabel('x (nm)');
         ylabel('y (nm)');
         hold off
         if ~isempty(obj.Fig_ext)
            print(sprintf('%s_ROI%d_collapsed.%s', base, i, obj.Fig_ext), ...
                  ['-d', obj.Fig_ext]);
         else
            print(sprintf('%s_ROI%d_collapsed.pdf', base, i), '-dpdf');
         end
         close(gcf);

         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            if obj.Auto && obj.Veatch
               save(sprintf('%s_ROI%d_results.mat', base, i),            ...
                    'results_pccA',  'results_pac1A',  'results_pac2A',  ...
                    'results_VpccA', 'results_Vpac1A', 'results_Vpac2A', ...
                    'results_pccB',  'results_pac1B',  'results_pac2B',  ...
                    'results_VpccB', 'results_Vpac1B', 'results_Vpac2B', ...
                    'results1', 'results2', 'obj', 'SRc');
            else
               save(sprintf('%s_ROI%d_results.mat', base, i), ...
                    'results_pccA',  'results_pccB',          ...
                    'results1', 'results2', 'obj', 'SRc');
            end
         end
         % --- END Verbose ----------------------------------------------------
      end
   end

   if obj.ROIs
      fprintf('ROI combined\n\n');

      base_pre = sprintf('%s_pre',  basefile);
      resultsRC_pccA = ...
         obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_pre, ...
                                          2, n_ROIs, ROIs_pre)
      % --- BEGIN Auto --------------------------------------------------------
      if obj.Auto
         base_pre1 = sprintf('%s_L1_pre', basefile);
         base_pre2 = sprintf('%s_L2_pre', basefile);
         resultsRC_pac1A = ...
            obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_pre1, ...
                                             1, n_ROIs, ROIs_pre, 1)
         resultsRC_pac2A = ...
            obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_pre2, ...
                                             1, n_ROIs, ROIs_pre, 2)
      end
      % --- END Auto ----------------------------------------------------------
      if obj.HSET
         base_post = sprintf('%s_post', basefile);
         resultsRC_pccB = ...
            obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_post, ...
                                             2, n_ROIs, ROIs_post)
         % --- BEGIN Auto -----------------------------------------------------
         if obj.Auto
            base_post1 = sprintf('%s_L1_post', basefile);
            base_post2 = sprintf('%s_L2_post', basefile);
            resultsRC_pac1B = ...
               obj.pair_correlation_ROIcombined(obj.Hist_bin_size,     ...
                                                base_post1, 1, n_ROIs, ...
                                                ROIs_post, 1)
            resultsRC_pac2B = ...
               obj.pair_correlation_ROIcombined(obj.Hist_bin_size,     ...
                                                base_post2, 1, n_ROIs, ...
                                                ROIs_post, 2)
            save(sprintf('%s_results.mat', base),      ...
                 'resultsRC_pccA', 'resultsRC_pccB',   ...
                 'resultsRC_pac1A', 'resultsRC_pac2A', ...
                 'resultsRC_pac1B', 'resultsRC_pac2B');
         % --- END Auto -------------------------------------------------------
         else
            save(sprintf('%s_results.mat', base), ...
                 'resultsRC_pccA', 'resultsRC_pccB');
         end
      else
         % --- BEGIN Auto -----------------------------------------------------
         if obj.Auto
            save(sprintf('%s_results.mat', base), 'resultsRC_pccA', ...
                 'resultsRC_pac1A', 'resultsRC_pac2A');
         % --- END Auto -------------------------------------------------------
         else
            save(sprintf('%s_results.mat', base), 'resultsRC_pccA');
         end
      end
   end

end

% -----------------------------------------------------------------------------

function PCC_SR(obj, SRDR1, SRDR2, basefile)
% Perform pair correlation before and after H-SET collapse on the given input
% SRD.Results datasets.
%
% Cross-correlation is performed.
%
% Use the mouse to select ROIs (regions of interest):
%    left click  chooses the center of a fixed size (x_size x y_size) region
%    right click chooses an adjustable rectangular size region
%    key press:
%       backspace or delete   deletes the previous region
%       anything else         terminates selection
%
% Inputs:
%    SRDR1      SRD.Results dataset label 1 
%    SRDR2      SRD.Results dataset label 2
%    basefile   base filename to which various text will be appended to
%               identify the results

   % Results directory.
   if ~isdir(obj.Results)
      mkdir(obj.Results);
   end
   base = fullfile(obj.Results, basefile);

   % ==========================================================================

   if obj.HSET
      SRc = SRcluster();
      SRc.Algorithm = obj.Algorithm;
      SRc.E = obj.E;
      SRc.minPts = obj.MinPts;
      SRc.PlotFigures = false;
      SRc.ShrinkFactor = obj.ShrinkFactor;
   end

   % --------------------------------------------------------------------------

   obj.Lines = true;

   % ==========================================================================

   RT = ROITools();
   if obj.ROIs
      % Choose some ROIs to analyze the data over.
      if ~isempty(obj.ROI)
         RT.ROI_sizes = [obj.ROI(3), obj.ROI(4)];
      elseif ~isempty(obj.ROI_size)
         RT.ROI_sizes = [obj.ROI_size, obj.ROI_size];
      else
         error('ROI sizes not defined!');
      end
      RT.pixel2nm = obj.Pixel2nm;

      file_display = regexprep(basefile, '_', '\\_');
      [n_ROIs, RoI, Sigma_Reg] = RT.getROI({SRDR1, SRDR2}, file_display);
      saveas(gcf, sprintf('%s_ROIs.fig', base));
      if ~isempty(obj.Fig_ext)
         print(sprintf('%s_ROIs.%s', base, obj.Fig_ext), ['-d', obj.Fig_ext]);
      end
      close
   else
      % Analyze the data over the entire region.
      n_ROIs = 1;

      Sigma_Reg{1} = [10, 10];
      Sigma_Reg{2} = [10, 10];

      [XY1, XY_STD1, Sigma_Reg1] = RT.import_XY(SRDR1, obj.Pixel2nm);
      [XY2, XY_STD2, Sigma_Reg2] = RT.import_XY(SRDR2, obj.Pixel2nm);
      RoI{1}.X{1} = XY1(:, 1);
      RoI{1}.Y{1} = XY1(:, 2);
      RoI{1}.X{2} = XY2(:, 1);
      RoI{1}.Y{2} = XY2(:, 2);
      if ~isempty(XY_STD1)
         RoI{1}.X_STD{1} = XY_STD1(:, 1);
         RoI{1}.Y_STD{1} = XY_STD1(:, 2);
      else
         RoI{1}.X_STD{1} = zeros(size(XY1, 1), 1);
         RoI{1}.Y_STD{1} = zeros(size(XY1, 1), 1);
      end
      if ~isempty(XY_STD2)
         RoI{1}.X_STD{2} = XY_STD2(:, 1);
         RoI{1}.Y_STD{2} = XY_STD2(:, 2);
      else
         RoI{1}.X_STD{2} = zeros(size(XY2, 1), 1);
         RoI{1}.Y_STD{2} = zeros(size(XY2, 1), 1);
      end

      RoI{1}.ROI = [min([RoI{1}.X{1}; RoI{1}.X{2}]), ...
                    max([RoI{1}.X{1}; RoI{1}.X{2}]), ...
                    min([RoI{1}.Y{1}; RoI{1}.Y{2}]), ...
                    max([RoI{1}.Y{1}; RoI{1}.Y{2}])];
   end
   save(sprintf('%s_ROIs.mat', base), 'n_ROIs', 'RoI', 'Sigma_Reg');

   % Perform H-SET (optionally) and pair correlation over each ROI in turn,
   % combining the results at the end.
   ROIs_pre  = cell(n_ROIs, 1);
   ROIs_post = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      fprintf('ROI %d\n\n', i);

      obj.ROI = [RoI{i}.ROI(1), RoI{i}.ROI(3), ...
                 RoI{i}.ROI(2) - RoI{i}.ROI(1), RoI{i}.ROI(4) - RoI{i}.ROI(3)];

      XY1    = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      XY2    = [ RoI{i}.X{2}, RoI{i}.Y{2} ];
      sigma1 = [ RoI{i}.X_STD{1}, RoI{i}.Y_STD{1} ];
      sigma2 = [ RoI{i}.X_STD{2}, RoI{i}.Y_STD{2} ];

      base_file_pre  = sprintf('%s_ROI%d_pre',  basefile, i);
      if obj.HSET
         base_file_post = sprintf('%s_ROI%d_post', basefile, i);
         base1 = sprintf('%s_L1_ROI%d', base, i);
         base2 = sprintf('%s_L2_ROI%d', base, i);
      end

      ROIs_pre{i}.ROI = RoI{i}.ROI;
      ROIs_pre{i}.XY1 = XY1;
      ROIs_pre{i}.XY2 = XY2;

      results_pccA = obj.pair_correlation(obj.Hist_bin_size, base_file_pre, ...
                                          XY1, XY2)
      if obj.Veatch
         results_VpccA = obj.pair_correlation_Veatch(XY1, XY2, ...
                            obj.Hist_bin_size, base_file_pre, 'cross');
      end

      if obj.HSET
         [XY1_SR, sigma1_SR, combined1] = ...
            SRc.clusterSR(XY1, sigma1, Sigma_Reg{1});
         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            SRclusterFig1 = SRc.plotSRclusters();
            SRclusterFig1.Visible = 'on';
            saveas(SRclusterFig1, sprintf('%s_%s.fig', base1, 'SRcluster'));
            close(SRclusterFig1);
            SRcollapseFig1 = SRc.plotSRcollapse();
            SRcollapseFig1.Visible = 'on';
            saveas(SRcollapseFig1, sprintf('%s_%s.fig', base1, 'SRcollapse'));
            close(SRcollapseFig1);
            [results1, analysisFigs1] = SRc.analyzeSRclusters();
            for j = 1 : length(analysisFigs1)
               analysisFigs1{j}.Visible = 'on';
               saveas(analysisFigs1{j}, ...
                      sprintf('%s_%s%1d.fig', base1, 'SR', j));
               close(analysisFigs1{j});
            end
         end
         % --- END Verbose ----------------------------------------------------

         [XY2_SR, sigma2_SR, combined2] = ...
            SRc.clusterSR(XY2, sigma2, Sigma_Reg{2});
         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            SRclusterFig2 = SRc.plotSRclusters();
            SRclusterFig2.Visible = 'on';
            saveas(SRclusterFig2, sprintf('%s_%s.fig', base2, 'SRcluster'));
            close(SRclusterFig2);
            SRcollapseFig2 = SRc.plotSRcollapse();
            SRcollapseFig2.Visible = 'on';
            saveas(SRcollapseFig2, sprintf('%s_%s.fig', base2, 'SRcollapse'));
            close(SRcollapseFig2);
            [results2, analysisFigs2] = SRc.analyzeSRclusters();
            for j = 1 : length(analysisFigs2)
               analysisFigs2{j}.Visible = 'on';
               saveas(analysisFigs2{j}, ...
                      sprintf('%s_%s%1d.fig', base2, 'SR', j));
               close(analysisFigs2{j});
            end
         end
         % --- END Verbose ----------------------------------------------------

         ROIs_post{i}.ROI = RoI{i}.ROI;
         ROIs_post{i}.XY1 = XY1_SR;
         ROIs_post{i}.XY2 = XY2_SR;

         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            save(sprintf('%s_ROI%d_XYsigma.mat', base, i), ...
                 'XY1', 'sigma1', 'XY2', 'sigma2',    ...
                 'XY1_SR', 'sigma1_SR', 'XY2_SR', 'sigma2_SR');
         end
         % --- END Verbose ----------------------------------------------------

         results_pccB = obj.pair_correlation(obj.Hist_bin_size, ...
                           base_file_post, XY1_SR, XY2_SR)
         if obj.Veatch
            results_VpccB = obj.pair_correlation_Veatch(XY1_SR, XY2_SR, ...
                               obj.Hist_bin_size, base_file_post, 'cross');
         end

         figure();
         hold on
         plot(XY1_SR(:, 1), XY1_SR(:, 2), 'k.', 'MarkerSize', 10);
         plot(XY2_SR(:, 1), XY2_SR(:, 2), 'b.', 'MarkerSize', 10);
         title(sprintf('Collapsed Two Color Plot ROI %d', i));
         xlabel('x (nm)');
         ylabel('y (nm)');
         hold off
         if ~isempty(obj.Fig_ext)
            print(sprintf('%s_ROI%d_collapsed.%s', base, i, obj.Fig_ext), ...
                  ['-d', obj.Fig_ext]);
         else
            print(sprintf('%s_ROI%d_collapsed.pdf', base, i), '-dpdf');
         end
         close(gcf);

         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            if obj.Veatch
               save(sprintf('%s_ROI%d_results.mat', base, i), ...
                    'results_pccA', 'results_VpccA',          ...
                    'results_pccB', 'results_VpccB',          ...
                    'results1', 'results2', 'obj', 'SRc');
            else
               save(sprintf('%s_ROI%d_results.mat', base, i), ...
                    'results_pccA', 'results_pccB',           ...
                    'results1', 'results2', 'obj', 'SRc');
            end
         end
         % --- END Verbose ----------------------------------------------------
      end
   end

   if obj.ROIs
      fprintf('ROI combined\n\n');

      base_pre  = sprintf('%s_pre',  basefile);
      resultsRC_pccA = ...
         obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_pre, ...
                                          2, n_ROIs, ROIs_pre)
      if obj.HSET
         base_post = sprintf('%s_post', basefile);
         resultsRC_pccB = ...
            obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_post, ...
                                             2, n_ROIs, ROIs_post)
         save(sprintf('%s_results.mat', base), ...
              'resultsRC_pccA', 'resultsRC_pccB');
      else
         save(sprintf('%s_results.mat', base), 'resultsRC_pccA');
      end
   end

end

% -----------------------------------------------------------------------------

function PAC_SR(obj, SRDR1, basefile)
% Perform pair correlation before and after H-SET collapse on the given input
% SRD.Results dataset.
%
% Auto-correlation is performed.
%
% Use the mouse to select ROIs (regions of interest):
%    left click  chooses the center of a fixed size (x_size x y_size) region
%    right click chooses an adjustable rectangular size region
%    key press:
%       backspace or delete   deletes the previous region
%       anything else         terminates selection
%
% Inputs:
%    SRDR1      SRD.Results dataset label 1 
%    SRDR2      SRD.Results dataset label 2
%    basefile   base filename to which various text will be appended to
%               identify the results

   % Results directory.
   if ~isdir(obj.Results)
      mkdir(obj.Results);
   end
   base = fullfile(obj.Results, basefile);

   % ==========================================================================

   if obj.HSET
      SRc = SRcluster();
      SRc.Algorithm = obj.Algorithm;
      SRc.E = obj.E;
      SRc.minPts = obj.MinPts;
      SRc.PlotFigures = false;
      SRc.ShrinkFactor = obj.ShrinkFactor;
   end

   % --------------------------------------------------------------------------

   obj.Lines = true;

   % ==========================================================================

   RT = ROITools();
   if obj.ROIs
      % Choose some ROIs to analyze the data over.
      if ~isempty(obj.ROI)
         RT.ROI_sizes = [obj.ROI(3), obj.ROI(4)];
      elseif ~isempty(obj.ROI_size)
         RT.ROI_sizes = [obj.ROI_size, obj.ROI_size];
      else
         error('ROI sizes not defined!');
      end
      RT.pixel2nm = obj.Pixel2nm;

      file_display = regexprep(basefile, '_', '\\_');
      [n_ROIs, RoI, Sigma_Reg] = RT.getROI({SRDR1}, file_display);
      saveas(gcf, sprintf('%s_ROIs.fig', base));
      if ~isempty(obj.Fig_ext)
         print(sprintf('%s_ROIs.%s', base, obj.Fig_ext), ['-d', obj.Fig_ext]);
      end
      close
   else
      % Analyze the data over the entire region.
      n_ROIs = 1;

      Sigma_Reg{1} = [10, 10];

      [XY1, XY_STD1, Sigma_Reg1] = RT.import_XY(SRDR1, obj.Pixel2nm);
      RoI{1}.X{1} = XY1(:, 1);
      RoI{1}.Y{1} = XY1(:, 2);
      if ~isempty(XY_STD1)
         RoI{1}.X_STD{1} = XY_STD1(:, 1);
         RoI{1}.Y_STD{1} = XY_STD1(:, 2);
      else
         RoI{1}.X_STD{1} = zeros(size(XY1, 1), 1);
         RoI{1}.Y_STD{1} = zeros(size(XY1, 1), 1);
      end

      RoI{1}.ROI = [min(RoI{1}.X{1}), max(RoI{1}.X{1}), ...
                    min(RoI{1}.Y{1}), max(RoI{1}.Y{1})];
   end
   save(sprintf('%s_ROIs.mat', base), 'n_ROIs', 'RoI', 'Sigma_Reg');

   % Perform H-SET (optionally) and pair correlation over each ROI in turn,
   % combining the results at the end.
   ROIs_pre  = cell(n_ROIs, 1);
   ROIs_post = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      fprintf('ROI %d\n\n', i);

      obj.ROI = [RoI{i}.ROI(1), RoI{i}.ROI(3), ...
                 RoI{i}.ROI(2) - RoI{i}.ROI(1), RoI{i}.ROI(4) - RoI{i}.ROI(3)];

      XY1    = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      sigma1 = [ RoI{i}.X_STD{1}, RoI{i}.Y_STD{1} ];

      base_file_pre1  = sprintf('%s_ROI%d_pre',  basefile, i);
      if obj.HSET
         base_file_post = sprintf('%s_ROI%d_post', basefile, i);
         base_file_post1 = sprintf('%s_ROI%d_post', basefile, i);
         base1 = sprintf('%s_ROI%d', base, i);
      end

      ROIs_pre{i}.ROI = RoI{i}.ROI;
      ROIs_pre{i}.XY1 = XY1;

      results_pac1A = obj.pair_correlation(obj.Hist_bin_size, ...
                                           base_file_pre1, XY1)
      if obj.Veatch
         results_Vpac1A = obj.pair_correlation_Veatch(XY1, [], ...
                             obj.Hist_bin_size, base_file_pre1, 'auto');
      end

      if obj.HSET
         [XY1_SR, sigma1_SR, combined1] = ...
            SRc.clusterSR(XY1, sigma1, Sigma_Reg{1});
         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            SRclusterFig1 = SRc.plotSRclusters();
            SRclusterFig1.Visible = 'on';
            saveas(SRclusterFig1, sprintf('%s_%s.fig', base1, 'SRcluster'));
            close(SRclusterFig1);
            SRcollapseFig1 = SRc.plotSRcollapse();
            SRcollapseFig1.Visible = 'on';
            saveas(SRcollapseFig1, sprintf('%s_%s.fig', base1, 'SRcollapse'));
            close(SRcollapseFig1);
            [results1, analysisFigs1] = SRc.analyzeSRclusters();
            for j = 1 : length(analysisFigs1)
               analysisFigs1{j}.Visible = 'on';
               saveas(analysisFigs1{j}, ...
                      sprintf('%s_%s%1d.fig', base1, 'SR', j));
               close(analysisFigs1{j});
            end
         end
         % --- END Verbose ----------------------------------------------------

         ROIs_post{i}.ROI = RoI{i}.ROI;
         ROIs_post{i}.XY1 = XY1_SR;

         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            save(sprintf('%s_ROI%d_XYsigma.mat', base, i), ...
                 'XY1', 'sigma1', 'XY1_SR', 'sigma1_SR');
         end
         % --- END Verbose ----------------------------------------------------

         results_pac1B = obj.pair_correlation(obj.Hist_bin_size, ...
                                              base_file_post1, XY1_SR)
         if obj.Veatch
            results_Vpac1B = ...
               obj.pair_correlation_Veatch(XY1_SR, [], obj.Hist_bin_size, ...
                                           base_file_post1, 'auto');
         end

         figure();
         hold on
         plot(XY1_SR(:, 1), XY1_SR(:, 2), 'k.', 'MarkerSize', 10);
         title(sprintf('Collapsed Two Color Plot ROI %d', i));
         xlabel('x (nm)');
         ylabel('y (nm)');
         hold off
         if ~isempty(obj.Fig_ext)
            print(sprintf('%s_ROI%d_collapsed.%s', base, i, obj.Fig_ext), ...
                  ['-d', obj.Fig_ext]);
         else
            print(sprintf('%s_ROI%d_collapsed.pdf', base, i), '-dpdf');
         end
         close(gcf);

         % --- BEGIN Verbose --------------------------------------------------
         if obj.Verbose
            if obj.Veatch
               save(sprintf('%s_ROI%d_results.mat', base, i), ...
                    'results_pac1A', 'results_Vpac1A',        ...
                    'results_pac1B', 'results_Vpac1B',        ...
                    'results1', 'obj', 'SRc');
            else
               save(sprintf('%s_ROI%d_results.mat', base, i), ...
                    'results1', 'obj', 'SRc');
            end
         end
         % --- END Verbose ----------------------------------------------------
      end
   end

   if obj.ROIs
      fprintf('ROI combined\n\n');

      base_pre  = sprintf('%s_pre',  basefile);
      resultsRC_pacA = ...
         obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_pre, ...
                                          1, n_ROIs, ROIs_pre)
      if obj.HSET
         base_post = sprintf('%s_post', basefile);
         resultsRC_pacB = ...
            obj.pair_correlation_ROIcombined(obj.Hist_bin_size, base_post, ...
                                             1, n_ROIs, ROIs_post)
         save(sprintf('%s_results.mat', base), ...
              'resultsRC_pacA', 'resultsRC_pacB');
      else
         save(sprintf('%s_results.mat', base), 'resultsRC_pacA');
      end
   end

end

% =============================================================================

function results = pair_correlation_ROIcombined(obj, hist_bin_size,  ...
                                                base_name, n_labels, ...
                                                n_ROIs, ROIs, label_num)
% Combine ROIs while performing pair correlation.
% Modified from code originally written by Carolyn Pehlke.
%
% Inputs:
%    hist_bin_size   histogram bin size sometims known as pixel size (nm)
%    base_name       descriptive name for the results files
%    n_labels        number of different labeled particles:
%                    1 -> auto-correlation, 2 -> cross_correlation
%    n_ROIs          number of ROIs to combine
%    ROIs            contains (x, y) coordinates of the two datasets (nm)
%    label_num       [OPTIONAL] if n_labels = 1, then whether XY1 or XY2 should
%                    be used for the coordinates contained in ROIs
% Output:
%    results         structure containing various results from the algorithm

   if n_labels == 2
      corr_type = 'C';   % cross-correlation
   else
      corr_type = 'A';   % auto-correlation
   end

   % Add in the DERIVESTsuite path.
   %addpath('DERIVESTsuite', filesep);

   ROI_size = zeros(n_ROIs, 1);
   IM1 = cell(n_ROIs, 1);
   if corr_type == 'C'
      IM2 = cell(n_ROIs, 1);
   end
   for j = 1 : n_ROIs
      x_min = ROIs{j}.ROI(1);
      x_max = ROIs{j}.ROI(2);
      y_min = ROIs{j}.ROI(3);
      y_max = ROIs{j}.ROI(4);
      ROI_size(j) = min(x_max - x_min, y_max - y_min);

      % Compute the number of pixels in x and y.
      imszX = round((x_max - x_min) / hist_bin_size);
      imszY = round((y_max - y_min) / hist_bin_size);
      % Create a blank image.
      im1 = zeros(imszX, imszY);
      if corr_type == 'C'
         im2 = zeros(imszX, imszY);
      end
      % Convert (x, y) coordinates into pixel units.
      if ~exist('label_num', 'var')
         label_num = 1;
      end
      if label_num == 1
         x1 = round((ROIs{j}.XY1(:, 1) - x_min) / hist_bin_size) + 1;
         y1 = round((ROIs{j}.XY1(:, 2) - y_min) / hist_bin_size) + 1;
      elseif label_num == 2
         x1 = round((ROIs{j}.XY2(:, 1) - x_min) / hist_bin_size) + 1;
         y1 = round((ROIs{j}.XY2(:, 2) - y_min) / hist_bin_size) + 1;
      else
         error('Invalid label_num: %d!', label_num);
      end
      if corr_type == 'C'
         x2 = round((ROIs{j}.XY2(:, 1) - x_min) / hist_bin_size) + 1;
         y2 = round((ROIs{j}.XY2(:, 2) - y_min) / hist_bin_size) + 1;
      end
      % Get the pixels within the image size.
      mask1 = (x1 > 0) & (y1 > 0) & (x1 <= imszX) & (y1 <= imszY);
      x1 = x1(mask1);
      y1 = y1(mask1);
      if corr_type == 'C'
         mask2 = (x2 > 0) & (y2 > 0) & (x2 <= imszX) & (y2 <= imszY);
         x2 = x2(mask2);
         y2 = y2(mask2);
      end
      % Make a histogram image.
      for i = 1 : size(x1, 1)
         im1(x1(i), y1(i)) = im1(x1(i), y1(i)) + 1;
      end
      if corr_type == 'C'
         for i = 1 : size(x2, 1)
            im2(x2(i), y2(i)) = im2(x2(i), y2(i)) + 1;
         end
      end
      IM1{j} = im1;
      if corr_type == 'C'
         IM2{j} = im2;
      end
   end
   % Establish rmax as half the size of the ROI in pixels.
   rmax = round(min(ROI_size) / (2 * hist_bin_size));

   if corr_type == 'C'
      % Pair crosscorrelation using Veatch method.
      [G, r, g, dg, rmax] = PairCorr.get_corr(n_ROIs, rmax, IM1, IM2);
   else
      % Pair autocorrelation using Veatch method.
      [G, r, g, dg, rmax] = PairCorr.get_corr(n_ROIs, rmax, IM1);
   end

   % Convert back to nm
   r_nm = r * hist_bin_size;

   if corr_type == 'C'
      i1 = 0;   i2 = 0;   % intensity sums
      p1 = 0;   p2 = 0;   % pixel sums
      for i = 1 : n_ROIs
         i1 = i1 + sum(sum(IM1{i}));
         i2 = i2 + sum(sum(IM2{i}));
         p1 = p1 + prod(size(IM1{i}));
         p2 = p2 + prod(size(IM2{i}));
      end
      rho1 = i1 / p1;
      rho2 = i2 / p2;
      %rho1 = mean(mean(im1));
      %rho2 = mean(mean(im2));
      rho = sqrt(rho1 * rho2);
      %paircorr = [{[im1, im2]}, {[r', g', dg']}, {rho}, {G}];
   else
      i1 = 0;   % intensity sum
      p1 = 0;   % pixel sum
      for i = 1 : n_ROIs
         i1 = i1 + sum(sum(IM1{i}));
         p1 = p1 + prod(size(IM1{i}));
      end
      rho = i1 / p1;
      %rho = mean(mean(im1));
      %paircorr = [{im1}, {[r', g', dg']}, {rho}, {G}];
   end

   [estimates, errors, model] = PairCorr.PC_GaussFit(r', g', rmax, rho);
   %  PairCorr.PC_GaussFit(paircorr{2}(:,1), paircorr{2}(:,2), rmax, ...
   %                       paircorr{3});
   estimates = abs(estimates);

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:});
   hold on
   if corr_type == 'C'
      name = fullfile(obj.Results, [base_name, '_crosscorr']);
   else
      name = fullfile(obj.Results, [base_name, '_autocorr']);
   end
   if obj.Lines
      plot(r_nm(2:end),g(2:end),'k.-','MarkerSize',20,'LineWidth',2)
      %plot(r_nm(2:end),paircorr{2}(2:end,2),'k.-', ...
      %     'MarkerSize',20,'LineWidth',2)
   else
      plot(r_nm(2:end),g(2:end),'k.','MarkerSize',10,'LineWidth',2)
      %plot(r_nm(2:end),paircorr{2}(2:end,2),'k.','MarkerSize',10,'LineWidth',2)
   end
   plot(r_nm(2:end),ones(1, numel(r) - 1),'b:','LineWidth',3)
   plot(r_nm(2:end),model(2:end),'--r','LineWidth',3)
   if corr_type == 'C'
      flegend{1} = 'Cross-Correlation';
   else
      flegend{1} = 'Auto-Correlation';
   end
   flegend{2} = 'g(r) Random';
   flegend{3} = 'Fit';
   axis tight
   if obj.Rmax_axis > 0
      xlim([0, obj.Rmax_axis]);
   end
   title(regexprep(base_name, '_', '\\_'));
   xlabel('r (nm)');
   ylabel('g(r)');
   legend(flegend, 'Location', 'Best');
   hold off

   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
      saveas(gcf, name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   Afound=estimates(1);
   s_d_found_pixels=estimates(2);
   Bfound=estimates(3);
   s_l_found_pixels=estimates(4);

   Afound_SE=errors(1);
   s_d_found_pixels_SE=errors(2);
   Bfound_SE=errors(3);
   s_l_found_pixels_SE=errors(4);

   s_d_found=s_d_found_pixels*hist_bin_size;
   s_l_found=s_l_found_pixels*hist_bin_size;
   s_d_found_SE=s_d_found_pixels_SE*hist_bin_size;
   s_l_found_SE=s_l_found_pixels_SE*hist_bin_size;

   % display results
   if corr_type == 'C'
      fprintf('Pair Cross-correlation for %s:\n\n', base_name);
   else
      fprintf('Pair Auto-correlation for %s:\n\n', base_name);
   end

   fprintf('Objects per domain ');
   fprintf('found:\t%5.3g +/- %5.3g\n',Afound,Afound_SE);

   fprintf('Domain size sigma ');
   fprintf('found:\t%5.3g +/- %5.3g\n',s_d_found,s_d_found_SE);

   fprintf('Localizations per object ');
   fprintf('found:\t%5.3g +/- %5.3g\n',Bfound,Bfound_SE);

   fprintf('Localization precision sigma ');
   fprintf('found:\t%5.3g +/- %5.3g\n',s_l_found,s_l_found_SE);

   results.G  = G;
   results.r  = r_nm;
   results.g  = g;
   results.dg = dg;
   results.objs_per_domain          = Afound;
   results.objs_per_domain_SE       = Afound_SE;
   results.sigma_domain             = s_d_found;
   results.sigma_domain_SE          = s_d_found_SE;
   results.localizations_per_obj    = Bfound;
   results.localizations_per_obj_SE = Bfound_SE;
   results.sigma_localization       = s_l_found;
   results.sigma_localization_SE    = s_l_found_SE;
   results.model = model;

end

% -----------------------------------------------------------------------------

function results = pair_correlation(obj, hist_bin_size, base_name, XY1, XY2)
% Perform pair correlation.
% Modified from code originally written by Carolyn Pehlke.
%
% Input:
%    hist_bin_size    pixel size (nm)
%    base_name        descriptive name for the results files
%    XY1 and XY2      (x, y) coordinates of the two datasets (nm) [N x 2]
%                     XY2 is optional and if omitted requests auto-correlation
%                     rather than cross-correlation
% Output:
%    results         structure containing various results from the algorithm

   if exist('XY2', 'var')
      corr_type = 'C';   % cross-correlation
   else
      corr_type = 'A';   % auto-correlation
   end

   if ~isempty(obj.ROI)
      x_min = obj.ROI(1);
      y_min = obj.ROI(2);
      x_max = x_min + obj.ROI(3);
      y_max = y_min + obj.ROI(4);
      ROI_size = min(obj.ROI(3 : 4));
   else
      if corr_type == 'C'
         x_min = min([XY1(:, 1); XY2(:, 1)]);
         y_min = min([XY1(:, 2); XY2(:, 2)]);
      else
         x_min = min(XY1(:, 1));
         y_min = min(XY1(:, 2));
      end
      if ~isempty(obj.ROI_size)
         x_max = x_min + obj.ROI_size;
         y_max = y_min + obj.ROI_size;
         ROI_size = obj.ROI_size;
      else
         if corr_type == 'C'
            x_max = max([XY1(:, 1); XY2(:, 1)]);
            y_max = max([XY1(:, 2); XY2(:, 2)]);
         else
            x_max = max(XY1(:, 1));
            y_max = max(XY1(:, 2));
         end
         ROI_size = min(x_max - x_min, y_max - y_min);
      end
   end

   % Set to 1 to get a figure of g(r) with error bars.
   flag = 0;

   % Add in the DERIVESTsuite path.
   %addpath('DERIVESTsuite', filesep);

   % Compute the number of pixels in x and y.
   imszX = round((x_max - x_min) / hist_bin_size);
   imszY = round((y_max - y_min) / hist_bin_size);
   % Create a blank image.
   im1 = zeros(imszX, imszY);
   if corr_type == 'C'
      im2 = zeros(imszX, imszY);
   end
   % Convert (x, y) coordinates into pixel units.
   x1 = round((XY1(:, 1) - x_min) / hist_bin_size) + 1;
   y1 = round((XY1(:, 2) - y_min) / hist_bin_size) + 1;
   if corr_type == 'C'
      x2 = round((XY2(:, 1) - x_min) / hist_bin_size) + 1;
      y2 = round((XY2(:, 2) - y_min) / hist_bin_size) + 1;
   end
   % Get the pixels within the image size.
   mask1 = (x1 > 0) & (y1 > 0) & (x1 <= imszX) & (y1 <= imszY);
   x1 = x1(mask1);
   y1 = y1(mask1);
   if corr_type == 'C'
      mask2 = (x2 > 0) & (y2 > 0) & (x2 <= imszX) & (y2 <= imszY);
      x2 = x2(mask2);
      y2 = y2(mask2);
   end
   % Make a histogram image.
   for i = 1 : size(x1, 1)
      im1(x1(i), y1(i)) = im1(x1(i), y1(i)) + 1;
   end
   if corr_type == 'C'
      for i = 1 : size(x2, 1)
         im2(x2(i), y2(i)) = im2(x2(i), y2(i)) + 1;
      end
   end
   % Establish rmax as half the size of the ROI in pixels.
   rmax = round(ROI_size / (2 * hist_bin_size));
   if corr_type == 'C'
      % Pair crosscorrelation using Veatch method.
      [G, r, g, dg, ~, rmax] = ...
         PairCorr.get_crosscorr(im1, im2, ones(size(im1)), rmax, flag);
   else
      % Pair autocorrelation using Veatch method.
      [G, r, g, dg, ~, rmax] = ...
         PairCorr.get_autocorr(im1, ones(size(im1)), rmax, flag);
   end
   % Convert back to nm
   r_nm = r * hist_bin_size;

   if corr_type == 'C'
      rho1 = mean(mean(im1));
      rho2 = mean(mean(im2));
      rho = sqrt(rho1 * rho2);
      paircorr = [{[im1, im2]}, {[r', g', dg']}, {rho}, {G}];
   else
      rho = mean(mean(im1));
      paircorr = [{im1}, {[r', g', dg']}, {rho}, {G}];
   end

   [estimates, errors, model] = ...
      PairCorr.PC_GaussFit(paircorr{2}(:,1), paircorr{2}(:,2), rmax, ...
                           paircorr{3});
   estimates = abs(estimates);

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:});
   hold on
   if corr_type == 'C'
      name = fullfile(obj.Results, [base_name, '_crosscorr']);
   else
      name = fullfile(obj.Results, [base_name, '_autocorr']);
   end
   if obj.Lines
      plot(r_nm(2:end),paircorr{2}(2:end,2),'k.-', ...
           'MarkerSize',20,'LineWidth',2)
   else
      plot(r_nm(2:end),paircorr{2}(2:end,2),'k.','MarkerSize',10,'LineWidth',2)
   end
   plot(r_nm(2:end),ones(1, numel(r) - 1),'b:','LineWidth',3)
   plot(r_nm(2:end),model(2:end),'--r','LineWidth',3)
   if corr_type == 'C'
      flegend{1} = 'Cross-Correlation';
   else
      flegend{1} = 'Auto-Correlation';
   end
   flegend{2} = 'g(r) Random';
   flegend{3} = 'Fit';
   axis tight
   if obj.Rmax_axis > 0
      xlim([0, obj.Rmax_axis]);
   end
   title(regexprep(base_name, '_', '\\_'));
   xlabel('r (nm)');
   ylabel('g(r)');
   legend(flegend, 'Location', 'Best');
   hold off

   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
      saveas(gcf, name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   Afound=estimates(1);
   s_d_found_pixels=estimates(2);
   Bfound=estimates(3);
   s_l_found_pixels=estimates(4);

   Afound_SE=errors(1);
   s_d_found_pixels_SE=errors(2);
   Bfound_SE=errors(3);
   s_l_found_pixels_SE=errors(4);

   s_d_found=s_d_found_pixels*hist_bin_size;
   s_l_found=s_l_found_pixels*hist_bin_size;
   s_d_found_SE=s_d_found_pixels_SE*hist_bin_size;
   s_l_found_SE=s_l_found_pixels_SE*hist_bin_size;

   % display results
   if corr_type == 'C'
      fprintf('Pair Cross-correlation for %s:\n\n', base_name);
   else
      fprintf('Pair Auto-correlation for %s:\n\n', base_name);
   end

   fprintf('Objects per domain ');
   fprintf('found:\t%5.3g +/- %5.3g\n',Afound,Afound_SE);

   fprintf('Domain size sigma ');
   fprintf('found:\t%5.3g +/- %5.3g\n',s_d_found,s_d_found_SE);

   fprintf('Localizations per object ');
   fprintf('found:\t%5.3g +/- %5.3g\n',Bfound,Bfound_SE);

   fprintf('Localization precision sigma ');
   fprintf('found:\t%5.3g +/- %5.3g\n',s_l_found,s_l_found_SE);

   results.G  = G;
   results.r  = r_nm;
   results.g  = g;
   results.dg = dg;
   results.objs_per_domain          = Afound;
   results.objs_per_domain_SE       = Afound_SE;
   results.sigma_domain             = s_d_found;
   results.sigma_domain_SE          = s_d_found_SE;
   results.localizations_per_obj    = Bfound;
   results.localizations_per_obj_SE = Bfound_SE;
   results.sigma_localization       = s_l_found;
   results.sigma_localization_SE    = s_l_found_SE;
   results.model = model;

end

% =============================================================================

function results = pair_correlation_Veatch(obj, XY1, XY2, hist_bin_size, ...
                                           base_name, correlation, fit)
% Pair correlation as originally written by Sarah L. Veatch.
%
% Input:
%    XY1, XY2         (x, y) coordinates of the two datasets (nm) [N x 2]
%    hist_bin_size    pixel size (nm)
%    base_name        descriptive name for the results files
%    correlation      ['auto'] 'auto' or 'cross'
%    fit              'exponential_and_gaussian', 'exponential_and_cosine' or
%                     'exponential'
% Output:
%    results         structure containing various results from the algorithm

   if ~isempty(obj.ROI)
      x_min = obj.ROI(1);
      y_min = obj.ROI(2);
      x_max = x_min + obj.ROI(3);
      y_max = y_min + obj.ROI(4);
      ROI_size = min(obj.ROI(3 : 4));
   else
      if ~isempty(XY2)
         x_min = min([XY1(:, 1); XY2(:, 1)]);
         y_min = min([XY1(:, 2); XY2(:, 2)]);
      else
         x_min = min(XY1(:, 1));
         y_min = min(XY1(:, 2));
      end
      if ~isempty(obj.ROI_size)
         x_max = x_min + obj.ROI_size;
         y_max = y_min + obj.ROI_size;
         ROI_size = obj.ROI_size;
      else
         if ~isempty(XY2)
            x_max = max([XY1(:, 1); XY2(:, 1)]);
            y_max = max([XY1(:, 2); XY2(:, 2)]);
         else
            x_max = max(XY1(:, 1));
            y_max = max(XY1(:, 2));
         end
         ROI_size = min(x_max - x_min, y_max - y_min);
      end
   end

   dataA(:, 1) = XY1(:, 1) - x_min;
   dataA(:, 2) = XY1(:, 2) - y_min;
   if ~isempty(XY2)
      dataB(:, 1) = XY2(:, 1) - x_min;
      dataB(:, 2) = XY2(:, 2) - y_min;
   else
      dataB = [];
   end

%% Filter the localisations
%rule = 1: precision or intensity
%rule = 2: both
% rule = 1;
% limits = [10 50]; %[min max]
% dataCol = 7; % variable to be used for filtering. Use a vector for two parameters eg [1 3]
% dataA = filter_localisations(dataA,dataCol,rule,limits);
% if ~isempty(dataB)
%     dataB = filter_localisations(dataB,dataCol,rule,limits);
% end
%% set some variables:
camPixSize = hist_bin_size; %pixel size on ccd in nm
originalX = round((x_max - x_min)/hist_bin_size);
originalY = round((y_max - y_min)/hist_bin_size); % image size in pixels
xScale = originalX .* camPixSize;
yScale = originalY .* camPixSize;
% nmPixSizeX = xScale / image_resolution;
% nmPixSizeY = yScale / image_resolution;
%nmPixSize = sqrt(nmPixSizeX^2 + nmPixSizeY^2); % pixel size in 2D histogram
nmPixSize = hist_bin_size;
image_resolution = [xScale/nmPixSize, yScale/nmPixSize]; % resolution for 2D histogram of localisation data
%image_resolution(1)
%% apply channel alignment?
transformation = []; %enter filename to apply transformation
calc_new = 0;
t_params = {transformation, calc_new};

%% set the type of correlation and the function to fit to the data
%%'auto' for auto-correlation, 'cross' for cross-correlation
if ~exist('correlation', 'var')
   correlation = 'auto';
   %correlation = 'cross';
end
if ~exist('fit', 'var')
   fit = obj.Veatch_fit;
   %fit = 'exponential_and_gaussian'; %name should match available fit functions
   %fit = 'exponential_and_cosine'; %name should match available fit functions
   %fit = 'exponential'; %name should match available fit functions
end

%radius = 1000; %in pixels
% Establish rmax as half the size of the ROI in pixels.
radius = round(ROI_size / (2 * hist_bin_size));

%% extract the (possibly filtered) x-y coordinates
% this step converts the coordinates pixels to nm; remove multiplication by
% camPixSize to work with data already in nm
Ax = dataA(:, 1);
Ay = dataA(:, 2);
Axy = [Ax Ay];
%Axy = [Ax Ay].*camPixSize;
if ~isempty(dataB)
    Bx = dataB(:, 1);
    By = dataB(:, 2);
    Bxy = [Bx By];
    %Bxy = [Bx By].*camPixSize;
else
    Bxy = [];
end

%% calculate the correlation and the fit
[correlation_data, params] = obj.run_correlation_and_fit(Axy, Bxy, image_resolution, nmPixSize, t_params, [xScale yScale], correlation, fit, radius, base_name);

fprintf('Veatch %s-correlation using an %s fit for %s:', ...
        correlation, fit, base_name);
%correlation_data
params{1}
if numel(params) > 1
   params{2}
end

results.correlation_data = correlation_data;
results.params = params{1};

end

% -----------------------------------------------------------------------------

function [corrData, params] = run_correlation_and_fit(obj, varargin) %(input1, input2, Xcol, Ycol, res, nm pix size, range, correlation type, fittype, radius)

    channel1 = varargin{1};
    channel2 = varargin{2};
    res = varargin{3};
    nmpixSize = varargin{4};
    t_params = varargin{5};
    range = varargin{6};
    Xrange = [0 range(1)];
    Yrange = [0 range(2)];
    correlation = varargin{7};
    fit = varargin{8};
    maxrad1 = varargin{9};
    base_name = varargin{10};

    %calculate histograms
    if ~isempty(channel2)
        density(:,:,1) = PairCorr.hist2d(channel1,res(1), res(2),Xrange,Yrange);
        density(:,:,2) = PairCorr.hist2d(channel2,res(1), res(2),Xrange,Yrange);
        density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
    else
        density = PairCorr.hist2d(channel1,res(1), res(2),Xrange,Yrange);
    end
    N = {};
    data = {};
    N{1} = size(channel1,1);
    if isempty(channel2)
        numChannels = 1;
        data{1} = channel1;
        data{2} = [];
    else
        numChannels = 2;
        N{2} = size(channel2,1);
        if ~isempty(t_params{1}) && t_params{2} == 0
            tformData = load(t_params{1});
            TFORM = tformData.TFORM;
            data{1} = tformfwd(TFORM,channel1);
        elseif isempty(t_params{1}) && t_params{2} == 0
            data{1} = channel1;
        elseif isempty(t_params{1}) && t_params{2} == 1
            [in_points,base_points,~,TFORM] = transformChannels(nmpixSize,density);
            data{1} = tformfwd(TFORM,channel1);
        end
        data{2} = channel2;
    end
        %calculate histograms
    if numChannels == 2
        density(:,:,1) = PairCorr.hist2d(data{1},res(1), res(2),Xrange,Yrange);
        density(:,:,2) = PairCorr.hist2d(data{2},res(1), res(2),Xrange,Yrange);
        density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
    else
        density = PairCorr.hist2d(data{1},res(1), res(2),Xrange,Yrange);
        cmax = 5;
        cmin = min(min(density));
    end

%   %display for ROI definition
%   LUT = RedMap;
%   r = LUT(:,2);
%   g = LUT(:,1);
%   b = LUT(:,3);
%   k = 1:numel(r);
%   map(k,:)=[r(k) g(k) b(k)];
%   figure;hIm = imagesc(Xrange,Yrange,density,[cmin,cmax]); axis equal tight off; colormap(map)
%   %mask = roipoly;
%   h = imrect;
%   pos = getPosition(h)
%   mask = createMask(h,hIm);
%   figure;imagesc(mask);axis equal tight off;colormap('gray')
%   figure;imagesc(density.*mask,[cmin,cmax]);axis equal tight off; colormap(map)
    params = {};
    corrData = repmat(struct('twoDcorr',[],'radius',[],'correlation',[],'error',[],'mask',[],'type',[],'L',[]),numChannels,1);
    L = {};
%     for ii = 1:numChannels
%         [corrData(ii).L,vq] = PairCorr.Ripley(h,data{ii},maxrad1*nmpixSize);
%         fname = sprintf('clustermap0%d.tif',ii);
%         imwrite(uint16(vq),fname,'tif','compression','lzw')
%     end
    % the input parameters for the fit depend on the fit type
    switch fit
        case 'exponential_and_gaussian'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
            P(3) = 20; %sqrt(2)* the PSF of the image. (the half with of the gaussian) in nm
            P(4) = 1; % the surface density of the probe (in 1/um^2)
        case 'exponential_and_cosine'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
            P(3) = 20; %sqrt(2)* the PSF of the image. (the half with of the gaussian) in nm
        case 'exponential'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
    end

    %calculate correlation
    % Set to 1 to get a figure of g(r) with error bars.
    flag = 0;
    switch correlation
        case 'auto'
            for iChan = 1: numChannels
                Imsize1=min(size(density(:,:,1)));
                if Imsize1<1.25*maxrad1
                    maxrad1=round(Imsize1/1.25);
                end
                mask = ones(size(density(:,:,iChan)));
                [G, r, g, dg, maskout] = PairCorr.get_autocorr(im2bw(density(:,:,iChan)), mask, maxrad1, flag);
                if isnan(g)
                    errordlg('auto-correlation calculation failed. try changing radius','modal');
                    return
                else
                    params{iChan} = obj.fitData(r,g,dg,nmpixSize,fit,P,base_name,correlation);
                    corrData(iChan).twoDcorr = G;
                    corrData(iChan).radius = r;
                    corrData(iChan).correlation = g;
                    corrData(iChan).error = dg;
                    corrData(iChan).mask = maskout;
                    corrData(iChan).type = 'auto';
                end
            end
        case 'cross'
            Imsize1=min(size(density(:,:,1)));
            if Imsize1<1.25*maxrad1
                maxrad1=round(Imsize1/1.25);
            end

            if numChannels == 2
                mask = ones(size(density(:,:,1)));
                [C, r, c, dc, maskout] = PairCorr.get_crosscorr(density(:,:,1), density(:,:,2), mask, maxrad1, flag);
            else
                errordlg('two channels are required to calculate cross-correlation','modal');
                return
            end

            if isnan(c)
                errordlg('cross-correlation calculation failed. try changing radius','modal');
                return
            else
                params{1} = obj.fitData(r,c,dc,nmpixSize,fit,P,base_name,correlation);
                corrData(1).twoDcorr = C;
                corrData(1).radius = r;
                corrData(1).correlation = c;
                corrData(1).error = dc;
                corrData(1).mask = maskout;
                corrData(1).type = 'cross';
            end
    end
end

% -----------------------------------------------------------------------------

function params = fitData(obj,x,y,err,pixSize,type,Pin,base_name,correlation)
%fit the data

    switch type
        case 'exponential_and_gaussian'
            params = repmat(struct('cluster_size',[],'magnitude',[],'density',[],'sigma',[]),1,1);
        case 'exponential_and_cosine'
            params = repmat(struct('cluster_size',[],'magnitude',[],'r0',[]),1,1);
        case 'exponential'
            params = repmat(struct('cluster_size',[],'magnitude',[]),1,1);
    end

    x = x .* pixSize;
    if ~isempty(obj.Fig_ext)
       figure('Visible', 'off');
    else
       figure;
    end
    axes(obj.Font_props{:});
    hold on
    errorbar(x(2:end), y(2:end), err(2:end), 'k.', 'LineWidth', 2)
    plot(x(2:end),ones(1, numel(x) - 1), 'b:', 'LineWidth', 3)
    switch type
        case 'exponential_and_gaussian'
            P0 = [Pin(1), Pin(2), Pin(3), Pin(4)];

            P = lsqcurvefit(@(P, r) PairCorr.exponential_and_gaussian(P, r), P0, x(2:end), y(2:end)-1);
            %P = lsqcurvefit('PairCorr.exponential_and_gaussian', P0, x(2:end), y(2:end)-1);

            %hold on
            plot(x, PairCorr.exponential_and_gaussian(P, x)+1, 'r', 'LineWidth', 3)
            %hold off

            %legend ('data', 'fit')

            params.cluster_size = P(1); %characteristic size of structure
            params.magnitude = P(2); % magnitude of clustering
            params.sigma = P(3)/sqrt(2); %in nm
            params.density = P(4); %in 1/um^2
        case 'exponential_and_cosine'
            P0 = [Pin(1), Pin(2), Pin(3)];
            P = lsqcurvefit(@(P, r) PairCorr.exponential_and_cosine(P, r), P0, x(2:end), y(2:end)-1);
            %P = lsqcurvefit('PairCorr.exponential_and_cosine', P0, x(2:end), y(2:end)-1);

            %hold on
            plot(x, PairCorr.exponential_and_cosine(P, x)+1, 'r', 'LineWidth', 3)
            %hold off

            %legend ('data', 'fit')

            params.cluster_size = P(1); %characteristic size of structure
            params.magnitude = P(2); % magnitude of clustering
            params.r0 = P(3); %in nm
        case 'exponential'
            P0 = [Pin(1), Pin(2)];
            P = lsqcurvefit(@(P, r) PairCorr.exponential(P, r), P0, x(2:end), y(2:end)-1);
            %P = lsqcurvefit('PairCorr.exponential', P0, x(2:end), y(2:end)-1);

            %hold on
            plot(x, PairCorr.exponential(P, x)+1, 'r', 'LineWidth', 3)
            %hold off

            %legend ('data', 'fit')

            params.cluster_size = P(1); %characteristic size of structure
            params.magnitude = P(2); % magnitude of clustering
    end
    legend(sprintf('%s-correlation data', correlation), 'g(r) Random', 'fit');
    axis tight
    xlabel('r (nm)');
    ylabel('g(r)');
    title(sprintf('%s (Veatch)', regexprep(base_name, '_', '\\_')));
    hold off
    name = fullfile(obj.Results, sprintf('%s_%scorrV', base_name, correlation));
    if ~isempty(obj.Fig_ext)
       print(['-d', obj.Fig_ext], name);
       saveas(gcf, name);
    else
       saveas(gcf, name);
       delete(gcf);
    end
end

% =============================================================================
end % methods

methods(Static)
% =============================================================================

function [C, r, c, dc, rmax] = get_corr(n_ROIs, rmax, II1, II2)
% function [G, r, g, dg] = get_crosscorr(I1, I2, rmax)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% I1 = First  image to be autocorrelated
% I2 = Second image to be autocorrelated
% rmax = maximum r value to correlate in units of pixels.
%
% OUTPUTS
% C = two dimensional cross-correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% c = angularly averaged autocorrelation function.
% dc = errors on angularly averaged g
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 2D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.24.11 by Sarah Veatch.

if exist('II2', 'var')
   corr_type = 'C';
else
   corr_type = 'A';
end

L1min = 1e+10;
L2min = 1e+10;
L1max = -1;
L2max = -1;
for i = 1 : n_ROIs
   L1 = size(II1{i}, 1)+rmax; % size of fft2 (for zero padding)
   L2 = size(II1{i}, 2)+rmax; % size of fft2 (for zero padding)
   L1min = min(L1min, L1);
   L2min = min(L2min, L2);
   L1max = max(L1max, L1);
   L2max = max(L2max, L2);
end
Lmax_min = min(L1max, L2max);

% Adjust rmax if it is too big and would cause problems cropping C1:
%    Need floor(Lmax_min/2+1) - rmax >= 1
if floor(Lmax_min/2) < rmax
   fprintf('rmax adjusted from %.1f to ', rmax);
   rmax = floor(Lmax_min/2);
   fprintf('%.1f\n', rmax);
end

sum_A  = 0;
sum_N1 = 0;
sum_N2 = 0;
sum_FF = 0;
sum_NP = 0;
for i = 1 : n_ROIs
   I1 = II1{i};
   if corr_type == 'C'
      I2 = II2{i};
      if sum(size(I1)==size(I2))<2,
         disp('images are not the same size')
         return
      end
   end

   N1 = sum(sum(I1));        % Total intensity in I1
   if corr_type == 'C'
      N2 = sum(sum(I2));     % Total intensity in I2
   end
   A = prod(size(I1));       % area of mask
   mask = ones(size(I1));

   I1 = double(I1);          % convert to double
   if corr_type == 'C'
      I2 = double(I2);       % convert to double
   end

   %L1 = size(I1, 1)+rmax;    % size of fft2 (for zero padding)
   %L2 = size(I1, 2)+rmax;    % size of fft2 (for zero padding)

   % Normalization for correct boundary conditions
   % Center the (L1, L2) size mask within the maximal (L1max, L2max) extents.
   NP = real(fftshift(ifft2(abs(fft2(mask, L1max, L2max)).^2)));
   %NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2)));
   % Collect the image areas (C1_A), image intensities (C1_N1 and C1_N2), FFT
   % transforms (C1_FF) and normalizations (C1_NP) separately, the last two
   % centered within the maximal rectangular extents (L1max, L2max), computing
   % the averaged C1 after the end of the loop.
   if corr_type == 'C'
      C1_A  = A;
      C1_N1 = N1;
      C1_N2 = N2;
      C1_FF = real(fftshift(ifft2(fft2(I1, L1max, L2max).* ...
                             conj(fft2(I2, L1max, L2max)))));
      C1_NP = NP;
      %C1 = A^2/N1/N2*real(fftshift(ifft2(fft2(I1, L1, L2).* ...
      %                              conj(fft2(I2, L1, L2)))))./NP;
   else
      % 2D G(r) with proper normalization
      C1_A  = A;
      C1_N1 = N1;
      C1_N2 = N1;
      C1_FF = real(fftshift(ifft2(abs(fft2(I1, L1max, L2max)).^2)));
      C1_NP = NP;
      %C1 = A^2/N1^2*real(fftshift(ifft2(abs(fft2(I1, L1, L2)).^2)))./NP;
   end

   sum_A  = sum_A  + C1_A ;
   sum_N1 = sum_N1 + C1_N1;
   sum_N2 = sum_N2 + C1_N2;
   sum_FF = sum_FF + C1_FF;
   sum_NP = sum_NP + C1_NP;
end
C1 = sum_A^2 / (sum_N1 * sum_N2) * sum_FF ./ sum_NP;

%only return valid part of G:
% a square with sides 2*rmax+1 about the center pixel, that is, the minimal
% rectangular region corresponding to the overlap of all the ROIs.
C = imcrop(C1, [floor(L2max/2+1)-rmax, floor(L1max/2+1)-rmax, 2*rmax, 2*rmax]);
%C = imcrop(C1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);

xvals = ones(1, 2*rmax+1)'*(-rmax:rmax); %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1); %map to y positions with center y=0
zvals = C;

% convert x, y to polar coordinates
[theta,r,v] = cart2pol(xvals, yvals, zvals);

% Label each pixel in the minimal square region by its radius r from the
% center, then sort the r's and bin them.  Collect the pixels with radii in
% each bin and average the results to compute c(r).
Ar = reshape(r,1, (2*rmax+1)^2);
Avals = reshape(v,1, (2*rmax+1)^2);
[rr,ind] = sort(Ar);                  % sort by r values
vv = Avals(ind);                      % reindex g
r = 0:floor(max(rr));                 % the radii you want to extract
[n bin] = histc(rr, r-.5);            % bin by radius

for j = 1:rmax+1;                     % now get averages
   m = bin==j;
   n2 = sum(m);                       % the number of pixels in that bin
   if n2==0, vals(j)=0; er(j)=0;      % if no bins, no data
   else
      c(j) = sum(m.*vv)/n2;           % the average G values in this bin
                                      % the variance of the mean
      dc(j) = sqrt(sum(m.*(vv-c(j)).^2))/n2;
   end
end

r = 0:rmax;

end

% -----------------------------------------------------------------------------

function M2 = ext(M1, m_max, n_max)
% Place m x n M1 into the center of m_max x n_max M2, where all dimensions are
% assumed to be odd positive integers.

   [m, n] = size(M1);
   m_ext = floor(m/2);
   n_ext = floor(n/2);

   m_max_c = ceil(m_max / 2);   % center pixel
   n_max_c = ceil(n_max / 2);   % center pixel

   M2 = zeros(m_max, n_max);
   M2(m_max_c - m_ext : m_max_c + m_ext, ...
      n_max_c - n_ext : n_max_c + n_ext) = M1;

end

% -----------------------------------------------------------------------------

function [C, r, c, dc, mask, rmax] = get_crosscorr(I1, I2, mask, rmax, flag)
% function [G, r, g, dg, mask] = get_crosscorr(I1, I2, mask, rmax, flag)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% I1 = First image to be autocorrelated
% I2 = Second image to be autocorrelated
% mask = region of interest, If none or [] specificed, user will be asked
%   to define.  To autocorreate the entire images, use mask = ones(size(I1))
% rmax = maximum r value to correlate in units of pixels. default is 100;
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation.
%
% OUTPUTS
% C = two dimensional cross-correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% c = angularly averaged autocorrelation function.
% dc = errors on angularly averaged g
% mask = masked used for calculation
%
%
% Last updated 01.24.11 by Sarah Veatch.



if nargin<5, flag = 0; end  % flag for display
if (nargin<4 || isempty(rmax)), rmax=100; end  % distance of maximum correlation recorded
if (nargin<3 || isempty(mask)),    %% draw a mask if needed
    hIm = imshow(I1); axis equal tight off;
    %mask = roipoly;
    h = imrect;
    mask = createMask(h,hIm);
end

if sum(size(I1)==size(I2))<2,
    disp('images are not the same size')
    return
end

%mask = ones(size(I1)); end

N1 = sum(sum(I1.*mask));  % Total intensity within mask in I1
N2 = sum(sum(I2.*mask));  % Total intensity within mask in I2
A = sum(sum(mask));      % area of mask

I1 = double(I1);         % convert to double
I2 = double(I2);         % convert to double

L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)

% Adjust rmax if it is too big and would cause problems cropping C1:
%    Need 2*rmax + 1 >= min(size(I1, 1), size(I1, 2)) + rmax
% This fix is for rectangular regions which are treated somewhat differently in
% get_corr, so rectangular regions may produce different results in the two
% routines---use with caution!  (Square region results are identical.)
Lmin = min(L1, L2);
if 2*rmax + 1 > Lmin
   fprintf('rmax adjusted from %.1f to ', rmax);
   rmax = min(size(I1, 1), size(I1, 2)) - 1;
   fprintf('%.1f\n', rmax);

   L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
   L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)
end

NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2))); % Normalization for correct boundary conditions
C1 = A^2/N1/N2*real(fftshift(ifft2(fft2(I1.*mask,L1, L2).*conj(fft2(I2.*mask, L1, L2)))))./NP;
%G1 = A^2/N^2*real(fftshift(ifft2(abs(fft2(I1.*mask,L1, L2)).^2)))./NP; % 2D G(r) with proper normalization
C = imcrop(C1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);  %only return valid part of G


xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    %map to y positions with center y=0
zvals = C;

[theta,r,v] = cart2pol(xvals,yvals, zvals);  % convert x, y to polar coordinates

Ar = reshape(r,1, (2*rmax+1)^2);
Avals = reshape(v,1, (2*rmax+1)^2);
[rr,ind] = sort(Ar);                         % sort by r values
vv = Avals(ind);                             % reindex g
r = 0:floor(max(rr));                        % the radii you want to extract
[n bin] = histc(rr, r-.5);                   % bin by radius

for j = 1:rmax+1;                            % now get averages
    m = bin==j;
    n2 = sum(m);                             % the number of pixels in that bin
    if n2==0, vals(j)=0; er(j)=0;            % if no bins, no data
    else
        c(j) = sum(m.*vv)/n2;               % the average G values in this bin
        dc(j) = sqrt(sum(m.*(vv-c(j)).^2))/n2; % the variance of the mean
    end
end

r = 0:rmax;

%end


if flag,
    r = 0:rmax;
    figure;errorbar(r, c, dc);
    axis tight
end

end

% -----------------------------------------------------------------------------

function [G, r, g, dg, mask, rmax] = get_autocorr(I1 , mask, rmax, flag)
% function [G, r, g, dg, mask] = get_autocorr(I1 , mask, rmax, flag)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% I1 = image to be autocorrelated
% mask = region of interest, If none or [] specificed, user will be asked
%   to define.  To autocorreate the entire images, use mask = ones(size(I1))
% rmax = maximum r value to correlate in units of pixels. default is 100;
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation.
%
% OUTPUTS
% G = two dimensional correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% g = angularly averaged autocorrelation function.
% dg = errors on angularly averaged g
% mask = masked used for calculation
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 2D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.26.10 by Sarah Veatch.



if nargin<4, flag = 0; end  % flag for display
if (nargin<3 || isempty(rmax)), rmax=100; end  % distance of maximum correlation recorded
if (nargin<2 || isempty(mask)),    %% draw a mask if needed
    imagesc(I1); axis equal tight off;
    mask = roipoly;
end

%mask = ones(size(I1)); end

N = sum(sum(I1.*mask));  % number of particles within mask
A = sum(sum(mask));      % area of mask

I1 = double(I1);         % convert to double

L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)

% Adjust rmax if it is too big and would cause problems cropping C1:
%    Need 2*rmax + 1 >= min(size(I1, 1), size(I1, 2)) + rmax
% This fix is for rectangular regions which are treated somewhat differently in
% get_corr, so rectangular regions may produce different results in the two
% routines---use with caution!  (Square region results are identical.)
Lmin = min(L1, L2);
if 2*rmax + 1 > Lmin
   fprintf('rmax adjusted from %.1f to ', rmax);
   rmax = min(size(I1, 1), size(I1, 2)) - 1;
   fprintf('%.1f\n', rmax);

   L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
   L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)
end

NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2))); % Normalization for correct boundary conditions
G1 = A^2/N^2*real(fftshift(ifft2(abs(fft2(I1.*mask,L1, L2)).^2)))./NP; % 2D G(r) with proper normalization
G = imcrop(G1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);  %only return valid part of G


xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    %map to y positions with center y=0
zvals = G;

[theta,r,v] = cart2pol(xvals,yvals, zvals);  % convert x, y to polar coordinates

Ar = reshape(r,1, (2*rmax+1)^2);
Avals = reshape(v,1, (2*rmax+1)^2);
[rr,ind] = sort(Ar);                         % sort by r values
vv = Avals(ind);                             % reindex g
r = 0:floor(max(rr));                        % the radii you want to extract
[n bin] = histc(rr, r-.5);                   % bin by radius

for j = 1:rmax+1;                            % now get averages
    m = bin==j;
    n2 = sum(m);                             % the number of pixels in that bin
    if n2==0, vals(j)=0; er(j)=0;            % if no bins, no data
    else
        g(j) = sum(m.*vv)/n2;               % the average G values in this bin
        dg(j) = sqrt(sum(m.*(vv-g(j)).^2))/n2; % the variance of the mean
    end
end

r = 0:rmax;

%end

G(rmax+1, rmax+1) = 0;

if flag,
    r = 0:rmax;
    errorbar(r(2:length(r)), g(2:length(r)), dg(2:length(r)));
    axis tight
end

end

% =============================================================================

% Fitting pair correlation results written by Keith Lidke
    function [ estimates errors model ] = PC_GaussFit( r,g_r,rmax,rho )
    %PC_GaussFit Fit Pair Correlation to Gaussian domain model
    %
    % INPUTS
    %   r:      correlation length
    %   g_r:    radially averaged correlation
    %   rmax:   limit of correlation fit region
    %   rho:    localization density
    %
    % OUTPUTS
    %   estimates:  estimated parameter values
    %       estimates(1):   objects per domain
    %       estimates(2):   sigma for 2D gaussian domain size
    %       estimages(3):   observations per object
    %       estimates(4):   sigma for Gaussian localization precision
    %   errors:     standard error on estimates calculated as
    %       sqrt(diag(inv(hessian)))
    %   model:      model caculated at estimated value
    %
    %   NOTES:
    %
    %   Requires DERIVESTsuite for hessian calculation
    %
    %   r(1) and g_r(1) not used in fit
    A = sqrt(g_r(2));
    B = A;
    X0=abs([A 20 B 20]);
    fitfunc=@PairCorr.Gauss2D;

    rfit=r(2:rmax);
    g_rfit=g_r(2:rmax);

    opts=optimset('MaxFunEvals',1e4,'Display','off','TolX',10^-6,'TolFun',10^-6);
    [estimates sse]= fminsearch(fitfunc,X0,opts,rfit,g_rfit,rho);

    anonfunc=@(X) PairCorr.Gauss2D(X,rfit,g_rfit,rho);
    [H Herrs]=PairCorr.hessian(anonfunc,estimates);
    errors=sqrt(diag(inv(H)))';

    [out model]=PairCorr.Gauss2D(estimates,r,g_r,rho);

    end

    function [out model]=Gauss2D(X0,r,g_r,rho)

    gnorm=inline('1/(2*pi*s^2)*exp(-r.^2/(2*s^2))','s','r');

    A=X0(1);
    sigma_dom=X0(2);
    B=X0(3);
    sigma_loc=X0(4);

    gr_psf=1/rho*gnorm(sqrt(2)*sigma_loc,r)*B;
    seff=sqrt(2*sigma_loc^2+2*sigma_dom^2);
    gr_dom=1/rho*gnorm(seff,r)*A*B;
    model=gr_psf+gr_dom+1;
    out=mean((model-g_r).^2); %mse(model,g_r);

    end

% =============================================================================

function [der,errest,finaldelta] = derivest(fun,x0,varargin)
% DERIVEST: estimate the n'th derivative of fun at x0, provide an error estimate
% usage: [der,errest] = DERIVEST(fun,x0)  % first derivative
% usage: [der,errest] = DERIVEST(fun,x0,prop1,val1,prop2,val2,...)
%
% Derivest will perform numerical differentiation of an
% analytical function provided in fun. It will not
% differentiate a function provided as data. Use gradient
% for that purpose, or differentiate a spline model.
%
% The methods used by DERIVEST are finite difference
% approximations of various orders, coupled with a generalized
% (multiple term) Romberg extrapolation. This also yields
% the error estimate provided. DERIVEST uses a semi-adaptive
% scheme to provide the best estimate that it can by its
% automatic choice of a differencing interval.
%
% Finally, While I have not written this function for the
% absolute maximum speed, speed was a major consideration
% in the algorithmic design. Maximum accuracy was my main goal.
%
%
% Arguments (input)
%  fun - function to differentiate. May be an inline function,
%        anonymous, or an m-file. fun will be sampled at a set
%        of distinct points for each element of x0. If there are
%        additional parameters to be passed into fun, then use of
%        an anonymous function is recommended.
%
%        fun should be vectorized to allow evaluation at multiple
%        locations at once. This will provide the best possible
%        speed. IF fun is not so vectorized, then you MUST set
%        'vectorized' property to 'no', so that derivest will
%        then call your function sequentially instead.
%
%        Fun is assumed to return a result of the same
%        shape as its input x0.
%
%  x0  - scalar, vector, or array of points at which to
%        differentiate fun.
%
% Additional inputs must be in the form of property/value pairs.
%  Properties are character strings. They may be shortened
%  to the extent that they are unambiguous. Properties are
%  not case sensitive. Valid property names are:
%
%  'DerivativeOrder', 'MethodOrder', 'Style', 'RombergTerms'
%  'FixedStep', 'MaxStep'
%
%  All properties have default values, chosen as intelligently
%  as I could manage. Values that are character strings may
%  also be unambiguously shortened. The legal values for each
%  property are:
%
%  'DerivativeOrder' - specifies the derivative order estimated.
%        Must be a positive integer from the set [1,2,3,4].
%
%        DEFAULT: 1 (first derivative of fun)
%
%  'MethodOrder' - specifies the order of the basic method
%        used for the estimation.
%
%        For 'central' methods, must be a positive integer
%        from the set [2,4].
%
%        For 'forward' or 'backward' difference methods,
%        must be a positive integer from the set [1,2,3,4].
%
%        DEFAULT: 4 (a second order method)
%
%        Note: higher order methods will generally be more
%        accurate, but may also suffere more from numerical
%        problems.
%
%        Note: First order methods would usually not be
%        recommended.
%
%  'Style' - specifies the style of the basic method
%        used for the estimation. 'central', 'forward',
%        or 'backwards' difference methods are used.
%
%        Must be one of 'Central', 'forward', 'backward'.
%
%        DEFAULT: 'Central'
%
%        Note: Central difference methods are usually the
%        most accurate, but sometiems one must not allow
%        evaluation in one direction or the other.
%
%  'RombergTerms' - Allows the user to specify the generalized
%        Romberg extrapolation method used, or turn it off
%        completely.
%
%        Must be a positive integer from the set [0,1,2,3].
%
%        DEFAULT: 2 (Two Romberg terms)
%
%        Note: 0 disables the Romberg step completely.
%
%  'FixedStep' - Allows the specification of a fixed step
%        size, preventing the adaptive logic from working.
%        This will be considerably faster, but not necessarily
%        as accurate as allowing the adaptive logic to run.
%
%        DEFAULT: []
%
%        Note: If specified, 'FixedStep' will define the
%        maximum excursion from x0 that will be used.
%
%  'Vectorized' - Derivest will normally assume that your
%        function can be safely evaluated at multiple locations
%        in a single call. This would minimize the overhead of
%        a loop and additional function call overhead. Some
%        functions are not easily vectorizable, but you may
%        (if your matlab release is new enough) be able to use
%        arrayfun to accomplish the vectorization.
%
%        When all else fails, set the 'vectorized' property
%        to 'no'. This will cause derivest to loop over the
%        successive function calls.
%
%        DEFAULT: 'yes'
%
%
%  'MaxStep' - Specifies the maximum excursion from x0 that
%        will be allowed, as a multiple of x0.
%
%        DEFAULT: 100
%
%  'StepRatio' - Derivest uses a proportionally cascaded
%        series of function evaluations, moving away from your
%        point of evaluation. The StepRatio is the ratio used
%        between sequential steps.
%
%        DEFAULT: 2.0000001
%
%        Note: use of a non-integer stepratio is intentional,
%        to avoid integer multiples of the period of a periodic
%        function under some circumstances.
%
%
% See the document DERIVEST.pdf for more explanation of the
% algorithms behind the parameters of DERIVEST. In most cases,
% I have chosen good values for these parameters, so the user
% should never need to specify anything other than possibly
% the DerivativeOrder. I've also tried to make my code robust
% enough that it will not need much. But complete flexibility
% is in there for your use.
%
%
% Arguments: (output)
%  der - derivative estimate for each element of x0
%        der will have the same shape as x0.
%
%  errest - 95% uncertainty estimate of the derivative, such that
%
%        abs(der(j) - f'(x0(j))) < erest(j)
%
%  finaldelta - The final overall stepsize chosen by DERIVEST
%
%
% Example usage:
%  First derivative of exp(x), at x == 1
%   [d,e]=derivest(@(x) exp(x),1)
%   d =
%       2.71828182845904
%
%   e =
%       1.02015503167879e-14
%
%  True derivative
%   exp(1)
%   ans =
%       2.71828182845905
%
% Example usage:
%  Third derivative of x.^3+x.^4, at x = [0,1]
%   derivest(@(x) x.^3 + x.^4,[0 1],'deriv',3)
%   ans =
%       6       30
%
%  True derivatives: [6,30]
%
%
% See also: gradient
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 12/27/2006

par.DerivativeOrder = 1;
par.MethodOrder = 4;
par.Style = 'central';
par.RombergTerms = 2;
par.FixedStep = [];
par.MaxStep = 100;
% setting a default stepratio as a non-integer prevents
% integer multiples of the initial point from being used.
% In turn that avoids some problems for periodic functions.
par.StepRatio = 2.0000001;
par.NominalStep = [];
par.Vectorized = 'yes';

na = length(varargin);
if (rem(na,2)==1)
  error 'Property/value pairs must come as PAIRS of arguments.'
elseif na>0
  par = PairCorr.parse_pv_pairs(par,varargin);
end
par = PairCorr.check_params(par);

% Was fun a string, or an inline/anonymous function?
if (nargin<1)
  help derivest
  return
elseif isempty(fun)
  error 'fun was not supplied.'
elseif ischar(fun)
  % a character function name
  fun = str2func(fun);
end

% no default for x0
if (nargin<2) || isempty(x0)
  error 'x0 was not supplied'
end
par.NominalStep = max(x0,0.02);

% was a single point supplied?
nx0 = size(x0);
n = prod(nx0);

% Set the steps to use.
if isempty(par.FixedStep)
  % Basic sequence of steps, relative to a stepsize of 1.
  delta = par.MaxStep*par.StepRatio .^(0:-1:-25)';
  ndel = length(delta);
else
  % Fixed, user supplied absolute sequence of steps.
  ndel = 3 + ceil(par.DerivativeOrder/2) + ...
     par.MethodOrder + par.RombergTerms;
  if par.Style(1) == 'c'
    ndel = ndel - 2;
  end
  delta = par.FixedStep*par.StepRatio .^(-(0:(ndel-1)))';
end

% generate finite differencing rule in advance.
% The rule is for a nominal unit step size, and will
% be scaled later to reflect the local step size.
fdarule = 1;
switch par.Style
  case 'central'
    % for central rules, we will reduce the load by an
    % even or odd transformation as appropriate.
    if par.MethodOrder==2
      switch par.DerivativeOrder
        case 1
          % the odd transformation did all the work
          fdarule = 1;
        case 2
          % the even transformation did all the work
          fdarule = 2;
        case 3
          % the odd transformation did most of the work, but
          % we need to kill off the linear term
          fdarule = [0 1]/PairCorr.fdamat(par.StepRatio,1,2);
        case 4
          % the even transformation did most of the work, but
          % we need to kill off the quadratic term
          fdarule = [0 1]/PairCorr.fdamat(par.StepRatio,2,2);
      end
    else
      % a 4th order method. We've already ruled out the 1st
      % order methods since these are central rules.
      switch par.DerivativeOrder
        case 1
          % the odd transformation did most of the work, but
          % we need to kill off the cubic term
          fdarule = [1 0]/PairCorr.fdamat(par.StepRatio,1,2);
        case 2
          % the even transformation did most of the work, but
          % we need to kill off the quartic term
          fdarule = [1 0]/PairCorr.fdamat(par.StepRatio,2,2);
        case 3
          % the odd transformation did much of the work, but
          % we need to kill off the linear & quintic terms
          fdarule = [0 1 0]/PairCorr.fdamat(par.StepRatio,1,3);
        case 4
          % the even transformation did much of the work, but
          % we need to kill off the quadratic and 6th order terms
          fdarule = [0 1 0]/PairCorr.fdamat(par.StepRatio,2,3);
      end
    end
  case {'forward' 'backward'}
    % These two cases are identical, except at the very end,
    % where a sign will be introduced.

    % No odd/even trans, but we already dropped
    % off the constant term
    if par.MethodOrder==1
      if par.DerivativeOrder==1
        % an easy one
        fdarule = 1;
      else
        % 2:4
        v = zeros(1,par.DerivativeOrder);
        v(par.DerivativeOrder) = 1;
        fdarule = v/PairCorr.fdamat(par.StepRatio,0,par.DerivativeOrder);
      end
    else
      % par.MethodOrder methods drop off the lower order terms,
      % plus terms directly above DerivativeOrder
      v = zeros(1,par.DerivativeOrder + par.MethodOrder - 1);
      v(par.DerivativeOrder) = 1;
      fdarule = v/PairCorr.fdamat(par.StepRatio,0,par.DerivativeOrder+par.MethodOrder-1);
    end

    % correct sign for the 'backward' rule
    if par.Style(1) == 'b'
      fdarule = -fdarule;
    end

end % switch on par.style (generating fdarule)
nfda = length(fdarule);

% will we need fun(x0)?
if (rem(par.DerivativeOrder,2) == 0) || ~strncmpi(par.Style,'central',7)
  if strcmpi(par.Vectorized,'yes')
    f_x0 = fun(x0);
  else
    % not vectorized, so loop
    f_x0 = zeros(size(x0));
    for j = 1:numel(x0)
      f_x0(j) = fun(x0(j));
    end
  end
else
  f_x0 = [];
end

% Loop over the elements of x0, reducing it to
% a scalar problem. Sorry, vectorization is not
% complete here, but this IS only a single loop.
der = zeros(nx0);
errest = der;
finaldelta = der;
for i = 1:n
  x0i = x0(i);
  h = par.NominalStep(i);

  % a central, forward or backwards differencing rule?
  % f_del is the set of all the function evaluations we
  % will generate. For a central rule, it will have the
  % even or odd transformation built in.
  if par.Style(1) == 'c'
    % A central rule, so we will need to evaluate
    % symmetrically around x0i.
    if strcmpi(par.Vectorized,'yes')
      f_plusdel = fun(x0i+h*delta);
      f_minusdel = fun(x0i-h*delta);
    else
      % not vectorized, so loop
      f_minusdel = zeros(size(delta));
      f_plusdel = zeros(size(delta));
      for j = 1:numel(delta)
        f_plusdel(j) = fun(x0i+h*delta(j));
        f_minusdel(j) = fun(x0i-h*delta(j));
      end
    end

    if ismember(par.DerivativeOrder,[1 3])
      % odd transformation
      f_del = (f_plusdel - f_minusdel)/2;
    else
      f_del = (f_plusdel + f_minusdel)/2 - f_x0(i);
    end
  elseif par.Style(1) == 'f'
    % forward rule
    % drop off the constant only
    if strcmpi(par.Vectorized,'yes')
      f_del = fun(x0i+h*delta) - f_x0(i);
    else
      % not vectorized, so loop
      f_del = zeros(size(delta));
      for j = 1:numel(delta)
        f_del(j) = fun(x0i+h*delta(j)) - f_x0(i);
      end
    end
  else
    % backward rule
    % drop off the constant only
    if strcmpi(par.Vectorized,'yes')
      f_del = fun(x0i-h*delta) - f_x0(i);
    else
      % not vectorized, so loop
      f_del = zeros(size(delta));
      for j = 1:numel(delta)
        f_del(j) = fun(x0i-h*delta(j)) - f_x0(i);
      end
    end
  end

  % check the size of f_del to ensure it was properly vectorized.
  f_del = f_del(:);
  if length(f_del)~=ndel
    error 'fun did not return the correct size result (fun must be vectorized)'
  end

  % Apply the finite difference rule at each delta, scaling
  % as appropriate for delta and the requested DerivativeOrder.
  % First, decide how many of these estimates we will end up with.
  ne = ndel + 1 - nfda - par.RombergTerms;

  % Form the initial derivative estimates from the chosen
  % finite difference method.
  der_init = PairCorr.vec2mat(f_del,ne,nfda)*fdarule.';

  % scale to reflect the local delta
  der_init = der_init(:)./(h*delta(1:ne)).^par.DerivativeOrder;

  % Each approximation that results is an approximation
  % of order par.DerivativeOrder to the desired derivative.
  % Additional (higher order, even or odd) terms in the
  % Taylor series also remain. Use a generalized (multi-term)
  % Romberg extrapolation to improve these estimates.
  switch par.Style
    case 'central'
      rombexpon = 2*(1:par.RombergTerms) + par.MethodOrder - 2;
    otherwise
      rombexpon = (1:par.RombergTerms) + par.MethodOrder - 1;
  end
  [der_romb,errors] = PairCorr.rombextrap(par.StepRatio,der_init,rombexpon);

  % Choose which result to return

  % first, trim off the
  if isempty(par.FixedStep)
    % trim off the estimates at each end of the scale
    nest = length(der_romb);
    switch par.DerivativeOrder
      case {1 2}
        trim = [1 2 nest-1 nest];
      case 3
        trim = [1:4 nest+(-3:0)];
      case 4
        trim = [1:6 nest+(-5:0)];
    end

    [der_romb,tags] = sort(der_romb);

    der_romb(trim) = [];
    tags(trim) = [];
    errors = errors(tags);
    trimdelta = delta(tags);

    [errest(i),ind] = min(errors);

    finaldelta(i) = h*trimdelta(ind);
    der(i) = der_romb(ind);
  else
    [errest(i),ind] = min(errors);
    finaldelta(i) = h*delta(ind);
    der(i) = der_romb(ind);
  end
end

end % mainline end

% ============================================
% subfunction - romberg extrapolation
% ============================================
function [der_romb,errest] = rombextrap(StepRatio,der_init,rombexpon)
% do romberg extrapolation for each estimate
%
%  StepRatio - Ratio decrease in step
%  der_init - initial derivative estimates
%  rombexpon - higher order terms to cancel using the romberg step
%
%  der_romb - derivative estimates returned
%  errest - error estimates
%  amp - noise amplification factor due to the romberg step

srinv = 1/StepRatio;

% do nothing if no romberg terms
nexpon = length(rombexpon);
rmat = ones(nexpon+2,nexpon+1);
switch nexpon
  case 0
    % rmat is simple: ones(2,1)
  case 1
    % only one romberg term
    rmat(2,2) = srinv^rombexpon;
    rmat(3,2) = srinv^(2*rombexpon);
  case 2
    % two romberg terms
    rmat(2,2:3) = srinv.^rombexpon;
    rmat(3,2:3) = srinv.^(2*rombexpon);
    rmat(4,2:3) = srinv.^(3*rombexpon);
  case 3
    % three romberg terms
    rmat(2,2:4) = srinv.^rombexpon;
    rmat(3,2:4) = srinv.^(2*rombexpon);
    rmat(4,2:4) = srinv.^(3*rombexpon);
    rmat(5,2:4) = srinv.^(4*rombexpon);
end

% qr factorization used for the extrapolation as well
% as the uncertainty estimates
[qromb,rromb] = qr(rmat,0);

% the noise amplification is further amplified by the Romberg step.
% amp = cond(rromb);

% this does the extrapolation to a zero step size.
ne = length(der_init);
rhs = PairCorr.vec2mat(der_init,nexpon+2,max(1,ne - (nexpon+2)));
rombcoefs = rromb\(qromb.'*rhs);
der_romb = rombcoefs(1,:).';

% uncertainty estimate of derivative prediction
s = sqrt(sum((rhs - rmat*rombcoefs).^2,1));
rinv = rromb\eye(nexpon+1);
cov1 = sum(rinv.^2,2); % 1 spare dof
errest = s.'*12.7062047361747*sqrt(cov1(1));

end % rombextrap


% ============================================
% subfunction - vec2mat
% ============================================
function mat = vec2mat(vec,n,m)
% forms the matrix M, such that M(i,j) = vec(i+j-1)
[i,j] = ndgrid(1:n,0:m-1);
ind = i+j;
mat = vec(ind);
if n==1
  mat = mat.';
end

end % vec2mat


% ============================================
% subfunction - fdamat
% ============================================
function mat = fdamat(sr,parity,nterms)
% Compute matrix for fda derivation.
% parity can be
%   0 (one sided, all terms included but zeroth order)
%   1 (only odd terms included)
%   2 (only even terms included)
% nterms - number of terms

% sr is the ratio between successive steps
srinv = 1./sr;

switch parity
  case 0
    % single sided rule
    [i,j] = ndgrid(1:nterms);
    c = 1./factorial(1:nterms);
    mat = c(j).*srinv.^((i-1).*j);
  case 1
    % odd order derivative
    [i,j] = ndgrid(1:nterms);
    c = 1./factorial(1:2:(2*nterms));
    mat = c(j).*srinv.^((i-1).*(2*j-1));
  case 2
    % even order derivative
    [i,j] = ndgrid(1:nterms);
    c = 1./factorial(2:2:(2*nterms));
    mat = c(j).*srinv.^((i-1).*(2*j));
end

end % fdamat



% ============================================
% subfunction - check_params
% ============================================
function par = check_params(par)
% check the parameters for acceptability
%
% Defaults
% par.DerivativeOrder = 1;
% par.MethodOrder = 2;
% par.Style = 'central';
% par.RombergTerms = 2;
% par.FixedStep = [];

% DerivativeOrder == 1 by default
if isempty(par.DerivativeOrder)
  par.DerivativeOrder = 1;
else
  if (length(par.DerivativeOrder)>1) || ~ismember(par.DerivativeOrder,1:4)
    error 'DerivativeOrder must be scalar, one of [1 2 3 4].'
  end
end

% MethodOrder == 2 by default
if isempty(par.MethodOrder)
  par.MethodOrder = 2;
else
  if (length(par.MethodOrder)>1) || ~ismember(par.MethodOrder,[1 2 3 4])
    error 'MethodOrder must be scalar, one of [1 2 3 4].'
  elseif ismember(par.MethodOrder,[1 3]) && (par.Style(1)=='c')
    error 'MethodOrder==1 or 3 is not possible with central difference methods'
  end
end

% style is char
valid = {'central', 'forward', 'backward'};
if isempty(par.Style)
  par.Style = 'central';
elseif ~ischar(par.Style)
  error 'Invalid Style: Must be character'
end
ind = find(strncmpi(par.Style,valid,length(par.Style)));
if (length(ind)==1)
  par.Style = valid{ind};
else
  error(['Invalid Style: ',par.Style])
end

% vectorized is char
valid = {'yes', 'no'};
if isempty(par.Vectorized)
  par.Vectorized = 'yes';
elseif ~ischar(par.Vectorized)
  error 'Invalid Vectorized: Must be character'
end
ind = find(strncmpi(par.Vectorized,valid,length(par.Vectorized)));
if (length(ind)==1)
  par.Vectorized = valid{ind};
else
  error(['Invalid Vectorized: ',par.Vectorized])
end

% RombergTerms == 2 by default
if isempty(par.RombergTerms)
  par.RombergTerms = 2;
else
  if (length(par.RombergTerms)>1) || ~ismember(par.RombergTerms,0:3)
    error 'RombergTerms must be scalar, one of [0 1 2 3].'
  end
end

% FixedStep == [] by default
if (length(par.FixedStep)>1) || (~isempty(par.FixedStep) && (par.FixedStep<=0))
  error 'FixedStep must be empty or a scalar, >0.'
end

% MaxStep == 10 by default
if isempty(par.MaxStep)
  par.MaxStep = 10;
elseif (length(par.MaxStep)>1) || (par.MaxStep<=0)
  error 'MaxStep must be empty or a scalar, >0.'
end

end % check_params


% ============================================
% Included subfunction - parse_pv_pairs
% ============================================
function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params =
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end

if ~isstruct(params)
  error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  p_i = lower(pv_pairs{2*i-1});
  v_i = pv_pairs{2*i};

  ind = strmatch(p_i,lpropnames,'exact');
  if isempty(ind)
    ind = find(strncmp(p_i,lpropnames,length(p_i)));
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  p_i = propnames{ind};

  % override the corresponding default in params
  params = setfield(params,p_i,v_i); %#ok

end

end % parse_pv_pairs

% -----------------------------------------------------------------------------

function [grad,err,finaldelta] = gradest(fun,x0)
% gradest: estimate of the gradient vector of an analytical function of n variables
% usage: [grad,err,finaldelta] = gradest(fun,x0)
%
% Uses derivest to provide both derivative estimates
% and error estimates. fun needs not be vectorized.
%
% arguments: (input)
%  fun - analytical function to differentiate. fun must
%        be a function of the vector or array x0.
%
%  x0  - vector location at which to differentiate fun
%        If x0 is an nxm array, then fun is assumed to be
%        a function of n*m variables.
%
% arguments: (output)
%  grad - vector of first partial derivatives of fun.
%        grad will be a row vector of length numel(x0).
%
%  err - vector of error estimates corresponding to
%        each partial derivative in grad.
%
%  finaldelta - vector of final step sizes chosen for
%        each partial derivative.
%
%
% Example:
%  [grad,err] = gradest(@(x) sum(x.^2),[1 2 3])
%  grad =
%      2     4     6
%  err =
%      5.8899e-15    1.178e-14            0
%
%
% Example:
%  At [x,y] = [1,1], compute the numerical gradient
%  of the function sin(x-y) + y*exp(x)
%
%  z = @(xy) sin(diff(xy)) + xy(2)*exp(xy(1))
%
%  [grad,err ] = gradest(z,[1 1])
%  grad =
%       1.7183       3.7183
%  err =
%    7.537e-14   1.1846e-13
%
%
% Example:
%  At the global minimizer (1,1) of the Rosenbrock function,
%  compute the gradient. It should be essentially zero.
%
%  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%  [g,err] = gradest(rosen,[1 1])
%  g =
%    1.0843e-20            0
%  err =
%    1.9075e-18            0
%
%
% See also: derivest, gradient
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/9/2007

% get the size of x0 so we can reshape
% later.
sx = size(x0);

% total number of derivatives we will need to take
nx = numel(x0);

grad = zeros(1,nx);
err = grad;
finaldelta = grad;
for ind = 1:nx
  [grad(ind),err(ind),finaldelta(ind)] = PairCorr.derivest( ...
    @(xi) fun(PairCorr.swapelement(x0,ind,xi)), ...
    x0(ind),'deriv',1,'vectorized','no', ...
    'methodorder',2);
end

end % mainline function end

% -----------------------------------------------------------------------------

function [HD,err,finaldelta] = hessdiag(fun,x0)
% HESSDIAG: diagonal elements of the Hessian matrix (vector of second partials)
% usage: [HD,err,finaldelta] = hessdiag(fun,x0)
%
% When all that you want are the diagonal elements of the hessian
% matrix, it will be more efficient to call HESSDIAG than HESSIAN.
% HESSDIAG uses DERIVEST to provide both second derivative estimates
% and error estimates. fun needs not be vectorized.
%
% arguments: (input)
%  fun - SCALAR analytical function to differentiate.
%        fun must be a function of the vector or array x0.
%
%  x0  - vector location at which to differentiate fun
%        If x0 is an nxm array, then fun is assumed to be
%        a function of n*m variables.
%
% arguments: (output)
%  HD  - vector of second partial derivatives of fun.
%        These are the diagonal elements of the Hessian
%        matrix, evaluated at x0.
%        HD will be a row vector of length numel(x0).
%
%  err - vector of error estimates corresponding to
%        each second partial derivative in HD.
%
%  finaldelta - vector of final step sizes chosen for
%        each second partial derivative.
%
%
% Example usage:
%  [HD,err] = hessdiag(@(x) x(1) + x(2)^2 + x(3)^3,[1 2 3])
%  HD =
%     0     2    18
%
%  err =
%     0     0     0
%
%
% See also: derivest, gradient, gradest
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/9/2007

% get the size of x0 so we can reshape
% later.
sx = size(x0);

% total number of derivatives we will need to take
nx = numel(x0);

HD = zeros(1,nx);
err = HD;
finaldelta = HD;
for ind = 1:nx
  [HD(ind),err(ind),finaldelta(ind)] = PairCorr.derivest( ...
    @(xi) fun(PairCorr.swapelement(x0,ind,xi)), ...
    x0(ind),'deriv',2,'vectorized','no');
end

end % mainline function end

% =======================================
%      sub-functions
% =======================================
function vec = swapelement(vec,ind,val)
% swaps val as element ind, into the vector vec
vec(ind) = val;

end % sub-function end

% -----------------------------------------------------------------------------

function [hess,err] = hessian(fun,x0)
% hessian: estimate elements of the Hessian matrix (array of 2nd partials)
% usage: [hess,err] = hessian(fun,x0)
%
% Hessian is NOT a tool for frequent use on an expensive
% to evaluate objective function, especially in a large
% number of dimensions. Its computation will use roughly
% O(6*n^2) function evaluations for n parameters.
%
% arguments: (input)
%  fun - SCALAR analytical function to differentiate.
%        fun must be a function of the vector or array x0.
%        fun does not need to be vectorized.
%
%  x0  - vector location at which to compute the Hessian.
%
% arguments: (output)
%  hess - nxn symmetric array of second partial derivatives
%        of fun, evaluated at x0.
%
%  err - nxn array of error estimates corresponding to
%        each second partial derivative in hess.
%
%
% Example usage:
%  Rosenbrock function, minimized at [1,1]
%  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
%  [h,err] = hessian(rosen,[1 1])
%  h =
%           842         -420
%          -420          210
%  err =
%    1.0662e-12   4.0061e-10
%    4.0061e-10   2.6654e-13
%
%
% Example usage:
%  cos(x-y), at (0,0)
%  Note: this hessian matrix will be positive semi-definite
%
%  hessian(@(xy) cos(xy(1)-xy(2)),[0 0])
%  ans =
%           -1            1
%            1           -1
%
%
% See also: derivest, gradient, gradest, hessdiag
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/10/2007

% parameters that we might allow to change
params.StepRatio = 2.0000001;
params.RombergTerms = 3;

% get the size of x0 so we can reshape
% later.
sx = size(x0);

% was a string supplied?
if ischar(fun)
  fun = str2func(fun);
end

% total number of derivatives we will need to take
nx = length(x0);

% get the diagonal elements of the hessian (2nd partial
% derivatives wrt each variable.)
[hess,err] = PairCorr.hessdiag(fun,x0);

% form the eventual hessian matrix, stuffing only
% the diagonals for now.
hess = diag(hess);
err = diag(err);
if nx<2
  % the hessian matrix is 1x1. all done
  return
end

% get the gradient vector. This is done only to decide
% on intelligent step sizes for the mixed partials
[grad,graderr,stepsize] = PairCorr.gradest(fun,x0);

% Get params.RombergTerms+1 estimates of the upper
% triangle of the hessian matrix
dfac = params.StepRatio.^(-(0:params.RombergTerms)');
for i = 2:nx
  for j = 1:(i-1)
    dij = zeros(params.RombergTerms+1,1);
    for k = 1:(params.RombergTerms+1)
      dij(k) = fun(x0 + PairCorr.swap2(zeros(sx),i, ...
        dfac(k)*stepsize(i),j,dfac(k)*stepsize(j))) + ...
        fun(x0 + PairCorr.swap2(zeros(sx),i, ...
        -dfac(k)*stepsize(i),j,-dfac(k)*stepsize(j))) - ...
        fun(x0 + PairCorr.swap2(zeros(sx),i, ...
        dfac(k)*stepsize(i),j,-dfac(k)*stepsize(j))) - ...
        fun(x0 + PairCorr.swap2(zeros(sx),i, ...
        -dfac(k)*stepsize(i),j,dfac(k)*stepsize(j)));

    end
    dij = dij/4/prod(stepsize([i,j]));
    dij = dij./(dfac.^2);

    % Romberg extrapolation step
    [hess(i,j),err(i,j)] =  PairCorr.rombextrap(params.StepRatio,dij,[2 4]);
    hess(j,i) = hess(i,j);
    err(j,i) = err(i,j);
  end
end


end % mainline function end

% =======================================
%      sub-functions
% =======================================
function vec = swap2(vec,ind1,val1,ind2,val2)
% swaps val as element ind, into the vector vec
vec(ind1) = val1;
vec(ind2) = val2;

end % sub-function end

% =============================================================================

% function [density, n, m, X, Y] = calcHistogram(N, data, dataCols, res)
%     %N = num rows in data set
%     Xcol = dataCols(1);
%     Ycol = dataCols(2);
%     subset = true(N,1);
%     [minX, maxX, minY, maxY] = PairCorr.getBounds(data, Xcol, Ycol);
%     NP = sum(subset);
%     XPosition = getX(data, Xcol);
%     YPosition = getY(data, Ycol);
%     n=linspace(minX,maxX,res); m=linspace(minY,maxY,res);
%
%     if NP>0
%         RR = round((res-1)*[...
%             (XPosition(subset)-minX)/(maxX-minX) ...
%             (YPosition(subset)-minY)/(maxY-minY) ])+1;
%         RRok = all(RR<=res,2) & all(RR>=1,2) ;
%         pxx=(maxX-minX)/res; pxy=(maxY-minY)/res;
%         pxArea=pxx * pxy;
%         density = accumarray(RR(RRok,:),1,[res,res])/pxArea;
%         X = repmat(n,res,1); Y = repmat(m',1,res);
%     else
%         density = zeros(res+1);
%         X = repmat(n,res+1,1); Y = repmat(m',1,res+1);
%     end

function [minX, maxX, minY, maxY] = getBounds(data, Xcol, Ycol)
    minX = min(data(:,Xcol));
    maxX = max(data(:,Xcol));
    minY = min(data(:,Ycol));
    maxY = max(data(:,Ycol));
end

function [L,vq] = Ripley(h,data,radius)
    [locs, range, box] = PairCorr.getCoords(h,data);
    maxrad = radius;
    dr = maxrad / 200;
    r = dr:dr:maxrad;
    [K,L] = RipleysK(locs,r,box,0);
    [Knew,Lnew] = ripleykfunction(locs,r,box,0);
    Lrand = PairCorr.ripleysimulation(100,r,locs,box);
    figure;
    hold on
    plot(r,L,'b-')
    plot(radius,mean(Lrand,1)+2*std(Lrand,1),'r-')
    plot(radius,mean(Lrand,1)-2*std(Lrand,1),'r-')
    hold off
    [distnew,~] = ginput(1);
    N = ripleykperpoint(locs, distnew, box, 0);
    figure;scatter(locs(:,1),locs(:,2),20,N);

    [ripX,ripY] = meshgrid(range{1},range{2});
    vq = griddata(locs(:,1),locs(:,2),N,ripX,ripY,'v4');
    figure;imagesc(vq); axis xy
end

function Lrand = ripleysimulation(numSimulations,radius,coords,box)

    Nobjects = size(coords,1);

    Lrand = zeros(numSimulations,length(radius));
    for r = 1: numSimulations
        randX = random('unif',min(coords(:,1)),max(coords(:,1)),[Nobjects,1]);
        randY = random('unif',min(coords(:,2)),max(coords(:,2)),[Nobjects,1]);
        rand_centers = [randX,randY];
        if ~isempty(rand_centers)
            [~,Lrand(r,:)] = ripleykfunction(rand_centers,radius,box,0);
        else
            Lrand(r,:) = NaN;
        end
    end
    % This gives a 95% C.I.
    % LrandOut{iROI} = Lrand;
    RadiusOut = radius;
    SDPlus = (mean(Lrand,1)+2*std(Lrand,1))';
    SDMinus = (mean(Lrand,1)-2*std(Lrand,1))';
end

function [coords, range, box] = getCoords(h,data)
    xyLim = getPosition(h);
    xmin = xyLim(1);
    ymin = xyLim(2);
    width = xyLim(3);
    height = xyLim(4);
    xmax = xmin + width;
    ymax = ymin + height;

    dx = 10; %nm
    dy = 10; %nm
    Xrange = xmin:dx:xmax;
    Yrange = ymin:dy:ymax;
    range = {};
    range{1} = Xrange;
    range{2} = Yrange;
    XPosition = data(:,1);
    YPosition = data(:,2);
    isROI = XPosition>xmin & XPosition<xmax & YPosition>ymin & YPosition<ymax;
    X = XPosition(isROI);
    Y = YPosition(isROI);
    coords = [X Y];
    box = [xmin, xmax, ymin, ymax];
end


function N = localisation_density(locs, r)
    DIST = createDistanceMatrix(locs,locs);
    DIST = sort(DIST);
    N = zeros(size(locs,1),1);
    for j = 1: size(locs,1)
        N(j) = (length(find(DIST(2:end,j)<r)));
    end
end

% -----------------------------------------------------------------------------

function [Hout Xbins Ybins] = hist2d(D, varargin) %Xn, Yn, Xrange, Yrange)
%HIST2D 2D histogram
%
% [H XBINS YBINS] = HIST2D(D, XN, YN, [XLO XHI], [YLO YHI])
% [H XBINS YBINS] = HIST2D(D, 'display' ...)
%
% HIST2D calculates a 2-dimensional histogram and returns the histogram
% array and (optionally) the bins used to calculate the histogram.
%
% Inputs:
%     D:         N x 2 real array containing N data points or N x 1 array
%                 of N complex values
%     XN:        number of bins in the x dimension (defaults to 20)
%     YN:        number of bins in the y dimension (defaults to 20)
%     [XLO XHI]: range for the bins in the x dimension (defaults to the
%                 minimum and maximum of the data points)
%     [YLO YHI]: range for the bins in the y dimension (defaults to the
%                 minimum and maximum of the data points)
%     'display': displays the 2D histogram as a surf plot in the current
%                 axes
%
% Outputs:
%     H:         2D histogram array (rows represent X, columns represent Y)
%     XBINS:     the X bin edges (see below)
%     YBINS:     the Y bin edges (see below)
%
% As with histc, h(i,j) is the number of data points (dx,dy) where
% x(i) <= dx < x(i+1) and y(j) <= dx < y(j+1). The last x bin counts
% values where dx exactly equals the last x bin value, and the last y bin
% counts values where dy exactly equals the last y bin value.
%
% If D is a complex array, HIST2D splits the complex numbers into real (x)
% and imaginary (y) components.
%
% Created by Amanda Ng on 5 December 2008

% Modification history
%   25 March 2009 - fixed error when min and max of ranges are equal.
%   22 November 2009 - added display option; modified code to handle 1 bin

    % PROCESS INPUT D
    if nargin < 1 %check D is specified
        error 'Input D not specified'
    end

    Dcomplex = false;
    if ~isreal(D) %if D is complex ...
        if isvector(D) %if D is a vector, split into real and imaginary
            D=[real(D(:)) imag(D(:))];
        else %throw error
            error 'D must be either a complex vector or nx2 real array'
        end
        Dcomplex = true;
    end

    if (size(D,1)<size(D,2) && size(D,1)>1)
        D=D';
    end

    if size(D,2)~=2;
        error('The input data matrix must have 2 rows or 2 columns');
    end

    % PROCESS OTHER INPUTS
    var = varargin;

    % check if DISPLAY is specified
    index = find(strcmpi(var,'display'));
    if ~isempty(index)
        display = true;
        var(index) = [];
    else
        display = false;
    end

    % process number of bins
    Xn = 20; %default
    Xndefault = true;
    if numel(var)>=1 && ~isempty(var{1}) % Xn is specified
        if ~isscalar(var{1})
            error 'Xn must be scalar'
        elseif var{1}<1 %|| ~isinteger(var{1})
            error 'Xn must be an integer greater than or equal to 1'
        else
            Xn = var{1};
            Xndefault = false;
        end
    end

    Yn = 20; %default
    Yndefault = true;
    if numel(var)>=2 && ~isempty(var{2}) % Yn is specified
        if ~isscalar(var{2})
            error 'Yn must be scalar'
        elseif var{2}<1 %|| ~isinteger(var{2})
            error 'Xn must be an integer greater than or equal to 1'
        else
            Yn = var{2};
            Yndefault = false;
        end
    end

    % process ranges
    if numel(var) < 3 || isempty(var{3}) %if XRange not specified
        Xrange=[min(D(:,1)),max(D(:,1))]; %default
    else
        if nnz(size(var{3})==[1 2]) ~= 2 %check is 1x2 array
            error 'XRange must be 1x2 array'
        end
        Xrange = var{3};
    end
    if Xrange(1)==Xrange(2) %handle case where XLO==XHI
        if Xndefault
            Xn = 1;
        else
            Xrange(1) = Xrange(1) - floor(Xn/2);
            Xrange(2) = Xrange(2) + floor((Xn-1)/2);
        end
    end

    if numel(var) < 4 || isempty(var{4}) %if XRange not specified
        Yrange=[min(D(:,2)),max(D(:,2))]; %default
    else
        if nnz(size(var{4})==[1 2]) ~= 2 %check is 1x2 array
            error 'YRange must be 1x2 array'
        end
        Yrange = var{4};
    end
    if Yrange(1)==Yrange(2) %handle case where YLO==YHI
        if Yndefault
            Yn = 1;
        else
            Yrange(1) = Yrange(1) - floor(Yn/2);
            Yrange(2) = Yrange(2) + floor((Yn-1)/2);
        end
    end

    % SET UP BINS
    Xlo = Xrange(1) ; Xhi = Xrange(2) ;
    Ylo = Yrange(1) ; Yhi = Yrange(2) ;
    if Xn == 1
        XnIs1 = true;
        Xbins = [Xlo Inf];
        Xn = 2;
    else
        XnIs1 = false;
        Xbins = linspace(Xlo,Xhi,Xn) ;
    end
    if Yn == 1
        YnIs1 = true;
        Ybins = [Ylo Inf];
        Yn = 2;
    else
        YnIs1 = false;
        Ybins = linspace(Ylo,Yhi,Yn) ;
    end

    Z = linspace(1, Xn+(1-1/(Yn+1)), Xn*Yn);

    % split data
    Dx = floor((D(:,1)-Xlo)/(Xhi-Xlo)*(Xn-1))+1;
    Dy = floor((D(:,2)-Ylo)/(Yhi-Ylo)*(Yn-1))+1;
    Dz = Dx + Dy/(Yn) ;

    % calculate histogram
    h = reshape(histc(Dz, Z), Yn, Xn);

    if nargout >=1
        Hout = h;
    end

    if XnIs1
        Xn = 1;
        Xbins = Xbins(1);
        h = sum(h,1);
    end
    if YnIs1
        Yn = 1;
        Ybins = Ybins(1);
        h = sum(h,2);
    end

    % DISPLAY IF REQUESTED
    if ~display
        return
    end

    [x y] = meshgrid(Xbins,Ybins);
    dispH = h;

    % handle cases when Xn or Yn
    if Xn==1
        dispH = padarray(dispH,[1 0], 'pre');
        x = [x x];
        y = [y y];
    end
    if Yn==1
        dispH = padarray(dispH, [0 1], 'pre');
        x = [x;x];
        y = [y;y];
    end

    surf(x,y,dispH);
    colormap(jet);
    if Dcomplex
        xlabel real;
        ylabel imaginary;
    else
        xlabel x;
        ylabel y;
    end

end

% -----------------------------------------------------------------------------

function y = exponential_and_gaussian(P, r)
xi = P(1);
A = P(2);
s = P(3);
rho = P(4)*1e-6;
if length(P)<5, C = 0; else C = P(5); end

A2 = 1/2/pi/s^2/rho;
y = A2*exp(-r.^2/2/s^2) + A*exp(-r/xi)+C;

end

% -----------------------------------------------------------------------------

function y = exponential_and_cosine(P, r)
xi = P(1);
A = P(2);
r0 = P(3);
if length(P)<4, C = 1; else C = P(4); end

y = C + A*exp(-r/xi).*cos((pi*r)/2/r0);

end

% -----------------------------------------------------------------------------

function y = exponential(P, r)
xi = P(1);
A = P(2);
if length(P)<3, C = 0; else C = P(3); end

y = A*exp(-r/xi)+C;

end

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
