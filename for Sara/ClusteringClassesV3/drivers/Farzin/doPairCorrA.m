function [results_pcc, resultsRC_pcc] = ...
   doPairCorrA(n_ROIs, RoI, desc, results, xy1, xy2)
% Pair cross-correlation for each ROI.

   pc = PairCorr();

   pc.Results = results; % results directory---this needs to exist beforehand
   pc.Fig_ext = 'png';   % figure extension
   pc.ROIs = true;       % use ROIs
   pc.ROI_size = 3000;   % x and y width of square fixed sized ROI (nm)
   pc.Rmax_axis = 500;   % Sets plotting limit if > 0
   % Conversion from pixels to nm (only needed when no ROIs)
   %pc.Pixel2nm = 104;
   % Histogram bin size for pairwise correlation---this is the number of pixels
   % per bin over which correlation statistics are collected.
   pc.Hist_bin_size = 104 / 20;
   pc.HSET = false;      % perform H-SET?
   pc.Verbose = false;   % verbose output and extra saved .mat files
   pc.Auto = false;      % perform auto- AND cross-correlations together(PC_SR)?
   pc.Veatch = false;    % in addition, use the Veatch code directly?

   ROIs_combined = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      fprintf('ROI %d\n\n', i);

      pc.ROI = [RoI{i}.ROI(1), RoI{i}.ROI(3), ...
                RoI{i}.ROI(2) - RoI{i}.ROI(1), RoI{i}.ROI(4) - RoI{i}.ROI(3)];
      txt = sprintf('%s_ROI%d', desc, i);
%     XY1 = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
%     XY2 = [ RoI{i}.X{2}, RoI{i}.Y{2} ];
      XY1 = xy1{i};
      XY2 = xy2{i};

      ROIs_combined{i}.ROI = RoI{i}.ROI;
      ROIs_combined{i}.XY1 = XY1;
      ROIs_combined{i}.XY2 = XY2;

      results_pcc{i} = pc.pair_correlation(pc.Hist_bin_size, txt, XY1, XY2);
   end

   fprintf('ROI combined\n\n');
   txt = sprintf('%s_RC', desc);
   resultsRC_pcc = ...
      pc.pair_correlation_ROIcombined(pc.Hist_bin_size, txt, ...
                                      2, n_ROIs, ROIs_combined)

end
