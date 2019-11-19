% Read in the RegistrationAnalysis object computed previously (by regis.m, for
% example), load the localization data (R488 and R647 here), map the right
% channel into the left, then perform pair correlation on the combination.
% Required: Clustering.m, PairCorr.m, ROITools.m, SRcluster.m
% Michael Wester (12/5/16)

clear all
close all

OoB = true;       % eliminate out-of-bounds coordinates
verbose = true;   % make a plot of original and mapped data

% The RA object is assumed saved in RegAnal.mat (see regis.m).
load(fullfile('mat', 'RegAnal.mat'));
SXS = RA.SensorXSplit;      % channel 2 x-coordinates should be > SXS
SX  = RA.SensorSizeX / 2;   % x-coordinates should approximately be in [0, SX]
SY  = RA.SensorSizeY;       % y-coordinates should approximately be in [0, SY]
% Construct the mapping from the right channel (647) to the left (488).
% Note that mapping occurs in absolute coordinates.
% Optimal map:
M = RA.getOptimalMapPixels();
% Maps by algorithm:
%    (1) Null
%    (2) GlobalAffine
%    (3) LocalAffine
%    (4) SmoothAffine
%    (5) LWM
%    (6) Polynomial
%    (7) NRS
%M = RA.maps(5).mapFunctionPixels;

% Input directory.
indir = 'data';

pc = PairCorr();

pc.Results = 'results'; % results directory---this needs to exist beforehand
pc.Fig_ext = 'png';   % figure extension
pc.ROIs = true;       % use ROIs
pc.ROI_size = 3000;   % x and y width of square fixed sized ROI (nm)
pc.Rmax_axis = 500;   % sets plotting limit if > 0
% Conversion from pixels to nm
pc.Pixel2nm = RA.pixelSize;  % consistent with RegistrationAnalysis
% Histogram bin size for pairwise correlation---this is the number of pixels
% per bin over which correlation statistics are collected.
pc.Hist_bin_size = 5; % 104 / 20;
pc.HSET = true;       % perform H-SET?
pc.Verbose = false;
pc.Auto = false;      % perform auto- AND cross-correlations together (PC_SR)?
pc.Veatch = false;    % in addition, use the Veatch code directly?

file = dir(fullfile(indir, 'ResultsOnly_*.mat'));
n_files = numel(file);

for i = 1 %: n_files
   infile1 = file(i).name;
   infile2 = '';
   basefile = regexprep(infile1, '.*ResultsOnly_(.*)\.mat', '$1');

   load(fullfile(indir, infile1));
   if exist('SRD', 'var')
      SRD1 = SRD.Results;
   else
      SRD1 = R488;
      SRD2 = R647;
   end
   if ~isempty(infile2)
      load(fullfile(indir, infile2));
      if exist('SRD', 'var')
         SRD2 = SRD.Results;
      end
   end

   fprintf('Total points in R488: %d;   R647: %d\n', ...
           numel(SRD1.X), numel(SRD2.X));
   if OoB
      % Eliminate out-of-bounds coordinates.
      indx = SRD1.X < 0 | SRD1.X > SX | SRD1.Y < 0 | SRD1.Y > SY;
      if sum(indx) > 0
         fprintf('Eliminated %d out-of-bounds points from R488.\n', sum(indx));
      end
      SRD1.X(indx) = [];
      SRD1.Y(indx) = [];
      SRD1.X_STD(indx) = [];
      SRD1.Y_STD(indx) = [];
      indx = SRD2.X < 0 | SRD2.X > SX | SRD2.Y < 0 | SRD2.Y > SY;
      if sum(indx) > 0
         fprintf('Eliminated %d out-of-bounds points from R647.\n', sum(indx));
      end
      SRD2.X(indx) = [];
      SRD2.Y(indx) = [];
      SRD2.X_STD(indx) = [];
      SRD2.Y_STD(indx) = [];
   end

   if verbose
      fprintf('(488): x = [%f, %f], y = [%f, %f]\n', ...
              min(SRD1.X), max(SRD1.X), min(SRD1.Y), max(SRD1.Y));
      fprintf('(647): x = [%f, %f], y = [%f, %f]\n', ...
              min(SRD2.X), max(SRD2.X), min(SRD2.Y), max(SRD2.Y));

      % Show original and mapped data.
      figure();
      hold on
      plot(SRD1.X, SRD1.Y, 'k.');
      plot(SRD2.X, SRD2.Y, 'g.');
   end
   % Adding SXS pixels shifts channel 2 x-coordinates from relative coordinates
   % with channel 2 in isolation back to the absolute coordinates used to
   % construct the map with the channels side-by-side.  The mapping then
   % transforms channel 2 into the coordinate space of channel 1.
   SRD2_XY = M([SRD2.X + SXS, SRD2.Y]);
   SRD2.X = SRD2_XY(:, 1);
   SRD2.Y = SRD2_XY(:, 2);
   if OoB
      % Eliminate out-of-bounds coordinates.
      indx = SRD2.X < 0 | SRD2.X > SX | SRD2.Y < 0 | SRD2.Y > SY;
      if sum(indx) > 0
        fprintf('Eliminated %d out-of-bounds points from mapped R647.\n', ...
                 sum(indx));
      end
      SRD2.X(indx) = [];
      SRD2.Y(indx) = [];
      SRD2.X_STD(indx) = [];
      SRD2.Y_STD(indx) = [];
   end

   if verbose
      plot(SRD2.X, SRD2.Y, 'r.');
      title(basefile);
      xlabel('x (pixels)');
      ylabel('y (pixels)');
      legend('488 original', '647 original', '647 mapped');
      hold off

      fprintf('Mapping:\n');
      fprintf('(647): x = [%f, %f], y = [%f, %f]\n', ...
              min(SRD2.X), max(SRD2.X), min(SRD2.Y), max(SRD2.Y));
   end

%  pc.PC_SR(SRD1, SRD2, basefile);          % cross- + auto-correlations
%  pc.PCC_SR(SRD1, SRD2, basefile);         % cross-correlations only
%  pc.PAC_SR(SRD1, [basefile, '_L1']);      % auto-correlations only
%  if exist('SRD2', 'var')
%     pc.PAC_SR(SRD2, [basefile, '_L2']);   % auto-correlations only
%  end
end
