clear all
close all

% Input directory.
indir = 'data';

pc = PairCorr();

pc.Results = 'results'; % results directory---this needs to exist beforehand
pc.Fig_ext = 'png';     % figure extension
pc.ROIs = true;         % use ROIs
pc.ROI_size = 3000;     % x and y width of square fixed sized ROI (nm)
pc.Rmax_axis = 500;     % Sets plotting limit if > 0
% Conversion from pixels to nm
pc.Pixel2nm = 104;
% Histogram bin size for pairwise correlation---this is the number of pixels
% per bin over which correlation statistics are collected.
pc.Hist_bin_size = 104 / 20;
pc.HSET = true;         % perform H-SET?
pc.Verbose = false;     % verbose output and extra saved .mat files
pc.Auto = false;        % perform auto- AND cross-correlations together(PC_SR)?
pc.Veatch = false;      % in addition, use the Veatch code directly?

%file = dir(fullfile(indir, 'ResultsOnly_*.mat'));
%n_files = numel(file);

%for i = 1 : n_files
   %infile1 = file(i).name;
   %infile2 = '';
   %basefile = regexprep(file(i).name, '.*ResultsOnly_(.*)\.mat', '$1');
   infile1 = 'NewSample_Cell2_Label1_secondtry#0001-2016-4-19-16-11-13_Results.mat';
   infile2 = 'Cell2_secondtry_label2#0001-2016-4-20-1-7-11_Results.mat';
   basefile = 'Farzin';

   load(fullfile(indir, infile1));
   if exist('SRD', 'var')
      SRD1 = SRD.Results;
   else
      SRD1 = R647;
      SRD2 = R488;
   end
   if ~isempty(infile2)
      load(fullfile(indir, infile2));
      if exist('SRD', 'var')
         SRD2 = SRD.Results;
      end
   end

   pc.PC_SR(SRD1, SRD2, basefile);
%  pc.PCC_SR(SRD1, SRD2, basefile);
%  pc.PAC_SR(SRD1, [basefile, '_L1']);
%  if exist('SRD2', 'var')
%     pc.PAC_SR(SRD2, [basefile, '_L2']);
%  end
%end
