clear all
close all

% Input directory.
indir = 'EMimages';

pc = PairCorr();

pc.Results = 'results'; % results directory---this needs to exist beforehand
%pc.Results = indir;     % results directory---this needs to exist beforehand
pc.Fig_ext = 'pdf';     % figure extension
pc.ROIs = true;         % use ROIs
pc.ROI_size = 500;      % x and y width of square fixed sized ROI (nm)
%pc.Rmax_axis = pc.ROI_size/2;   % Sets plotting limit if > 0
% Conversion from pixels to nm
pc.Pixel2nm = 104;   
% Histogram bin size for pairwise correlation---this is the number of pixels
% per bin over which correlation statistics are collected.
%pc.Hist_bin_size = 5;
pc.Hist_bin_size = 104 / 20;
pc.HSET = false;        % perform H-SET?
pc.Verbose = false;     % verbose output and extra saved .mat files
pc.Auto = true;         % perform auto- AND cross-correlations together(PC_SR)?
pc.Veatch = false;      % in addition, use the Veatch code directly?

file = dir(fullfile(indir, '*_6.txt'));
n_files = numel(file);

for i = 1 : n_files
   basefile = regexprep(file(i).name, '_6.txt', '');
   infile1 = [basefile, '_6.txt'];
   infile2 = [basefile, '_12.txt'];

   [x, y] = textread(fullfile(indir, infile1), '%*u %u %u', 'headerlines', 1);
   XY1 = [x, y];
   [x, y] = textread(fullfile(indir, infile2), '%*u %u %u', 'headerlines', 1);
   XY2 = [x, y];

   pc.PC_SR(XY1, XY2, basefile);
%  pc.PCC_SR(XY1, XY2, basefile);
%  pc.PAC_SR(XY1, [basefile, '_L1']);
%  pc.PAC_SR(XY2, [basefile, '_L2']);
end
