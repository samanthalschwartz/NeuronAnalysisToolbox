clear all
close all

addpath('Y:\MJW\IgE-Clathrin');
addpath('Y:\MJW\cluster');

data_dir = 'Y:\Farzin\IgE-Clathrin\ProjectReport_Sep30_2017\5min_Overlay_Sep30';
results_dir = 'Y:\Farzin\IgE-Clathrin\ProjectReport_Sep30_2017\5min_Overlay_Sep30\results_wholeROI';
desc = '5min_Overlay_Sep30-Cell10';

%data_dir = 'data';
%results_dir = 'results';
%desc = '1min_Overlay_Sep29-Cell6';

pixel2nm = 104;   % conversion factor from pixels to nm

% Create results_dir if it does not already exist.
if ~isdir(results_dir)
   mkdir(results_dir)
end

% Note:
%    results_pcc{1:n_ROIs}      pair cross-correlation results
%    resultsRC_pcc              as above for ROIs combined
%    results_c{1:n_ROIs}{1:2}   clustering results by ROI and label
%    results_cs{1:n_ROIs}       cluster separations between 2 labels

createROIs(data_dir, results_dir, desc, pixel2nm)
[results_pcc, resultsRC_pcc, results_c, results_cs] = ...
   doAnalysis(data_dir, results_dir, desc, pixel2nm);
