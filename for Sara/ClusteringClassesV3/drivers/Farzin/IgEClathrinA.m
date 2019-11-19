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

% Note:
%    results_pcc{1:n_ROIs}      pair cross-correlation results
%    resultsRC_pcc              as above for ROIs combined
%    results_c{1:n_ROIs}{1:2}   clustering results by ROI and label
%    results_cs{1:n_ROIs}       cluster separations between 2 labels

load(fullfile(results_dir, 'ROIs.mat'));
load(fullfile(results_dir, sprintf('%s_results.mat', desc)));

centers1 = cell(1, n_ROIs);
centers2 = cell(1, n_ROIs);
for i = 1 : n_ROIs
   centers1{i} = results_c{i}{1}.centers';
   centers2{i} = results_c{i}{2}.centers';
end
txt = sprintf('%s_centers', desc);
[results_pccA, resultsRC_pccA] = ...
   doPairCorrA(n_ROIs, RoI, txt, results_dir, centers1, centers2);

c = Clustering();
for i = 1 : n_ROIs
   xy = { [RoI{i}.X{1}, RoI{i}.Y{1}], [RoI{i}.X{2}, RoI{i}.Y{2}] };
   CC = { results_c{i}{1}.C, results_c{i}{2}.C };
   Centers = { results_c{i}{1}.centers results_c{i}{2}.centers };
   PtsI = { results_c{i}{1}.ptsI, results_c{i}{2}.ptsI };
   shrinkFactor = 0.5;
   txt = sprintf('%s_ROI%d_Lcombined', desc, i);

   clusterFig = c.plotLabelClusters(xy, CC, Centers, PtsI, txt, shrinkFactor);
   saveas(clusterFig, fullfile(results_dir, sprintf('%s.png', txt)));
end
