clear all
close all

%A_range = 1 : 1 : 10;
%N_range = 1 : 1 : 15;
A_range = 1 : 1 : 3;
N_range = 3 : 1 : 5;

nA_range = numel(A_range);
nN_range = numel(N_range);

density    = cell(nA_range, nN_range);
densities  = cell(nA_range, nN_range);
n_clusters = cell(nA_range, nN_range);
n_multi    = cell(nA_range, nN_range);
sing_fract = cell(nA_range, nN_range);
min_c2c_ds = cell(nA_range, nN_range);
c2c_dists  = cell(nA_range, nN_range);
frac_clust = cell(nA_range, nN_range);
radii      = cell(nA_range, nN_range);
n_mult_obj = cell(nA_range, nN_range);

n_collapsed_per_singlet = [];
cutoff_eps = [];

SRc = SRcluster();
%SRc.PvalueStatistics = true;
SRc.ShrinkFactor = 0.5;
SRc.LoS = 0.01;
SRc.A_ROI = (1000*1000);
SRc.PlotFigures = false;
SRc.Printing = false;
%SRc.Algorithm = 'Hierarchical';
%SRc.Algorithm = 'DBSCAN_Daszykowski';   SRc.minPts = 3;
SRc.Algorithm = 'Voronoi';   SRc.E = 0;   SRc.Valgorithm = 2;

load('ROI_collapsed.mat');

%for j = 1 : numel(DataDirs)
%  DataDir = DataDirs{j};
%  Files = dir(fullfile(DataDir, '*_ROI.mat'));
%  FileName = fullfile(DataDir, Files.name)
%  load(FileName);
%   
%  %if isnan(Sigma_Reg)
      Sigma_Reg = [0, 0];
%  %end
%
%  for i = 1 : numel(RoI)
%     X = double(RoI{i}.X);   % nm
%     Y = double(RoI{i}.Y);   % nm
%     X_STD = double(RoI{i}.X_STD);   % nm
%     Y_STD = double(RoI{i}.Y_STD);   % nm
%
%     SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
%
%     SRc.Cutoff = E_range;
%     cutoff_eps = [cutoff_eps, SRc.chooseCutoff()];
   for i = 1 : numel(XY_SR)

      SRc.set_XY_sigma(Sigma_Reg, XY_orig{i}, Sigma_orig{i}, XY_SR{i}, ...
                       Sigma_SR{i}, Combined{i});

      k = 0;
      for A = A_range
         k = k + 1;
         SRc.Alpha = A;   % alpha

         l = 0;
         for N = N_range
            l = l + 1;
            SRc.minPts = N;

            [results, analysisFigs] = SRc.analyzeSRclusters();

            if isempty(results.singlets)
               ncollapsed_per_singlet = NaN;
            else
               ncollapsed_per_singlet = cellfun(@numel, results.singlets);
            end
            n_collapsed_per_singlet = ...
               [n_collapsed_per_singlet, ncollapsed_per_singlet];

            ROI_cluster_density = ...
               sum(results.cluster_numobjs) / sum(results.cluster_areas);
            ROI_cluster_densities = ...
               results.cluster_numobjs ./ results.cluster_areas;
            ROI_cluster_density(isnan(ROI_cluster_density)) = [];
            ROI_cluster_densities(isnan(ROI_cluster_densities)) = [];
            density{k, l}    = [density{k, l}, ROI_cluster_density];
            densities{k, l}  = [densities{k, l}, ROI_cluster_densities];
            n_clusters{k, l} = [n_clusters{k, l}, sum(results.numclust)];
            n_multi{k, l}    = [n_multi{k, l}, results.numclust(3)];
            sing_fract{k, l} = [sing_fract{k, l}, ...
                                results.singles_per_total_clusters];
            min_c2c_ds{k, l} = [min_c2c_ds{k, l}, min(results.min_c2c_dists)];
            c2c_dists{k, l}  = [c2c_dists{k, l}, results.min_c2c_dists];
            frac_clust{k, l} = [frac_clust{k, l}, ...
               1 - results.numobjs(1) / sum(results.numobjs)];
            radii{k, l}      = [radii{k, l}, ...
                                results.radius_per_2orMore_cluster];
            n_mult_obj{k, l} = [n_mult_obj{k, l}, ...
                                results.numobjs_per_multiple_cluster];
%           fprintf('ROI %1d: A = %2d, minPts = %2d, n_multi = %2d\n', ...
%                   i, A, N, results.numclust(3));
         end
      end
   end
%end

n_multi_mean      = cellfun(@mean, n_multi);
n_multi_median    = cellfun(@median, n_multi);
densities_mean    = cellfun(@mean, densities);
densities_median  = cellfun(@median, densities);
frac_clust_mean   = cellfun(@mean, frac_clust);
frac_clust_median = cellfun(@median, frac_clust);
n_clusters_mean   = cellfun(@mean, n_clusters);
n_clusters_median = cellfun(@median, n_clusters);
n_mult_obj_mean   = cellfun(@mean, n_mult_obj);
n_mult_obj_median = cellfun(@median, n_mult_obj);

figure();
hold on
surf(A_range, N_range, n_multi_mean');
title('mean per ROI');
xlabel('\alpha');
ylabel('minPts');
zlabel('# multi-clusters');
hold off

%figure();
%hold on
%surf(A_range, N_range, n_multi_median');
%title('median per ROI');
%xlabel('\alpha');
%ylabel('minPts');
%zlabel('# multi-clusters');
%hold off

figure();
hold on
surf(A_range, N_range, n_clusters_mean');
title('mean per ROI');
xlabel('\alpha');
ylabel('minPts');
zlabel('# clusters including singlets');
hold off

%figure();
%hold on
%surf(A_range, N_range, n_clusters_median');
%title('median per ROI');
%xlabel('\alpha');
%ylabel('minPts');
%zlabel('# clusters including singlets');
%hold off

figure();
hold on
surf(A_range, N_range, n_mult_obj_mean');
title('mean per ROI');
xlabel('\alpha');
ylabel('minPts');
zlabel('# localizations per cluster');
hold off

%figure();
%hold on
%surf(A_range, N_range, n_mult_obj_median');
%title('median per ROI');
%xlabel('\alpha');
%ylabel('minPts');
%zlabel('# localizations per cluster');
%hold off

save('n_multi_ExperV.mat', 'A_range', 'N_range', ...
     'n_multi_mean', 'n_multi_median',           ...
     'densities_mean', 'densities_median',       ...
     'frac_clust_mean', 'frac_clust_median',     ...
     'n_clusters_mean', 'n_clusters_median',     ...
     'n_mult_obj_mean', 'n_mult_obj_median');
