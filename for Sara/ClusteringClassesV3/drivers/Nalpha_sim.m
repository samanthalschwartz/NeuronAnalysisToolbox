clear all
close all

Quiet = true;
%Random = true;
Collapse = false;

A_range = 1 : 1 : 3;
N_range = 3 : 1 : 5;

nA_range = numel(A_range);
nN_range = numel(N_range);

N_density    = cell(nA_range, nN_range);
R_density    = cell(nA_range, nN_range);
N_densities  = cell(nA_range, nN_range);
R_densities  = cell(nA_range, nN_range);
N_numclust_3 = cell(nA_range, nN_range);
R_numclust_3 = cell(nA_range, nN_range);
N_frac_clust = cell(nA_range, nN_range);
R_frac_clust = cell(nA_range, nN_range);
N_radii      = cell(nA_range, nN_range);
R_radii      = cell(nA_range, nN_range);
N_n_clusters = cell(nA_range, nN_range);
R_n_clusters = cell(nA_range, nN_range);
N_n_mult_obj = cell(nA_range, nN_range);
R_n_mult_obj = cell(nA_range, nN_range);

%n_sims = 100;
n_sims = 2;

Area = 1000^2;   % ROI area

for j = 1 : 2
   if j == 1
      Random = true;
   else
      Random = false;
   end

   k = 0;
   for A = A_range
      k = k + 1;

      l = 0;
      for N = N_range
         l = l + 1;

         %nC = [];
         %n_points = [];
         %singlet_fraction = [];
         %min_c2c_dist = [];
         %min_c2c_dists = [];
         %equiv_radii = [];
         %density = [];
         densities = [];
         %numclust_1 = [];
         numclust_3 = [];
         frac_clust = [];
         radii = [];
         n_clusters = [];
         n_mult_obj = [];
         for i = 1 : n_sims
            [results, results_o, results_c, results_t, results_simulated] = ...
               Nalpha_sim1(A, N, Quiet, Random, Collapse);
        
            %nC = [nC, results_o.nC];
            %n_points = [n_points, results_o.n_points];
            %singlet_fraction = [singlet_fraction, results_o.singlet_fraction];
            %equiv_radii = [equiv_radii, results_o.equiv_radii];
            %min_c2c_dist = [min_c2c_dist, results_o.min_c2c_dist];
            %min_c2c_dists = [min_c2c_dists, results_o.min_c2c_dists];
            %density = [density, sum(results_o.numclust) / Area];
            densities = [densities, results_o.n_pts_per_area];
            %numclust_1 = [numclust_1, results_o.numclust(1)];
            numclust_3 = [numclust_3, results_o.numclust(3)];
            frac_clust = [frac_clust, ...
                          results_o.n_clustered / results_o.n_points];
            radii      = [radii, results_o.equiv_radii];
            n_clusters = [n_clusters, sum(results_o.numclust)];
            n_mult_obj = [n_mult_obj, results_o.n_pts];
%           fprintf('Sim %3d: E = %2d, minPts = %2d, n_multi = %2d\n', ...
%                   i, E, N, results_o.numclust(3));
         end

%        fprintf('\n');
%        fprintf('mean(nC) = %.3f\n', mean(nC));
%        fprintf('mean(n_points) = %.3f\n', mean(n_points));
%        fprintf('mean(singlet_fraction) = %.3f\n', mean(singlet_fraction));
%        fprintf('mean(equiv_radii) = %.3f\n', mean(equiv_radii));
%        fprintf('mean(min_c2c_dist) = %.3f\n', mean(min_c2c_dist));
%        fprintf('mean(min_c2c_dists) = %.3f\n', mean(min_c2c_dists));

         if j == 1
            %R_singlet_fraction = singlet_fraction;
            %R_equiv_radii = equiv_radii;
            %R_density{k, l} = density;
            R_densities{k, l} = densities;
            %R_numdensity_1 = numclust_1 / Area;
            %R_numdensity_3 = numclust_3 / Area;
            R_numclust_3{k, l} = numclust_3;
            R_frac_clust{k, l} = frac_clust;
            R_radii{k, l} = radii;
            R_n_clusters{k, l} = n_clusters;
            R_n_mult_obj{k, l} = n_mult_obj;
         else
            %N_singlet_fraction = singlet_fraction;
            %N_equiv_radii = equiv_radii;
            %N_density{k, l} = density;
            N_densities{k, l} = densities;
            %N_numdensity_1 = numclust_1 / Area;
            %N_numdensity_3 = numclust_3 / Area;
            N_numclust_3{k, l} = numclust_3;
            N_frac_clust{k, l} = frac_clust;
            N_radii{k, l} = radii;
            N_n_clusters{k, l} = n_clusters;
            N_n_mult_obj{k, l} = n_mult_obj;
         end
      end
   end
    
end

N_numclust_3_mean   = cellfun(@mean,   N_numclust_3);
N_numclust_3_median = cellfun(@median, N_numclust_3);
R_numclust_3_mean   = cellfun(@mean,   R_numclust_3);
R_numclust_3_median = cellfun(@median, R_numclust_3);
N_densities_mean    = cellfun(@mean,   N_densities);
N_densities_median  = cellfun(@median, N_densities);
R_densities_mean    = cellfun(@mean,   R_densities);
R_densities_median  = cellfun(@median, R_densities);
N_frac_clust_mean   = cellfun(@mean,   N_frac_clust);
N_frac_clust_median = cellfun(@median, N_frac_clust);
R_frac_clust_mean   = cellfun(@mean,   R_frac_clust);
R_frac_clust_median = cellfun(@median, R_frac_clust);
N_n_clusters_mean   = cellfun(@mean,   N_n_clusters);
N_n_clusters_median = cellfun(@median, N_n_clusters);
R_n_clusters_mean   = cellfun(@mean,   R_n_clusters);
R_n_clusters_median = cellfun(@median, R_n_clusters);
N_n_mult_obj_mean   = cellfun(@mean,   N_n_mult_obj);
N_n_mult_obj_median = cellfun(@median, N_n_mult_obj);
R_n_mult_obj_mean   = cellfun(@mean,   R_n_mult_obj);
R_n_mult_obj_median = cellfun(@median, R_n_mult_obj);

figure();
hold on
surf(A_range, N_range, N_numclust_3_mean');
title('Non-random mean per ROI');
xlabel('\alpha');
ylabel('minPts');
zlabel('# multi-clusters');
hold off

%figure();
%hold on
%surf(A_range, N_range, N_numclust_3_median');
%title('Non-random median per ROI');
%xlabel('\alpha');
%ylabel('minPts');
%zlabel('# multi-clusters');
%hold off

figure();
hold on
surf(A_range, N_range, R_numclust_3_mean');
title('Random mean per ROI');
xlabel('\alpha');
ylabel('minPts');
zlabel('# multi-clusters');
hold off

%figure();
%hold on
%surf(A_range, N_range, R_numclust_3_median');
%title('Random median per ROI');
%xlabel('\alpha');
%ylabel('minPts');
%zlabel('# multi-clusters');
%hold off

save('n_multi_SimNV.mat', ...
     'A_range', 'N_range', 'N_numclust_3_mean', 'N_numclust_3_median', ...
     'N_densities_mean',  'N_densities_median',                        ...
     'N_frac_clust_mean', 'N_frac_clust_median',                       ...
     'N_n_clusters_mean', 'N_n_clusters_median',                       ...
     'N_n_mult_obj_mean', 'N_n_mult_obj_median');
save('n_multi_SimRV.mat', ...
     'A_range', 'N_range', 'R_numclust_3_mean', 'R_numclust_3_median', ...
     'R_densities_mean',  'R_densities_median',                        ...
     'R_frac_clust_mean', 'R_frac_clust_median',                       ...
     'R_n_clusters_mean', 'R_n_clusters_median',                       ...
     'R_n_mult_obj_mean', 'R_n_mult_obj_median');
