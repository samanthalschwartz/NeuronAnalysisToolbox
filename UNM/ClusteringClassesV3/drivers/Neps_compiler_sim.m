function mainSLMatt

clear all
close all

sextresults = load('sext.mat');

Quiet = true;
Random = true;
Collapse = false;

%n_sims = 100;
n_sims = 5;

A = 1000^2;   % ROI area
E = 26;       % epsilon
minPts = 3;

for j = 1 : 2
   if j == 1
      Random = true;
   else
      Random = false;
   end
    
   nC = [];
   n_points = [];
   singlet_fraction = [];
   min_c2c_dist = [];
   min_c2c_dists = [];
   equiv_radii = [];
   density = [];
   numclust_1 = [];
   numclust_3 = [];
   n_pts = [];
   for i = 1 : n_sims
      [results, results_o, results_c, results_t, results_simulated] = ...
          Neps_sim1(E, minPts, Quiet, Random, Collapse);
        
      nC = [nC, results_o.nC];
      n_points = [n_points, results_o.n_points];
      singlet_fraction = [singlet_fraction, results_o.singlet_fraction];
      equiv_radii = [equiv_radii, results_o.equiv_radii];
      min_c2c_dist = [min_c2c_dist, results_o.min_c2c_dist];
      min_c2c_dists = [min_c2c_dists, results_o.min_c2c_dists];
      density = [density, sum(results_o.numclust) / A];
      numclust_1 = [numclust_1, results_o.numclust(1)];
      numclust_3 = [numclust_3, results_o.numclust(3)];
      n_pts = [n_pts, results_o.n_pts];
   end
    
   fprintf('\n');
   fprintf('mean(nC) = %.3f\n', mean(nC));
   fprintf('mean(n_points) = %.3f\n', mean(n_points));
   fprintf('mean(singlet_fraction) = %.3f\n', mean(singlet_fraction));
   fprintf('mean(equiv_radii) = %.3f\n', mean(equiv_radii));
   fprintf('mean(min_c2c_dist) = %.3f\n', mean(min_c2c_dist));
   fprintf('mean(min_c2c_dists) = %.3f\n', mean(min_c2c_dists));
    
   if j == 1
      R_singlet_fraction = singlet_fraction;
      R_equiv_radii = equiv_radii;
      R_density = density;
      R_n_pts = n_pts;
      R_numdensity_1 = numclust_1 / A;
      R_numdensity_3 = numclust_3 / A;
   else
      N_singlet_fraction = singlet_fraction;
      N_equiv_radii = equiv_radii;
      N_density = density;
      N_n_pts = n_pts;
      N_numdensity_1 = numclust_1 / A;
      N_numdensity_3 = numclust_3 / A;
   end
end

% singles/total clust
SRplot({sextresults.singles_per_total_clusters, N_singlet_fraction, R_singlet_fraction}, 'Fractional amount of singles to all (%)')

% equivalent cluster radii
SRplot({sextresults.radius_per_2orMore_cluster, N_equiv_radii, R_equiv_radii}, 'Equivalent cluster radii (nm)')

% density of exposure/ROI
SRplot({sextresults.numclust_all, N_density, R_density}, 'Density of exposures per ROI (nm^{-2})')

% # of events per cluster
SRplot({sextresults.numobjs_per_multiple_cluster, N_n_pts, R_n_pts}, '# of localizations per multi-exposure')
 
% # of multiples/ROI
SRplot({sextresults.numdensity_3, N_numdensity_3, R_numdensity_3}, 'Density of multi-exposures per ROI (nm^{-2})')

% # of singles/ROI
SRplot({sextresults.numdensity_1, N_numdensity_1, R_numdensity_1}, 'Density of singles per ROI (nm^{-2})')
end

function SRplot(datums, y_text) 

   n = numel(datums);
   data   = [];
   groups = [];
   for i = 1 : n
      data = [data, datums{i}];
      groups = [groups, ones(size(datums{i})) + i - 1];
   end

   figure();
   G2diam = axes('LineWidth', 4, 'FontSize', 24);
   
   boxplot(data, groups, 'colors', [0 0 0], 'symbol', 'k+');
   ylabel(y_text);
   set(gca, 'xticklabel', {''}, 'box', 'off');
   set(gcf, 'color', [1 1 1]);

end
