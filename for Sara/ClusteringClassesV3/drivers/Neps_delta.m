clear all
close all

for i = 1 : 2
    
load('n_multi_Exper.mat');
densities_mean(isnan(densities_mean)) = 0;
densities_median(isnan(densities_median)) = 0;
E_E_range = E_range;
E_N_range = N_range;
if i == 1
   E_mean    = n_multi_mean;
   E_median  = n_multi_median;
elseif i == 2
   E_mean    = frac_clust_mean;
   E_median  = frac_clust_median;
elseif i == 3
   E_mean    = densities_mean;
   E_median  = densities_median;
end

load('n_multi_SimN.mat');
N_densities_mean(isnan(N_densities_mean)) = 0;
N_densities_median(isnan(N_densities_median)) = 0;
N_E_range = E_range;
N_N_range = N_range;
if i == 1
   N_mean    = N_numclust_3_mean;
   N_median  = N_numclust_3_median;
elseif i == 2
   N_mean    = N_frac_clust_mean;
   N_median  = N_frac_clust_median;
elseif i == 3
   N_mean    = N_densities_mean;
   N_median  = N_densities_median;
end

load('n_multi_SimR.mat');
R_densities_mean(isnan(R_densities_mean)) = 0;
R_densities_median(isnan(R_densities_median)) = 0;
R_E_range = E_range;
R_N_range = N_range;
if i == 1
   R_mean    = R_numclust_3_mean;
   R_median  = R_numclust_3_median;
elseif i == 2
   R_mean    = R_frac_clust_mean;
   R_median  = R_frac_clust_median;
elseif i == 3
   R_mean    = N_densities_mean;
   R_median  = N_densities_median;
end

if any(E_E_range ~= N_E_range) || any(E_E_range ~= R_E_range)
   E_E_range
   N_E_range
   R_E_range
   error('Inconsistencies in E_range!');
end

if any(E_N_range ~= N_N_range) || any(E_N_range ~= R_N_range)
   E_N_range
   N_N_range
   R_N_range
   error('Inconsistencies in R_range!');
end

dEN_mean   = + (E_mean   - N_mean);
dEN_median = + (E_median - N_median);
dER_mean   = + (E_mean   - R_mean);
dER_median = + (E_median - R_median);
%dEN_mean   = E_mean   ./ N_mean;
%dEN_median = E_median ./ N_median;
%dER_mean   = E_mean   ./ R_mean;
%dER_median = E_median ./ R_median;
dEN_mean(isnan(dEN_mean))     = 1;
dEN_median(isnan(dEN_median)) = 1;
dER_mean(isnan(dER_mean))     = 1;
dER_median(isnan(dER_median)) = 1;

if i == 1
   z_label = '# multi-clusters';
elseif i == 2
   z_label = 'fraction of localizations in clusters';
elseif i == 3
   z_label = 'multi-cluster average density';
end

% figure();
% hold on
% surf(E_range, N_range, dEN_mean');
% title('EN mean per ROI');
% xlabel('\epsilon (nm)');
% ylabel('minPts');
% zlabel(z_label);
% hold off
% 
% figure();
% hold on
% surf(E_range, N_range, dEN_median');
% title('EN median per ROI');
% xlabel('\epsilon (nm)');
% ylabel('minPts');
% zlabel(z_label);
% hold off

figure();
hold on
surf(E_range, N_range, dER_mean');
title('ER mean per ROI');
xlabel('\epsilon (nm)');
ylabel('minPts');
zlabel(z_label);
hold off

% figure();
% hold on
% surf(E_range, N_range, dER_median');
% title('ER median per ROI');
% xlabel('\epsilon (nm)');
% ylabel('minPts');
% zlabel(z_label);
% hold off

if i == 1
   fprintf('# multi-clusters\n\n');
elseif i == 2
   fprintf('fraction of localizations in clusters\n\n');
elseif i == 3
   fprintf('multi-cluster average density\n\n');
end
[i, j] = find(max(max(dER_mean)) == dER_mean);
fprintf('dER_mean   max at E = %2d, minPts = %2d [%2d, %2d]\n', ...
        E_range(i), N_range(j), i, j);
[i, j] = find(max(max(dER_median)) == dER_median);
fprintf('dER_median max at E = %2d, minPts = %2d [%2d, %2d]\n', ...
        E_range(i), N_range(j), i, j);
fprintf('\n');

end
