function [results, results_o, results_c, results_t, results_simulated] = ...
   testSingleLabel()

clear all
close all

SD = SimulateDomains();

SD.Printing = false;   % Print statistics

SD.X_nm = 1000;       % x-dimension of domain (nm)
SD.Y_nm = 1000;       % y-dimension of domain (nm)

SD.Rho_d = 1.0e-05;    % domain density (1 / nm^2)
SD.Sigma_dom = 20;    % 2D Gaussian sigma for domain size (nm)
SD.Sigma_dom = 100;    % 2D Gaussian sigma for domain size (nm)
%SD.MAXdomains = 5;     % maximum number of domains allowed
%SD.N_dom_part = 4 ;    % particles per domain
SD.Domain_sep = 80;   % domain center separation minimum (nm)
SD.N_observ = 3;      % observations per molecule
SD.N_observ = 10;      % observations per molecule
%SD.N_observ = 0;       % observations per molecule [SIMPLE MODEL -> 1]
SD.Sigma_loc = 20;     % localization error in each dimension (nm)

SD.N_dom_part = 4;    % particles per domain
% poissrnd(0) = 0; %=> 1 particle always
% poissrnd(1)     %=> mostly 1 particle, but sometimes more



[N_domains, domain, domain_center, pts_center, sigmas, N_observations] = ...
   SD.generateDomains();
fprintf('N_domains = %d\n', N_domains);
SD.plotDomains(N_domains, domain, domain_center);
results_simulated.N_observations = N_observations;

XY = vertcat(domain{:});
fprintf('N_points = %d\n', size(XY, 1));
XY_STD = SD.Sigma_loc .* ones(size(XY));
Sigma_Reg = [10, 10];
Sigma_Reg=[0,0];
pt_center = vertcat(pts_center{:});

SRc = SRcluster();
%SRc.LoS = 0.99;
SRc.clusterSR(XY, XY_STD, Sigma_Reg);
SRc.E = 50;
[results_simulated, SRclusterFig] = ...
   SRc.analyzeSRclusters_simulated(pt_center, SD.X_nm, SD.Y_nm);
SRc.plotSRcollapse();
[results, analysisFigs] = SRc.analyzeSRclusters();
results_simulated.N_observations = N_observations;

figure();
hold on
hist(N_observations, 1 : max(N_observations));
xlabel('# of observations');
ylabel('# of localizations');
title(sprintf('mean # of observations per localization = %.3f +- %.3f\n', ...
              mean(N_observations), std(N_observations)));
hold off

algorithm = 'DBSCAN_Daszykowski';
E = SRc.E;
minPts = 3;
c = Clustering();
% o = observations, c = collapsed, t = true
XY_o = XY;
XY_c = pt_center;
XY_t = results.centers';
[nC_o, C_o, center_o, ptsI_o] = c.cluster(algorithm, XY_o, E, minPts);
[nC_c, C_c, center_c, ptsI_c] = c.cluster(algorithm, XY_c, E, minPts);
[nC_t, C_t, center_t, ptsI_t] = c.cluster(algorithm, XY_t, E, minPts);

results_o = c.clusterStats(XY_o, C_o, center_o);
results_c = c.clusterStats(XY_c, C_c, center_c);
results_t = c.clusterStats(XY_t, C_t, center_t);

fprintf('\n');
fprintf('o = observations, c = collapsed, t = true\n');
fprintf('-----------------------------------------\n');
fprintf('(o, c, t) nC = %d %d %d\n', results_o.nC, results_c.nC, results_t.nC);
fprintf('(o, c, t) n_points = %d %d %d\n', ...
   results_o.n_points, results_c.n_points, results_t.n_points);
fprintf('(o, c, t) n_clustered = %d %d %d\n', ...
   results_o.n_clustered, results_c.n_clustered, results_t.n_clustered);
fprintf('(o, c, t) n_isolated = %d %d %d\n', ...
   results_o.n_isolated, results_c.n_isolated, results_t.n_isolated);
fprintf('(o, c, t) min_c2c_dist = %.3f %.3f %.3f\n', ...
   results_o.min_c2c_dist, results_c.min_c2c_dist, results_t.min_c2c_dist);
fprintf('(o, c, t) min_e2e_dist = %.3f %.3f %.3f\n', ...
   results_o.min_e2e_dist, results_c.min_e2e_dist, results_t.min_e2e_dist);
fprintf('(o, c, t) mean(min_c2c_dists) = %.3f %.3f %.3f\n', ...
   mean(results_o.min_c2c_dists), mean(results_c.min_c2c_dists), ...
   mean(results_t.min_c2c_dists));
fprintf('(o, c, t) mean(min_e2e_dists) = %.3f %.3f %.3f\n', ...
   mean(results_o.min_e2e_dists), mean(results_c.min_e2e_dists), ...
   mean(results_t.min_e2e_dists));
fprintf('(o, c, t) mean(sigma_actual) = %.3f %.3f %.3f\n', ...
   mean(results_o.sigma_actual), mean(results_c.sigma_actual), ...
   mean(results_t.sigma_actual));
fprintf('(o, c, t) mean(areas) = %.3f %.3f %.3f\n', ...
   mean(results_o.areas), mean(results_c.areas), mean(results_t.areas));
fprintf('(o, c, t) mean(equiv_radii) = %.3f %.3f %.3f\n', ...
   mean(results_o.equiv_radii), mean(results_c.equiv_radii), ...
   mean(results_t.equiv_radii));
fprintf('(o, c, t) mean(n_pts_per_area) = %.3f %.3f %.3f\n', ...
   mean(results_o.n_pts_per_area), mean(results_c.n_pts_per_area), ...
   mean(results_t.n_pts_per_area));
fprintf('\n');
fprintf('[min, max](nearest_true_dists) = [%.3f, %.3f]\n', ...
   min(results_simulated.nearest_true_dists), ...
   max(results_simulated.nearest_true_dists));
fprintf('mean/median(nearest_true_dists) = %.3f, %.3f\n', ...
   mean(results_simulated.nearest_true_dists), ...
   median(results_simulated.nearest_true_dists));
