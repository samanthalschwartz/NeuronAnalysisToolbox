function [results, results_o, results_c, results_t, results_simulated] = ...
   testSingleLabel()

clear all
close all

SD = SimulateDomains();

SD.Printing = false;   % print statistics

SD.X_nm = 10000;       % x-dimension of domain (nm)
SD.Y_nm = 10000;       % y-dimension of domain (nm)

SD.Rho_d = 1.0e-07;    % domain density (1 / nm^2)
SD.Sigma_dom = 500;    % 2D Gaussian sigma for domain size (nm)
%SD.MAXdomains = 5;     % maximum number of domains allowed
SD.N_dom_part = 20;    % particles per domain
SD.Domain_sep = 100;   % domain center separation minimum (nm)
SD.N_observ = 10;      % observations per molecule
%SD.N_observ = 0;       % observations per molecule [SIMPLE MODEL -> 1]
SD.Sigma_loc = 20;     % localization error in each dimension (nm)

% polygonal ROI boundary
xy_region = [[0 1 1 2 2 3 3 4 4 5 5 4 4 3 3 2 2 1 1 0 0] * 2000; ...
             [0 0 1 1 0 0 1 1 0 0 2 2 3 3 4 4 3 3 2 2 0] * 2500]';

if exist('xy_region', 'var')
   SD.Rho_d = 5.0e-07;
   SD.Sigma_dom_xy = [500, 50];   % 2D (x, y) Gaussian sigma (nm) for producing
                                  % elliptical domains if nonzero
   SD.Fract_elongated = 0.7;      % fraction of structures that will be
                                  % elongated if Sigma_dom_xy != [0, 0]
   SD.Rot_angle = [0, pi];        % rotation range to be applied to domains 
   [N_domains, domain, domain_center, pts_center, sigmas, N_observations] = ...
      SD.generateDomains(xy_region);
   fprintf('N_domains = %d\n', N_domains);
   SD.plotDomains(N_domains, domain, domain_center, xy_region);
else   % rectangular ROI
   [N_domains, domain, domain_center, pts_center, sigmas, N_observations] = ...
      SD.generateDomains();
   fprintf('N_domains = %d\n', N_domains);
   SD.plotDomains(N_domains, domain, domain_center);
end

% concatenate results together into a single collection of points
XY = vertcat(domain{:});
XY_STD = vertcat(sigmas{:});
fprintf('N_points = %d\n', size(XY, 1));
Sigma_Reg = [10, 10];
pt_center = vertcat(pts_center{:});

SRc = SRcluster();
[xy_SR, sigma_SR, combined] = SRc.clusterSR(XY, XY_STD, Sigma_Reg);
SRc.E = SD.Sigma_loc;
if exist('xy_region', 'var')
   SRcollapseFig = SRc.plotSRcollapse(xy_region);
   [results_simulated, SRclusterFig] = ...
      SRc.analyzeSRclusters_simulated(pt_center, SD.X_nm, SD.Y_nm, xy_region);
else
   SRcollapseFig = SRc.plotSRcollapse();
   [results_simulated, SRclusterFig] = ...
      SRc.analyzeSRclusters_simulated(pt_center, SD.X_nm, SD.Y_nm);
end
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

algorithm = 'Hierarchical';
E = SRc.E;
minPts = 1;
c = Clustering();
% o = observations, c = collapsed, t = true
XY_o = XY;
XY_c = xy_SR;
XY_t = pt_center;
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
