function [results, results_o, results_c, results_t, results_simulated] = ...
   Nalpha_sim1(A, minPts, Quiet, Random, Collapse)

%clear all
%close all

if exist('Quiet', 'var')
   quiet = Quiet;
else
   quiet = false;
end

if exist('Random', 'var')
   RANDOM = Random;
else
   RANDOM = false;
end
if exist('Collapse', 'var')
   COLLAPSE = Collapse;
else
   COLLAPSE = false;
end
if ~quiet
   fprintf('RANDOM = %d, COLLAPSE = %d\n\n', RANDOM, COLLAPSE);
end

SD = SimulateDomains();
SD.Printing = false;   % Print statistics
SD.X_nm = 1000;        % x-dimension of domain (nm)
SD.Y_nm = 1000;        % y-dimension of domain (nm)

%SD.MAXdomains = 5;     % maximum number of domains allowed
% poissrnd(0) = 0; % => 1 particle always
% poissrnd(1)      % => mostly 1 particle, but sometimes more

%n_pts = 2391;   % observations
n_pts = 671;   % localizations
if RANDOM
   %SD.Rho_d = 1.0e-03;   % domain density (1 / nm^2)
   %SD.Sigma_dom = 0;     % 2D Gaussian sigma for domain size (nm)
   %SD.Domain_sep = 0;    % domain center separation minimum (nm)
   %SD.N_dom_part = 0 ;   % particles per domain
   %SD.N_observ = 0;      % observations per molecule [SIMPLE MODEL -> 1]
   %SD.Sigma_loc = 20;    % localization error in each dimension (nm)

   [XY, XY_STD] = SD.generateRandoms(n_pts);
   N_observations = size(XY, 1);
   fprintf('RANDOM N_points = %d\n', N_observations);
   Sigma_Reg = [0, 0];
   pt_center = XY;
else
   SD.Rho_d = 3.7e-05;      % domain density (# / nm^2)
%  SD.Sigma_dom = 23.14;    % 2D Gaussian sigma for domain size (nm)
%  SD.Sigma_dom = 11.57;    % 2D Gaussian sigma for domain size (nm) = r / 2
   SD.Sigma_dom = 7.71;     % 2D Gaussian sigma for domain size (nm) = r / 3
%  SD.Domain_sep = 38.53;   % domain center separation minimum (nm)
   SD.Domain_sep = 19.26;   % domain center separation minimum (nm) [PLAY]
   SD.N_dom_part = 11.58;   % particles per domain
%  SD.N_observ = 2;         % observations per molecule
   SD.N_observ = 0;         % observations per molecule [SIMPLE MODEL -> 1]
   SD.Sigma_loc = 20;       % localization error in each dimension (nm)

   singlet_fract = 0.86;

   [N_domains, domain, domain_center, pts_center, sigmas, N_observations] = ...
      SD.generateDomains();
   fprintf('N_domains = %d\n', N_domains);

   XY = vertcat(domain{:});
   N_points = size(XY, 1);
   fprintf('N_points = %d\n', N_points);
   XY_STD = SD.Sigma_loc .* ones(size(XY));
   Sigma_Reg = [0, 0];
   pt_center = vertcat(pts_center{:});

   area_ROI = 1000^2;
   area_d   = N_domains * pi * (3 * SD.Sigma_dom)^2;
   f_d_ROI  = area_d / area_ROI;
   fudge = 1 / (1 - f_d_ROI);
%  N_singlets = n_pts - SD.Rho_d * area_ROI * SD.N_dom_part * fudge;
   % Compute the number of singlets to be produced from the specified singlet
   % fraction and the number of non-singlet domains generated.  Note that
   %    singlet_fraction = N_singlets / (N_domains + N_singlets)
   % However, some of the random singlets will fall on top of previously
   % generated domains (on average, f_d_ROI = area_d / area_ROI), so increase
   % the number of singlets generated to account for this loss.
   N_singlets = round(singlet_fract / (1 - singlet_fract) * N_domains * fudge);
   fprintf('N_singlets = %d\n', N_singlets);
   [XYs, XYs_STD] = SD.generateRandoms(N_singlets);
   XY = [XY; XYs];
   XY_STD = [XY_STD; XYs_STD];
   pt_center = [pt_center; XY];

   N_total = size(XY, 1);
   fprintf('N_total = %d\n', N_total);

%  SD.plotDomains(N_domains, domain, domain_center);
%  figure();
%  hold on
%  plot(XY(:, 1), XY(:, 2), 'k.');
%  title(sprintf('N_{total} = %d', N_total));
%  xlabel('x (nm)');
%  ylabel('y (nm)');
%  hold off
end

results = [];
results_simulated = [];
if COLLAPSE
   SRc = SRcluster();
   [xy_SR, sigma_SR, combined] = SRc.clusterSR(XY, XY_STD, Sigma_Reg);
   SRc.Alpha = A;
   SRc.E = 0;
   %SRc.PlotFigures = false; % !!!
   [results_simulated, SRclusterFig] = ...
      SRc.analyzeSRclusters_simulated(pt_center, SD.X_nm, SD.Y_nm);
   SRc.plotSRcollapse();
   SRc.plotSRclusters();
   SRc.PlotFigures = false;
   SRc.Printing = false;
   [results, analysisFigs] = SRc.analyzeSRclusters();
   results_simulated.N_observations = N_observations;

   fprintf('nC = %d\n', results.numclust(3));

   %figure();
   %hold on
   %hist(N_observations, 1 : max(N_observations));
   %xlabel('# of observations');
   %ylabel('# of localizations');
   %title(...
   %   sprintf('mean # of observations per localization = %.3f +- %.3f\n', ...
   %           mean(N_observations), std(N_observations)));
   %hold off
end

algorithm = 'Voronoi';
c = Clustering();
c.Valgorithm = 2;
c.Alpha = A;
E = 0;
% o = observations, c = collapsed, t = true
XY_o = XY;
if COLLAPSE
   %XY_c = results.centers';
   XY_c = xy_SR;
else
   XY_c = XY_o;
end
XY_t = pt_center;
[nC_o, C_o, center_o, ptsI_o] = c.cluster(algorithm, XY_o, E, minPts);
[nC_c, C_c, center_c, ptsI_c] = c.cluster(algorithm, XY_c, E, minPts);
[nC_t, C_t, center_t, ptsI_t] = c.cluster(algorithm, XY_t, E, minPts);

results_o = c.clusterStats(XY_o, C_o, center_o);
results_c = c.clusterStats(XY_c, C_c, center_c);
results_t = c.clusterStats(XY_t, C_t, center_t);

%clusterFig_o = c.plotClusters(XY_o, C_o, center_o, ptsI_o, algorithm);
%showm(clusterFig_o);
%clusterFig_c = c.plotClusters(XY_c, C_c, center_c, ptsI_c, algorithm);
%showm(clusterFig_c);
%clusterFig_t = c.plotClusters(XY_t, C_t, center_t, ptsI_t, algorithm);
%showm(clusterFig_t);

if ~quiet
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
fprintf('(o, c, t) singlet_fraction = %.3f %.3f %.3f\n', ...
   results_o.singlet_fraction, results_c.singlet_fraction, ...
   results_t.singlet_fraction);
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
if COLLAPSE
   fprintf('[min, max](nearest_true_dists) = [%.3f, %.3f]\n', ...
      min(results_simulated.nearest_true_dists), ...
      max(results_simulated.nearest_true_dists));
   fprintf('mean/median(nearest_true_dists) = %.3f, %.3f\n', ...
      mean(results_simulated.nearest_true_dists), ...
      median(results_simulated.nearest_true_dists));
end

end
