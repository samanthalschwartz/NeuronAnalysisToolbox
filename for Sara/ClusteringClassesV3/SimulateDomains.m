classdef SimulateDomains < handle

% SimulateDomains class written by Michael Wester and Keith Lidke (2/27/2018)
%    <wester@math.unm.edu>
% The New Mexico Center for the Spatiotemporal Modeling of Cell Signaling
% University of New Mexico Health Sciences Center
% Albuquerque, New Mexico, USA   87131
% Copyright (c) 2015-2018 by Michael J. Wester and Keith A. Lidke
%
% Example main program:
%
%    SD = SimulateDomains();
%
%    SD.Printing = false;   % Print statistics
%
%    SD.X_nm = 10000;       % x-dimension of domain (nm)
%    SD.Y_nm = 10000;       % y-dimension of domain (nm)
%
%    SD.Rho_d = 1.0e-07;    % domain density (1 / nm^2) [KEITH]
%    SD.Sigma_dom = 500;    % 2D Gaussian sigma for domain size (nm) [KEITH]
%    SD.N_dom_part = 32;    % particles per domain
%    SD.Domain_sep = 202;   % domain center separation minimum (nm)
%    SD.N_observ = 10;      % observations per molecule [KEITH]
%    %SD.N_observ = 0;       % observations per molecule [SIMPLE MODEL -> 1]
%    SD.Sigma_loc = 20;     % localization error in each dimension (nm) [KEITH]
%
%    [N_domains1, domain1, center1] = SD.generateDomains();
%    fprintf('N_domains = %d\n', N_domains1);
%    SD.plotDomains(N_domains1, domain1, center1);
%
%    SD.N_dom_part = 64;    % particles per domain
%    SD.Domain_sep = 404;   % domain center separation minimum (nm)
%
%    [N_domains2, domain2, center2] = SD.generateDomains();
%    fprintf('N_domains = %d\n', N_domains2);
%    SD.plotDomains(N_domains2, domain2, center2);
%
%    figure();
%    hold on
%    SD.plotDomainsOneColor(N_domains1, domain1, center1, 'k');
%    SD.plotDomainsOneColor(N_domains2, domain2, center2, 'r');
%    hold off
%
% See testSingleLabel.m:
%
%    SD.N_dom_part = 20;    % particles per domain
%    SD.Domain_sep = 100;   % domain center separation minimum (nm)
%
%    % polygonal ROI boundary 
%    xy_region = [[0 1 1 2 2 3 3 4 4 5 5 4 4 3 3 2 2 1 1 0 0] * 2000; ...
%                [0 0 1 1 0 0 1 1 0 0 2 2 3 3 4 4 3 3 2 2 0] * 2500]';
%
%   if exist('xy_region', 'var')
%      SD.Rho_d = 5.0e-07;
%      SD.Sigma_dom_xy = [500, 50];   % 2D (x, y) Gaussian sigma (nm) for
                                      % producing elliptical domains if nonzero
%      SD.Fract_elongated = 0.7;      % fraction of structures that will be
%                                     % elongated if Sigma_dom_xy != [0, 0]
%      SD.Rot_angle = [0, pi];        % rotation range to be applied to domains 
%      [N_domains, domain, domain_center, pts_center, sigmas, ...
%       N_observations] = SD.generateDomains(xy_region);
%      fprintf('N_domains = %d\n', N_domains);
%      SD.plotDomains(N_domains, domain, domain_center, xy_region);
%   else   % rectangular ROI
%      [N_domains, domain, domain_center, pts_center, sigmas, ...
%       N_observations] = SD.generateDomains();
%      fprintf('N_domains = %d\n', N_domains);
%      SD.plotDomains(N_domains, domain, domain_center);
%   end
%
% See D3.m
%
%    SD.Dim = 3;
%    SD.Z_nm = 10000;
%    SD.Rho_d = 1.0e-11;    % domain density (1 / nm^2)
%    %SD.Sigma_dom_xy = [500, 50, 5000];
%    [N_domains, domain, domain_center, pts_center, sigmas, ...
%     N_observations] = SD.generateDomains();
%    fprintf('N_domains = %d\n', N_domains);
%    SD.plotDomains3(N_domains, domain, domain_center);

% =============================================================================
properties
% =============================================================================

   Printing = false;    % Print statistics

   X_nm = 10000;        % x-dimension of domain (nm)
   Y_nm = 10000;        % y-dimension of domain (nm)
   Z_nm = -1;           % z-dimension of domain (nm)

   Rho_d = 10^(-7);     % domain density (1 / nm^2) [KEITH]
   Sigma_dom = 500;     % 2D Gaussian sigma for domain size (nm) [KEITH]
   Sigma_dom_xy = [0,0];% 2D (x, y) Gaussian sigma (nm) for producing
                        % elliptical domains if nonzero
   Rot_angle = [0, 0];  % rotation range to be applied to domains; useful with
                        % Sigma_dom_xy to produce linear domains oriented
                        % over a range of rotations like scattered needles
   Fract_elongated = 1; % fraction of structures that will be elongated if
                        % Sigma_dom_xy != [0, 0]
   MAXdomains = 100;    % maximum number of domains allowed

   N_dom_part = 25;     % mean particles per domain
   Domain_sep = 100;    % minimum domain center separation (nm)
   MAX_tries = 1000;    % maximum number of tries to keep domains separated

   N_observ = 10;       % observations per molecule [KEITH]
   %N_observ = 0;        % observations per molecule [SIMPLE MODEL -> 1]
   Sigma_loc = 20;      % localization error in each dimension (nm) [KEITH]

   Dim = 2;             % allowed values are 2 (2D) and 3 (3D)

% =============================================================================
end % properties

methods
% =============================================================================

function [N_domains, domain, domain_center, pts_center, sigmas, ...
          N_observations] = generateDomains(obj, xy_region)
% Generate domains.
%
% Input:
%   xy_region        OPTIONAL (x, y) coordinates of the domain vertices [N x 2]
% Output:
%   N_domains        number of domains produced
%   domain           cell array containing the points of each generated domain
%   domain_center    cell array of the coordinates of each original domain
%                    center
%   pts_center       cell array (per domain) of the centers of each observation
%                    (localization) cluster
%   sigmas           cell array containing the sigmas of the generated points
%   N_observations   array of the numbers of observations generated per
%                    localization

   % Gaussian distribution simulation, from SuperCluster, Carolyn Pehlke 2014
   % Create simuated data based on input parameters (originally written by
   % Keith Lidke).

   % Compute properties of the domain
   if exist('xy_region', 'var')
      [xy, minimum, delta, area] = obj.domainProperties(xy_region);
   else
      [xy, minimum, delta, area] = obj.domainProperties();
   end

   sd_nonzero = all(obj.Sigma_dom_xy ~= zeros(1, numel(obj.Sigma_dom_xy)));

   % === Choose the number of molecular domains via a Poisson distribution ===
   n_tries = 0;
   N_domains = 0;
   while N_domains == 0 && n_tries < obj.MAX_tries
      n_tries = n_tries + 1;
      N_domains = min(poissrnd(obj.Rho_d * area), obj.MAXdomains);
   end

   % Error if the number of domains is zero!
   if N_domains == 0
      error('N_domains = 0!');
   elseif obj.Printing
      fprintf('N_domains = %d (initial)\n', N_domains);
   end

   domain = cell(N_domains, 1); % create cell array for each domain
   N_observations = [];
   Nd = 0;
   while Nd < N_domains
      n_tries = 0;
      while n_tries < obj.MAX_tries
         n_tries = n_tries + 1;

         % Domain center
         X_c = delta .* rand(1, obj.Dim) + minimum;
         % If the new "center" is not within the domain, try again
         if ~obj.inside(X_c, xy)
            continue;
         end

         % Check that the new domain's center is well separated from previous
         % domain centers
         separated = true;
         i = 0;
         while i < Nd && separated
            i = i + 1;
            separated = pdist([X_c; domain_center{i}]) >= obj.Domain_sep;
         end
         if separated
            % Domain is well separated from previous domains, so done
            break;
         end
         % Domain not well separated from previous domains, so try again
      end

      if n_tries >= obj.MAX_tries && ~separated
         % Give up trying to add a new domain, so just exit the loop ...
         break;
      end

      % Domain center accepted!
      Nd = Nd + 1;

      if sd_nonzero && rand <= obj.Fract_elongated
         rot_angle = ...
            rand * (obj.Rot_angle(2) - obj.Rot_angle(1)) + obj.Rot_angle(1);
         [N_observations_D, pts_center_D, tmpclust, tmpsigma] =      ...
            obj.generateLocalizations(Nd, xy, X_c, obj.Sigma_dom_xy, ...
                                      rot_angle);
      else
         [N_observations_D, pts_center_D, tmpclust, tmpsigma] = ...
            obj.generateLocalizations(Nd, xy, X_c,              ...
                                      repmat(obj.Sigma_dom, [1, obj.Dim]), 0);
      end
      N_observations = [N_observations, N_observations_D];
      pts_center{Nd} = pts_center_D;

      % Save results
      domain_center{Nd} = X_c;
      domain{Nd} = tmpclust;
      sigmas{Nd} = tmpsigma;
   end

   % Reset N_domains as the actual number of domains may be less than was
   % originally asked for due to separation constraints
   N_domains = Nd;
   if obj.Printing
      fprintf('N_domains = %d (applying separation constraints)\n', N_domains);
   end
   domain_center = domain_center';

end

% -----------------------------------------------------------------------------

function [N_observations, pts_center, tmpclust, tmpsigma] = ...
   generateLocalizations(obj, Nd, xy, X_c, sigma_dom_xy, rot_angle)
% Generate localizations.
%
% Input:
%    Nd               domain number
%    xy               (x, y) coordinates of the ROI
%    X_c              domain center
%    sigma_dom_xy     (x, y) components of the domain sigma
%    rot_angle        rotation angle (radians, counterclockwise)
% Output:
%    N_observations   number of observations per localization
%    pts_center       center of each cluster of observations (i.e., the true
%                     localization)
%    tmpclust         observations generated around the localizations
%    tmpsigma         sigmas of the generated points

   N_observations = [];

   % Number of molecules: Poisson random numbers with N_dom_part = lambda
   N_mol = max(1, poissrnd(obj.N_dom_part));
   if obj.Printing
      fprintf('domain %d: %d (initial)\n', Nd, N_mol);
   end
   % Choose domain distortion
   pts_center = zeros(N_mol, obj.Dim);
   tmpclust = []; % temporary cluster
   tmpsigma = [];
   for nm = 1 : N_mol
      X_mol = repmat(-1e+10, [1, obj.Dim]);
      % Pick a point inside the sector
      while ~obj.inside(X_mol, xy)
         X_mol =                                                            ...
            SimulateDomains.rotateCoords(sigma_dom_xy .* randn(1, obj.Dim), ...
                                         rot_angle) + X_c;
      end
      pts_center(nm, :) = X_mol;
      [pts, sigma] = obj.generateObservations(X_mol);
      %pts = SimulateDomains.makeCluster(N_mol, X_c, sigma_dom_xy);
      N_obs = size(pts, 1);
      N_observations = [N_observations, N_obs];
      tmpclust = [tmpclust; pts];
      tmpsigma = [tmpsigma; sigma];
   end

   % Delete points that are outside the sector in either coordinate
   if obj.Dim == 2
      outside = ~arrayfun(@(xc, yc) obj.inside([xc, yc], xy), ...
                          tmpclust(:, 1), tmpclust(:, 2));
   else
      outside = ~arrayfun(@(xc, yc, zc) obj.inside([xc, yc, zc], xy), ...
                          tmpclust(:, 1), tmpclust(:, 2), tmpclust(:, 3));
   end
   tmpclust(outside, :) = [];
   tmpsigma(outside, :) = [];
   N_mol = size(tmpclust, 1);
   if obj.Printing
      fprintf('domain %d: %d (removed outside x,y)\n', Nd, N_mol);
   end

end

% -----------------------------------------------------------------------------

function [pts, sigmas] = generateObservations(obj, center)
% Generate observations about the given center point.
%
% Input:
%    center   center of the cluster of observations
% Output:
%    pts      observations generated around the center
%    sigmas   sigmas of the generated points

   % Number of observations (localizations)
   if obj.N_observ <= 0
      N_obs = 1;
   else
      N_obs = max(1, poissrnd(obj.N_observ));
   end

   sigmas = obj.generateSigmas(N_obs);
   pts = sigmas .* randn(N_obs, obj.Dim) + repmat(center, N_obs, 1);

end

% -----------------------------------------------------------------------------

function sigmas = generateSigmas(obj, N_obs)
% Generate N_obs(ervations) sigmas.
%
% Input:
%    N_obs    number of observations
% Output:
%    sigmas   sigmas of the generated points

   % sigma = 1.3 nm / sqrt(N) where N is chosen from an exponential
   % distribution of the number of photons with mean 500 photons.
   % Values of sigma > Sigma_loc nm are eliminated, so generate more values
   % than needed (2 * # points created) to be sure to have enough.
   ss = [];
   while length(ss) < obj.Dim*N_obs
      N_photons = exprnd(500 .* ones((obj.Dim + 1)*N_obs, 1));
      s = 1.3 * (16000/150) ./ sqrt(N_photons);
      s(s > obj.Sigma_loc) = [];
      ss = [ss; s];
   end
   sigmas = zeros(N_obs, obj.Dim);
   sigmas(:, 1) = ss(1 : N_obs);
   sigmas(:, 2) = ss(N_obs + 1 : 2*N_obs);
   if obj.Dim == 3
      sigmas(:, 3) = ss(2*N_obs + 1 : 3*N_obs);
   end

end

% -----------------------------------------------------------------------------

function [xy, minimum, delta, area] = domainProperties(obj, xy_region)
% Compute properties of the domain.
%
% Input:
%    xy_region   vertices of the domain supplied by the user (2D or 3D)
% Output:
%    xy          (x, y {, z}) coordinates defining the region
%    minimum     minimum (x, y {, z})
%    delta       maximum (x, y {, z}) - minimum (x, y {, z})
%    area        area (2D) or volume (3D) of the region

   if exist('xy_region', 'var');
      % Vertices of the domain supplied by the user
      x = xy_region(:, 1);
      y = xy_region(:, 2);
      xy = [x, y];
      if obj.Dim == 3
         z = xy_region(:, 3);
         xy = [x, y, z];
      end
   else
      if obj.Dim == 2
         % Default rectangular domain
         x = [0, obj.X_nm, obj.X_nm, 0, 0]';
         y = [0, 0, obj.Y_nm, obj.Y_nm, 0]';
         xy = [x, y];
      else
         if obj.Z_nm <= 0
            error('Z_nm not positive!');
         else
            % Default cuboid domain
            x = [0, obj.X_nm, obj.X_nm, 0, 0, obj.X_nm, obj.X_nm, 0, 0]';
            y = [0, 0, obj.Y_nm, obj.Y_nm, 0, 0, obj.Y_nm, obj.Y_nm, 0]';
            z = [0, 0, 0, 0, obj.Z_nm, obj.Z_nm, obj.Z_nm, obj.Z_nm, 0]';
            xy = [x, y, z];
         end
      end
   end
   min_x = min(x);
   max_x = max(x);
   min_y = min(y);
   max_y = max(y);
   delta_x = max_x - min_x;
   delta_y = max_y - min_y;
   if obj.Dim == 2
      minimum = [min_x, min_y];
      delta = [delta_x, delta_y];
      area = polyarea(x, y);
   else
      min_z = min(z);
      max_z = max(z);
      delta_z = max_z - min_z;
      minimum = [min_x, min_y, min_z];
      delta = [delta_x, delta_y, delta_z];
      area = delta_x * delta_y * delta_z;
   end

end

% =============================================================================

function [pts, sigmas] = generateRandoms(obj, n_pts, xy_region)
% Generate a set of approximately n_pts random points, where the number is
% chosen from a Poisson distribution centered on n_pts.
%
% Input:
%    n_pts       (approximate) number of point to generate
%    xy_region   vertices of the domain supplied by the user (2D or 3D) in
%                which the generated points will lie
% Output:
%    pts         coordinates of the generated points
%    sigmas      generated uncertainties (sigmas) associated with each point

   if exist('xy_region', 'var')
      [xy, minimum, delta, area] = obj.domainProperties(xy_region);
   else
      [xy, minimum, delta, area] = obj.domainProperties();
   end
   center = minimum + delta / 2;

   % Number of points
   if n_pts <= 0
      N_pts = 1;
   else
      N_pts = max(1, poissrnd(n_pts));
   end

   % Generate more values than needed as many likely will lay outside the
   % domain boundary
   N_pts = 10 * N_pts;

   % Generate points and sigmas
   sigmas = obj.generateSigmas(N_pts);
   pts = (2 * rand(N_pts, obj.Dim) - 1) * diag(delta) ...
         + repmat(center, N_pts, 1);

   % Delete points that are outside the sector in either coordinate
   if obj.Dim == 2
      outside = ~arrayfun(@(xc, yc) obj.inside([xc, yc], xy), ...
                          pts(:, 1), pts(:, 2));
   else
      outside = ~arrayfun(@(xc, yc, zc) obj.inside([xc, yc, zc], xy), ...
                          pts(:, 1), pts(:, 2), pts(:, 3));
   end
   pts(outside, :)    = [];
   sigmas(outside, :) = [];
   

   if obj.Printing
      fprintf('%d -> %d (removed outside x,y)\n', N_pts, size(pts, 1));
   end

   % Remove excess points
   N_pts = N_pts / 10;
   if size(pts, 1) > N_pts
      pts(N_pts + 1 : end, :)    = [];
      sigmas(N_pts + 1 : end, :) = [];
   end

end

% =============================================================================

function plotDomainsOneColor(obj, N_domains, domain, domain_center, color)
% 2D plot of the random domains.
%
% Input:
%    N_domains       number of domains
%    domain          domain structure containing cluster member coordinates
%    domain_center   domain center (x, y) coordinates [N_domains x 2]
%    color           line color to use

   for i = 1 : N_domains
      % Plot the points contained in each random domain
      xy = domain{i};
      plot(xy(:, 1), xy(:, 2), [color, '.'], 'MarkerSize', 10);
      % Plot the center of each random domain as a large circle
      c = domain_center{i};
      plot(c(:, 1), c(:, 2), [color, 'o'], 'MarkerSize', 25);
   end
   axis([0, obj.X_nm, 0, obj.Y_nm]);
   axis equal;
   xlabel('x (nm)');
   ylabel('y (nm)');

end

% -----------------------------------------------------------------------------

function plotDomainsOneColor3(obj, N_domains, domain, domain_center, color)
% 3D plot of the random domains.
%
% Input:
%    N_domains       number of domains
%    domain          domain structure containing cluster member coordinates
%    domain_center   domain center (x, y, z) coordinates [N_domains x 3]
%    color           line color to use

   for i = 1 : N_domains
      % Plot the points contained in each random domain
      xy = domain{i};
      plot3(xy(:, 1), xy(:, 2), xy(:, 3), [color, '.'], 'MarkerSize', 10);
      % Plot the center of each random domain as a large circle
      c = domain_center{i};
      plot3(c(:, 1), c(:, 2), c(:, 3), [color, 'o'], 'MarkerSize', 25);
   end
   axis([0, obj.X_nm, 0, obj.Y_nm, 0, obj.Z_nm]);
   axis equal;
   xlabel('x (nm)');
   ylabel('y (nm)');
   zlabel('z (nm)');

end

% -----------------------------------------------------------------------------

function plotDomains(obj, N_domains, domain, domain_center, xy_region)
% 2D plot of the random domains.
%
% Input:
%    N_domains       number of domains
%    domain          domain structure containing cluster member coordinates
%    domain_center   domain center (x, y) coordinates [N_domains x 2]
%    xy_region       [OPTIONAL] rectangular region defining the maximum extent
%                    of the domain

   color = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];
   n_colors = length(color);

   figure;
   hold on
   for i = 1 : N_domains
      % Plot the points contained in each random domain
      xy = domain{i};
      plot(xy(:, 1), xy(:, 2), ...
           [color(SimulateDomains.nMODm(i, n_colors)), '.'], 'MarkerSize', 10);
      % Plot the center of each random domain as a large circle
      c = domain_center{i};
      plot(c(:, 1), c(:, 2), ...
           [color(SimulateDomains.nMODm(i, n_colors)), 'o'], 'MarkerSize', 25);
   end
   if exist('xy_region', 'var')
      plot(xy_region(:, 1), xy_region(:, 2), 'k-', 'LineWidth', 3);
      axis([min([0; xy_region(:, 1)]), max([obj.X_nm; xy_region(:, 1)]), ...
            min([0; xy_region(:, 2)]), max([obj.Y_nm; xy_region(:, 2)])]);
   else
      axis([0, obj.X_nm, 0, obj.Y_nm]);
   end
   axis equal;
   title(sprintf('N_{domains} = %d', N_domains));
   xlabel('x (nm)');
   ylabel('y (nm)');
   hold off

end

% -----------------------------------------------------------------------------

function plotDomains3(obj, N_domains, domain, domain_center, xy_region)
% 3D plot of the random domains.
%
% Input:
%    N_domains       number of domains
%    domain          domain structure containing cluster member coordinates
%    domain_center   domain center (x, y, z) coordinates [N_domains x 3]
%    xy_region       [OPTIONAL] 3D rectangular box defining the maximum extent
%                    of the domain

   color = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];
   n_colors = length(color);

   figure;
   hold on
   for i = 1 : N_domains
      % Plot the points contained in each random domain
      xy = domain{i};
      plot3(xy(:, 1), xy(:, 2), xy(:, 3), ...
           [color(SimulateDomains.nMODm(i, n_colors)), '.'], 'MarkerSize', 10);
      % Plot the center of each random domain as a large circle
      c = domain_center{i};
      plot3(c(:, 1), c(:, 2), c(:, 3), ...
           [color(SimulateDomains.nMODm(i, n_colors)), 'o'], 'MarkerSize', 25);
   end
   if exist('xy_region', 'var')
      plot3(xy_region(:, 1), xy_region(:, 2), xy_region(:, 3), ...
            'k-', 'LineWidth', 3);
      axis([min([0; xy_region(:, 1)]), max([obj.X_nm; xy_region(:, 1)]), ...
            min([0; xy_region(:, 2)]), max([obj.Y_nm; xy_region(:, 2)]), ...
            min([0; xy_region(:, 3)]), max([obj.Z_nm; xy_region(:, 3)])]);
   else
      axis([0, obj.X_nm, 0, obj.Y_nm, 0, obj.Z_nm]);
   end
   axis equal;
   title(sprintf('N_{domains} = %d', N_domains));
   xlabel('x (nm)');
   ylabel('y (nm)');
   zlabel('z (nm)');
   hold off

end

% -----------------------------------------------------------------------------

function tf = inside(obj, X_c, xy)
% True if the point with (x, y) coordinates X_c is inside the polygonal region
% (rectangular box) defined by xy in 2D (3D), otherwise false.
%
% Input:
%    X_c   point to be tested [1 x 2]
%    xy    polygonal region or 3D rectangular box coordinates [N x 2]
% Output:
%    tf    true if inside, false if outside

   if obj.Dim == 2
      tf = inpolygon(X_c(1), X_c(2), xy(:, 1), xy(:, 2));
   else
      min_xy = min(xy);
      max_xy = max(xy);
      tf = all(arrayfun(@(i) all(min_xy(i) <= X_c(:, i) & ...
                                 X_c(:, i) <= max_xy(i)), 1:3));
   end

end

% =============================================================================
end % methods

methods(Static)
% =============================================================================

function r = nMODm(n, m)
% Modulus such that r is in [1, m] rather than [0, m - 1].
%
% Input:
%    n   integer
%    m   integer modulus
% Output:
%    r   remainder

   r = mod(n, m);
   if r <= 0
      r = r + m;
   end

end

% -----------------------------------------------------------------------------

function xy_r = rotateCoords(xy, rot_angle)
% Rotate the (x, y) coordinates xy by the rotation angle rot_angle (radians).
%
% Input:
%    xy          (x, y) coordinates [N x 2]
%    rot_angle   rotation angle (radians)
% Output:
%    xy_r        rotated (x, y) coordinates

   ca = cos(rot_angle);
   sa = sin(rot_angle);

   if numel(xy) == 2
      xy_r = ([ca, -sa; sa, ca] * xy')';
   else
      % 2D rotation
      xy_r = ([ca, -sa, 0; sa, ca, 0; 0, 0, 1] * xy')';
   end

end

% -----------------------------------------------------------------------------

function xy = makeCluster(n_points, domain_center, radius)
% Generate a series of normally distributed points
%
% Input:
%    n_points        number of points to be in the cluster
%    domain_center   location of the cluster's center [1 x 2]
%    radius          maximum distance of cluster points from the cluster center
% Output:
%    xy              (x, y) coordinates of cluster members [N x 2]

   xy = randn(n_points, 2) .* radius;

   % Remove points outside of [-r, r] in either coordinate
   %r = 2 * radius;
   %xy(find(xy(:, 1) < -r | xy(:, 1) > r | ...
   %        xy(:, 2) < -r | xy(:, 2) > r), :) = [];

   % Remove points separated by more than r from n_points - m other points
   % where m is a fraction of the original number of points
   r = radius;
   m = round(0.15 * n_points);
   D = sum(squareform(pdist(xy) > r));
   xy(find(D >= n_points - m), :) = [];

   % Translate to the specified center
   xy(:, 1) = xy(:, 1) + domain_center(1);
   xy(:, 2) = xy(:, 2) + domain_center(2);

end

% =============================================================================
end % methods(Static)
% =============================================================================
end % classdef
