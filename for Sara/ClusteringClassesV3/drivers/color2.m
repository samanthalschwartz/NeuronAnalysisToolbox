clear all
close all

SD = SimulateDomains();

SD.Printing = false;   % Print statistics

SD.X_nm = 10000;       % x-dimension of domain (nm)
SD.Y_nm = 10000;       % y-dimension of domain (nm)

SD.Rho_d = 1.0e-07;    % domain density (1 / nm^2)
SD.Sigma_dom = 500;    % 2D Gaussian sigma for domain size (nm)
SD.N_dom_part = 32;    % particles per domain
SD.Domain_sep = 202;   % domain center separation minimum (nm)
SD.N_observ = 10;      % observations per molecule
%SD.N_observ = 0;       % observations per molecule [SIMPLE MODEL -> 1]
SD.Sigma_loc = 20;     % localization error in each dimension (nm)

[N_domains1, domain1, d_center1, p_center1] = SD.generateDomains();
fprintf('N_domains = %d\n', N_domains1);
SD.plotDomains(N_domains1, domain1, d_center1);

SD.N_dom_part = 64;    % particles per domain
SD.Domain_sep = 404;   % domain center separation minimum (nm)

[N_domains2, domain2, d_center2, p_center2] = SD.generateDomains();
fprintf('N_domains = %d\n', N_domains2);
SD.plotDomains(N_domains2, domain2, d_center2);

figure
hold on
SD.plotDomainsOneColor(N_domains1, domain1, d_center1, 'k');
SD.plotDomainsOneColor(N_domains2, domain2, d_center2, 'r');
hold off
