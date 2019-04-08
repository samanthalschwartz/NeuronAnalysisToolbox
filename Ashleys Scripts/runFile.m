close all; clear all;
%% define file paths


%% set up the Ashley File
% file parts, shifting, trimming image etc goes here






%% Masking
% mask cell -- change these to modify image smoothing
aa.cellFill.lsig = [1 1 1];
aa.cellFill.gsig = [1 1 1];
aa.cellFill.mask_img;
aa.cellFill.viewMaskOverlayPerim;
% mask surface Cargo -- change these to modify image smoothing
aa.surfaceCargo.lsig = [1 1 0];
aa.surfaceCargo.gsig = [1 1 1];
aa.surfaceCargo.mask_img_highsens
aa.surfaceCargo.viewMaskOverlayPerim;
% select the soma and save Ashley File
aa.cellFill.selectSoma();
save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
h = aa.plot_cargo_minFrame();
saveas(h,fullfile(datafilepath,[savename '_timeHeatMap']),'fig');
saveas(h,fullfile(datafilepath,[savename '_timeHeatMap']),'png');
close(h);
%% calculate number of objects
aa.calculateSurfaceCargoDistances(plotflag,plotsavedir);
%% calculate object density
aa.plotDensityperTime(distances);