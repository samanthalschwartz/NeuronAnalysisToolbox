close all; clear all;
filename = 'G:\zapERtrap\Raw Data\GLOBAL RELEASE\NL1\041718\TIFF files\cell1_stitch_AshleyFile';
load(filename)
%-- set imaging parameters:
baselineframe_start = 1; % first frame number that baseline acquisition begins
baselineframe_end = 7; % last frame number of baseline
baselineframerate = 1; % frame rate in minutes/frame 
releaseframe = 7; % time in minutes, after release, that first frame of post release starts
postrelease(1).frame_start = 8; % first frame number of post release
postrelease(1).frame_end = 'end'; % last frame number of post release - or 'end' if post release goes until end of series
postrelease(1).framerate = 2; % frame rate in minutes/frame

aa.imagingparams.releaseframe = 7;

aa.cleanedcargomask = [];
aa.cellFill.mask_img
aa.surfaceCargo.mask_img_highsens
aa.cellFill.selectSoma();
aa.cleanSurfaceCargoMask_Manual(); % call this line instead if you want to start again

savename = filename;
close all;
aa.cargo_heatmap = [];
h = aa.plotCargoHeatMap;
% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
save(fullfile([savename '-AshleyFile.mat']), 'aa','-v7.3'); 

saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
saveas(h,fullfile([savename '_timeHeatMap']),'eps');
close(h);
