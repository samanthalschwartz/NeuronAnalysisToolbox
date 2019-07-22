close all; clear all;
filename = 'Y:\Lab Projects\zapERtrap\Raw Data\LOCAL RELEASE\050918_local_NL1\cell2dendrite.tif';
uiopen(filename); close all;
%-- set imaging parameters:
aa = AshleyAnalysis();
aa.imagingparams.releaseframe = 12;
aa.cellFill.image = image(:,:,:,1);
aa.surfaceCargo.image = image(:,:,:,2);
close all;

aa.cellFill.mask_img
aa.surfaceCargo.mask_img_highsens
aa.cellFill.selectSoma();

aa.c
aa.cleanSurfaceCargoMask_Manual(1)
close all;
aa.cargo_heatmap = [];
h = aa.plotCargoHeatMap;
% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
savename = filename(1:end-4);
save(fullfile([savename '-AshleyFile.mat']), 'aa','-v7.3'); 

saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
saveas(h,fullfile([savename '_timeHeatMap']),'eps');
close(h);
%%
close all; clear all;
filename = 'Y:\Lab Projects\zapERtrap\Raw Data\LOCAL RELEASE\040318_local_soma_NL1\cell4_AshleyFile.mat';
load(filename);
savename = filename(1:end-4);

%%
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];
if isempty(aa.M)
    M = aa.plotDensityperTime([distances]);
end
figure;
%plot intensity density as a function of time norm to the max somatic intensity
%density
frame_120min = 70;
plot(aa.M.areanormintensity'./aa.M.areanormintensity(1,frame_120min)')
M1 = aa.M.areanormintensity'./aa.M.areanormintensity(1,frame_120min)'

% % plot the intensity density per time norm to the max intensity
% density for each distance
figure;
M3 =  aa.M.areanormintensity';
plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,frame_120min)')
M2 = aa.M.areanormintensity'./aa.M.areanormintensity(:,frame_120min)';
xlswrite(fullfile([savename '_results_areanormintensity']),M3);
xlswrite(fullfile([savename '_results_norm2soma']),M1);
xlswrite(fullfile([savename '_results_norm2each']),M2)
close all



