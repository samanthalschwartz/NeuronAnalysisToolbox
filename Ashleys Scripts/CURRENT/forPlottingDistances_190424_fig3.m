close all;
clear all;
filename = '\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\GLOBAL RELEASE\NL1\050118\TIFF files\cell1_AshleyFile.mat';
load(filename);
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];

M = aa.plotDensityperTime([distances]);

% show distance 1 area
currmask = aa.distmask<=d1;
maskd1 = dipshow(currmask)

% show distance 2 area
currmask = aa.distmask>d1 & aa.distmask<=d2;
maskd2 = dipshow(currmask)

% show distance 3 area
currmask = aa.distmask>d2 & aa.distmask<=d3;
maskd3 = dipshow(currmask)
save(fullfile([savename '_AshleyFile.mat']), 'aa'); 

save(filename,'aa');
%% can just run this section if you've already calculated everything and just want to plot
% figure;
% %plot density as a function of time
% plot(aa.M.areanormintensity')
figure;
%plot intensity density as a function of time norm to the max somatic intensity
%density
plot(aa.M.areanormintensity'./aa.M.areanormintensity(1,end)')
M1 = aa.M.areanormintensity'./aa.M.areanormintensity(1,end)'

% % plot the intensity density per time norm to the max intensity
% density for each distance
figure;
plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,end)')
M2 = aa.M.areanormintensity'./aa.M.areanormintensity(:,end)';
