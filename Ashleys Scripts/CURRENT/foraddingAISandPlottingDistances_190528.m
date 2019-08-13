%% for adding in an AIS region to Ashley files and re-calculating distances
clear all;
filename = '\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\LOCAL RELEASE\030618_localexpt_NL1\cell1_AshleyFile.mat'
load(filename);
aa.cellFill.selectAIS(aa.surfaceCargo.image) %--- select AIS
%--- check if mask looks good AIS
aa.cellFill.AIS_mask.*aa.cleanedcargomask
%--- set and calculate distance
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];
% possible strings
possiblestrings = {'Total','No AIS','AIS only'};
%%
M = aa.plotDensityperTime([distances],possiblestrings{1});
save(filename,'aa');
%% to look at heatmap
aa.plotCargoHeatMap()
%%
clear M1 M2
figure;
plot(M.areanormintensity'./aa.M.areanormintensity(1,52)')
M1 = M.areanormintensity'./aa.M.areanormintensity(1,52)'

% % plot the intensity density per time norm to the max intensity
% density for each distance
figure;
plot(M.areanormintensity'./M.areanormintensity(:,52)')
M2 = M.areanormintensity'./M.areanormintensity(:,52)'

figure; hold on;
plot(M.areanormintensity'./M.areanormintensity(:,70)','r')
plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,70)','g')
%%
currmask = aa.distmask>d2 & aa.distmask<=d3;
maskd3 = dipshow(currmask)