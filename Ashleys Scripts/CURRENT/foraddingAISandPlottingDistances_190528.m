%% for adding in an AIS region to Ashley files and re-calculating distances
clear all;
close all;
filename = 'C:\Users\bourkea\Dropbox\Data from Cloud\091319_GluA1_activitydependence\TTX\3_AshleyFile.mat';
load(filename);
aa.cellFill.selectAIS(aa.cellFill.image(:,:,end)) %--- select AIS
%--- check if mask looks good AIS
aa.cellFill.AIS_mask.*aa.cleanedcargomask
%--- set and calculate distance
pxsize = 0.114374*2;
d1 = 5/pxsize;
d2 = 40/pxsize;
d3 = 200/pxsize;
distances = [d1 d2 d3];
aa.imagingparams.postrelease.framerate = 2.0;
% possible strings
possiblestrings = {'Total','No AIS','AIS only'};
%%
%M = aa.plotDensityperTime([distances],possiblestrings{1});
M = aa.plotDensityperTime([distances],possiblestrings{1});
save(filename,'aa','-v7.3');
%% to look at heatmap
aa.imagingparams.postrelease.framerate = 2.0;
aa.cargo_heatmap = [];
imgparam.maxtime = 120;
h = aa.plotCargoHeatMap(1,imgparam);
%% if you want to save raw intensity data to normalize at a later time
% aa.M_AIS.rawintensity
% aa.M_noAIS.rawintensity
aa.M.rawintensity
%%
clear M1 M2
figure;
plot(M.areanormintensity'./aa.M.areanormintensity(1,51)')
M1 = M.areanormintensity'./aa.M.areanormintensity(1,51)'

% % plot the intensity density per time norm to the max intensity
% density for each distance
figure;
plot(M.areanormintensity'./M.areanormintensity(:,51)')
M2 = M.areanormintensity'./M.areanormintensity(:,51)'


%figure; hold on;
%plot(M.areanormintensity'./M.areanormintensity(:,30)','r')
%plot(aa.M.areanormintensity'./aa.M.areanormintensity(:,30)','g')
%%
currmask1 = aa.distmask<=d1;
currmask = aa.distmask>d1 & aa.distmask<=d2;
currmask3 = aa.distmask>d2
currmask.*aa.cellFill.AIS_mask
maskd3 = dipshow(currmask.*aa.cellFill.AIS_mask)