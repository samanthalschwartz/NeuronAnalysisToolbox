% calling SEP analysis from SEP class
% filepath = 'E:\Matt Becker Data (For Review)\SEPGlua1_mch';
filepath = 'E:\Matt Becker Data (For Review)\SEPGlua1_mch2aTeNT';
savedir = fullfile(filepath,'SEPfiles');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
files = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);


for ff= 1:numel(files)
filename = files{ff};
display(['Analyzing File: ' filename]);
[path,name,ext] = fileparts(filename);
% load and calculate things
sp = SEP(); sp.filepath = path; sp.filename = [name, ext];
disp('Loading Images...');
sp.loadimages;
disp('Masking Images...');
sp.make_mask_cellfill;
sp.make_mask_sep;
sp.make_mask_sep_fixed;
disp('Calculating individual mask intensities...');
sp.calculate_sepintensities;
disp('Calculating overall mask intensities...');
sp.calculateSpineShaftIntensities;
sp.loadIJROIs1;
sp.loadIJROIs2;
% save and plot things
disp('Saving SEP File...');
sp.saveSEP(savepath);
h = sp.viewSEPMaskFixedwOldRois();
% saveas(h,fullfile(savepath,[sp.filename(1:end-4) '_allROIs']),'fig');
saveas(h,fullfile(savedir,[sp.filename(1:end-4) '_allROIs']),'png');
close(h);
f = sp.calculate_sepIHeatMap;
set(f,'Position', [2146,49,1193,935]);
% saveas(f,fullfile(savepath,[sp.filename(1:end-4) '_heatmap']),'fig');
saveas(f,fullfile(savedir,[sp.filename(1:end-4) '_heatmap']),'png');
close(f);
end

%%
% this section is to loop through all traces of a file and extract only
% those that respond above some threshold (1.1 fold). Count the fraction of
% these and then make a heatmap for it. 

% first load in all the data (this is slowish)
datadir = 'E:\Matt Becker Data (For Review)\SEPfiles\-TeNT';
% files = uipickfiles('Prompt','Pick Files','FilterSpec',datadir);
files_noT = dir(fullfile(datadir,'*.mat'));
alltraces_noT = cell(numel(files_noT),1);
for ff = 1:numel(files_noT)
sep = load(fullfile(datadir,files_noT(ff).name));
[f] = sep.obj.calculate_sepIHeatMap;
close(f)
alltraces_noT{ff} = sep.obj.ord_trace;
end

datadir = 'E:\Matt Becker Data (For Review)\SEPfiles\+TeNT';
% files = uipickfiles('Prompt','Pick Files','FilterSpec',datadir);
files_T = dir(fullfile(datadir,'*.mat'));
alltraces_T = cell(numel(files_T),1);
for ff = 1:numel(files_T)
sep = load(fullfile(datadir,files_T(ff).name));
[f] = sep.obj.calculate_sepIHeatMap;
close(f)
alltraces_T{ff} = sep.obj.ord_trace;
end
%%
% now do some things
minval = 1.1; timewindow = 1:25;

% + TeNT 
allvalsT = cell2mat(alltraces_T);
testT = sum(allvalsT(:,timewindow)>minval,2);
goodvalsT = allvalsT(testT>1,:);
goodvalsTm = movmean(goodvalsT,4,2);

allsum = sum(goodvalsTm,2);
[~, ordx] = sort(allsum, 'descend');
ord_traceT = goodvalsTm(ordx,:);
% % ord_traceT = movmean(goodvalsTm(ordx,:),3,2);


figure;
hm1 = heatmap(ord_traceT(500:end,:))
hm1.GridVisible = 'off';
hm1.Colormap = jet(50);
hm1.ColorLimits = [0 2.3];
hm1.Title = '+TeNT'

% - TeNT
allvals_noT = cell2mat(alltraces_noT);
testnoT = sum(allvals_noT(:,timewindow)>minval,2);
goodvals_noT = allvals_noT(testnoT>1,:);
goodvals_noTm = movmean(goodvals_noT,4,2);
goodvals_noTm(312,:) = [];

allsum = sum(goodvals_noTm,2);
[~, ordx] = sort(allsum, 'descend');
ord_trace_noT = goodvals_noTm(ordx,:);
% fix dropout frame: 13
% ord_trace_noT(:,13) = ord_trace_noT(:,12);
% ord_trace_noT(:,24) = ord_trace_noT(:,23);
% ord_trace_noT(:,34) = ord_trace_noT(:,33);
figure;
hm2 = heatmap(ord_trace_noT(500:end,:))
hm2.GridVisible = 'off';
hm2.Colormap = jet(50);
hm2.ColorLimits = [0 2.3];
hm2.Title = '-TeNT';
%%
hm1.Colormap = hot(50);
hm2.Colormap = hot(50);
hm1.ColorLimits = [0 2.2];
hm2.ColorLimits = [0 2.2];

saveas(1,'E:\Matt Becker Data (For Review)\Figures\-TenT_heatmap','tif');
saveas(2,'E:\Matt Becker Data (For Review)\Figures\+TenT_heatmap','tif');
%%
% calculate fraction of ROIs in different grouping
minval = 1.2;
allvalspreT = cell2mat(alltraces_T);
testT = sum(allvalspreT>minval,2);
allvalsT = allvalspreT(testT>1,:);

allvalsprenoT = cell2mat(alltraces_noT);
testT = sum(allvalsprenoT>minval,2);
allvalsnoT = allvalsprenoT(testT>1,:);




timewindow = 1:25;

newT = sort(allvalsT(:,timewindow),2,'descend');
newT(newT(:,1)>5,:) = [];
newnoT = sort(allvalsnoT(:,timewindow),2,'descend');
newnoT(newnoT(:,1)>5,:) = [];
figure; hold on;
ecdf(newT(:,3))
ecdf(newnoT(:,3))
xlabel('Max Fold Increase (within first 25 frames)','FontSize',16)
ylabel('Cumulative Probability','FontSize',16);

lastframewindow = 3;
newTmax = max(allvalsT,[],2);
newT = mean(allvalsT(:,end-lastframewindow:end),2);
newnoTmax = max(allvalsnoT,[],2);
newnoT = mean(allvalsnoT(:,end-lastframewindow:end),2);
badids=  newnoT(:,1)>5;
newnoT(badids,:) = [];
newnoTmax(badids,:) = [];

figure; hold on;
ecdf(newT-newTmax)
ecdf(newnoT-newnoTmax)
xlim([-0.5 0])
xlabel('Max Fold Increase (within first 25 frames)','FontSize',16)
ylabel('Cumulative Probability','FontSize',16);






%%
% plot curves;

tmean = mean(goodvalsTm,1);
tstd = std(goodvalsTm,1)/sqrt(size(goodvalsTm,1));
goodvals_noTm(isnan(goodvals_noTm(:,1)),:) = []


ntmean = mean(goodvals_noTm,1);
ntstd = std(goodvals_noTm,1)/sqrt(size(goodvals_noTm,1));


figure; hold on;
errorbar(tmean,tstd)
plot(tmean);hold on;
errorbar(ntmean,ntstd)
plot(ntmean);

