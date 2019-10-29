%% for pre/post things only --- don't really need anymore
%% select pre/post ROIs
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819');
for ff= 1:numel(filepath)
load(filepath{ff});
obj.selectPrePostROI;
% obj.selectPrePostROI(1); %for saving new set of ROIs as moreselectiveROIs
obj.save;
end
%% calculate results
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819');
% wb = waitbar(0,'Looping Through Files');
prepostdis_list_all = cell(1,numel(filepath));
numabeta_aroundsynapse_all = cell(1,numel(filepath));
plotvals_all = cell(1,numel(filepath));
for ff= 1:numel(filepath)
    clear obj;
    load(filepath{ff});
    obj.abetaDensityAlongPrePost('moreselectiveROIs');
    obj.save(filepath{ff});

    prepostdis_list_all{ff} =  obj.results.prepostdis_list;
    numabeta_aroundsynapse_all{ff} = obj.results.numabeta_aroundsynapse;
    plotvals_all{ff} = obj.results.plotvals;
end
%% Make some plots
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819');
prepostdis_list_all = cell(1,numel(filepath));
numabeta_aroundsynapse_all = cell(1,numel(filepath));
plotvals_all = cell(1,numel(filepath));
numabeta = [];
all_dists = [];
plotsvals = [];
figure;hold on;


for ff= 1:numel(filepath)
    clear obj;
    load(filepath{ff});
    prepostdis_list_all{ff} =  obj.results_moreselective.prepostdis_list;
    numabeta_aroundsynapse_all{ff} = obj.results_moreselective.numabeta_aroundsynapse;
    plotvals_all{ff} = obj.results_moreselective.plotvals;
%     prepostdis_list_all{ff} =  obj.results.prepostdis_list;
%     numabeta_aroundsynapse_all{ff} = obj.results.numabeta_aroundsynapse;
%     plotvals_all{ff} = obj.results.plotvals;
end

cols = lines(numel(plotvals_all));
for ii = 1:numel(plotvals_all)
   all_dists = [all_dists; prepostdis_list_all{ii}'];
   numabeta = [numabeta; numabeta_aroundsynapse_all{ii}'];
   plotsvals = cat(1,plotsvals,plotvals_all{ii}');
   scatter(plotvals_all{ii}(1,:),plotvals_all{ii}(2,:),'MarkerFaceColor',cols(ii,:),'MarkerEdgeColor',cols(ii,:))
end
figure; scatter(plotsvals(:,1).*obj.XYpxsize,plotsvals(:,2).*obj.XYpxsize,'*'); pbaspect([1 1 1])
ylabel('Distance from Center of POST-PRE')
xlabel('Distance Perpendicular to Synapse')
xlim([0 1]); ylim([-1 1]); pbaspect([1 1 1])

hold on; plot([0:40],repmat(0,1,41),'w')

plotsvals(plotsvals(:,1)==0,:) = [];
[gca,N] = scatter2heatmap(plotsvals(:,1).*obj.XYpxsize,plotsvals(:,2).*obj.XYpxsize,45);
% colormap('hot')

xlim([0 1]); ylim([-1 1]); pbaspect([1 1 1])
hold on; plot([0:40],repmat(0,1,41),'w','LineWidth',2)
ylabel('Distance from Center of POST-PRE')
xlabel('Distance Perpendicular to Synapse')
colorbar

figure; histogram(all_dists.*obj.XYpxsize,20); title('Distance from Pre to Post');
xlabel('Distance in microns')
figure; histogram(numabeta,20); title('Number of Abeta COM within region Analyzed');