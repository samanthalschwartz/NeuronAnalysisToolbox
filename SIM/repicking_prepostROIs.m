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