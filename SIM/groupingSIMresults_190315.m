% Grouping Binned Information
% Make Page for each protein
% - FileName
% - Bins
% - NumAbeta Raw
% - NumObj
% - NormNumAbeta

close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','F:\Hannah Dropbox SIM data\SIM_Files');
ch1vals = cell(2,numel(filepath));
ch2vals = cell(2,numel(filepath));

for ff = 1:numel(filepath)
load(filepath{ff});
ch1vals{1,ff} = obj.ch1.results.bins;
ch1vals{2,ff} = obj.ch1.results.numabeta/obj.ch1.results.numobj;
ch2vals{1,ff} = obj.ch2.results.bins;
ch2vals{2,ff} = obj.ch2.results.numabeta/obj.ch2.results.numobj;
end

figure; hold on;
for ff = 1:numel(filepath)
plot(ch1vals{1,ff},ch1vals{2,ff});
end

figure; hold on;
for ff = 1:numel(filepath)
plot(ch2vals{1,ff},ch2vals{2,ff});
end