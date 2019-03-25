
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\sammy\Desktop\Brooke SIM\SIM_Files\071917');

for ff= 1:numel(filepath)
    clear obj ch1_lab abeta_dist ch1_N ch2_lab ch2_N ch1_msr ch2_msr
load(filepath{ff});
ch1_lab = label(obj.ch1.mask,1);
ch2_lab = label(obj.ch2.mask,1);
abeta_dist = obj.abeta.distance_mask;

ch1_msr = measure(ch1_lab,abeta_dist,{'MinVal'});
ch2_msr = measure(ch2_lab,abeta_dist,{'MinVal'});


[ch1_N,ch1_edges] = histcounts(ch1_msr.MinVal,'BinWidth',1,'Normalization','cdf');
[ch2_N,ch2_edges] = histcounts(ch2_msr.MinVal,'BinWidth',1,'Normalization','cdf');

obj.ch1.results.mindist2abeta.N = ch1_N;
obj.ch1.results.mindist2abeta.edges = ch1_edges;
obj.ch2.results.mindist2abeta.N = ch2_N;
obj.ch2.results.mindist2abeta.edges = ch2_edges;
obj.save(filepath{ff});
end

%%


% close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\sammy\Desktop\Brooke SIM\SIM_Files\071917');
% allfilepath = filepath;
% filepath = allfilepath(4:6);
ch1vals = cell(2,numel(filepath));
ch2vals = cell(2,numel(filepath));

for ff = 1:numel(filepath)
load(filepath{ff});
ch1vals{1,ff} = obj.ch1.results.mindist2abeta.N;
ch1vals{2,ff} = obj.ch1.results.mindist2abeta.edges(1:end-1);
ch2vals{1,ff} = obj.ch2.results.mindist2abeta.N;
ch2vals{2,ff} = obj.ch2.results.mindist2abeta.edges(1:end-1);
end

figure; hold on;
for ff = 1:numel(filepath)
plot(ch1vals{2,ff},ch1vals{1,ff},'Color',[0 0 1]);
plot(ch2vals{2,ff},ch2vals{1,ff},'Color',[0 1 0]);
end
%%
figure; plot(obj.ch1.results.mindist2abeta.edges(1:end-1),obj.ch1.results.mindist2abeta.N,'Color',[0 1 0]);
hold on; plot(obj.ch2.results.mindist2abeta.edges(1:end-1),obj.ch2.results.mindist2abeta.N,'Color',[0 0 1]);
xlim([0 5]);ylim([0 1]);