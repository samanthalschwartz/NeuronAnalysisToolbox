close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
wb = waitbar(0);
for ff= 1:numel(filepath)
load(filepath{ff});
obj.make_cellmask;
obj.make_maskchAB;
obj.make_maskch1;
obj.make_maskch2;
obj.make_distancemasks;
obj.measurements = [{'size'},   {'sum'}, {'Gravity'}];
obj.measure_AB;
obj.simulationAbeta(20);
obj.calculateNumberDensityCOM;
obj.calculateNumberDensityCOM(0,1);
obj.save(filepath{ff});
waitbar(ff/numel(filepath),wb);
end
close(wb);
% 
figure; hold on;
plot(obj.ch1.results.bins,obj.ch1.results.numabeta./obj.ch1.results.volume,'Color',[0 0 1])
plot(obj.ch1.abetaSIM.results.bins,obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume,'Color',[1 0 0])

figure; hold on;
plot(obj.ch2.results.bins,obj.ch2.results.numabeta./obj.ch2.results.volume)
plot(obj.ch2.abetaSIM.results.bins,obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume)

%%
figure; hold on;
plot(obj.ch2.results.bins,movmean(obj.ch2.results.numabeta./obj.ch2.results.volume,3))
plot(obj.ch2.abetaSIM.results.bins,movmean(obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume,10));
title(obj.channelorderingstr{3});

figure; hold on;
plot(obj.ch1.results.bins,movmean(obj.ch1.results.numabeta./obj.ch1.results.volume,2))
plot(obj.ch1.abetaSIM.results.bins,movmean(obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume,10));
title(obj.channelorderingstr{2});
%%
%%
figure; hold on;
ch1dat = movmean(obj.ch1.results.numabeta./obj.ch1.results.volume,3);
ch1sim =  movmean(obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume,10);
ch1norm = ch1dat./ch1sim;
ch2dat = movmean(obj.ch2.results.numabeta./obj.ch2.results.volume,3);
ch2sim =  movmean(obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume,10);
ch2norm = ch2dat./ch2sim;

plot(obj.ch1.abetaSIM.results.bins,ch1dat./ch1sim); title(obj.channelorderingstr{2});
figure; hold on;
plot(obj.ch2.abetaSIM.results.bins,ch2dat./ch2sim); title(obj.channelorderingstr{3});

assert(isequal(obj.ch1.abetaSIM.results.bins,obj.ch2.abetaSIM.results.bins))
bins = obj.ch1.abetaSIM.results.bins;
A = [[obj.filepath ' Bins']; num2cell(bins)'];
B = [[obj.filepath ' Ch1']; num2cell(ch1norm)];
C = [[obj.filepath ' Ch2']; num2cell(ch2norm)];
matdata = [A,B,C];
xlswrite('G:\Hannah Dropbox SIM data\SIM_Files\071117\data.xlsx',matdata);
%%



%% ch2
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\061019');
wb = waitbar(0)
savename = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071_062819\Bungaro488_ch2';
for ff= 1:numel(filepath)
load(filepath{ff});
ch2dat = movmean(obj.ch2.results.numabeta./obj.ch2.results.volume,3);
ch2sim =  movmean(obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume,10);
ch2norm = ch2dat./ch2sim;

bins = obj.ch2.abetaSIM.results.bins;
out_bins = [[obj.filepath ' Bins']; num2cell(bins)'];
out_data = [[obj.filepath ' ch2']; num2cell(ch2norm)];
if ff == 1
    matdata = [out_bins,out_data];
else
    matdata = [matdata,out_data];
end
waitbar(ff/numel(filepath),wb);
end
close(wb)
xlswrite(savename,matdata);

%%
%% ch1
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071_062819');
savename = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071_062819\Bungaro488_ch1';
wb = waitbar(0)
for ff= 1:numel(filepath)
load(filepath{ff});
ch1dat = movmean(obj.ch1.results.numabeta./obj.ch1.results.volume,3);
ch1sim =  movmean(obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume,10);
ch1norm = ch1dat./ch1sim;

bins = obj.ch1.abetaSIM.results.bins;
out_bins = [[obj.filepath ' Bins']; num2cell(bins)'];
out_data = [[obj.filepath ' ch1']; num2cell(ch1norm)];
if ff == 1
    matdata = [out_bins,out_data];
else
    matdata = [matdata,out_data];
end
waitbar(ff/numel(filepath),wb);
end
close(wb)
xlswrite(savename,matdata);


%%
ch1data = obj.ch1.results.numabeta./obj.ch1.results.volume;
ch1sim = obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume;
ch1bins = obj.ch1.results.bins;
ch2  data = obj.ch2.results.numabeta./obj.ch2.results.volume;
ch2sim = obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume;
ch2bins = obj.ch2.results.bins;

figure; hold on;
plot(ch1bins,movmean(ch1data,3)./movmean(ch1sim,10),'Color',[0 0 1]);
% plot(ch1bins,movmean(ch1sim,10),'Color',[1 0 0]);
title('Ch1');

figure; hold on;
plot(ch2bins,movmean(ch2data,3)./movmean(ch2sim,10),'Color',[0 0 1]);
% plot(ch2bins,movmean(ch2sim,10),'Color',[1 0 0]);
title('Ch2');
