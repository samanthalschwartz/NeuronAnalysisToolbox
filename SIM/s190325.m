close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
for ff= 1:numel(filepath)
load(filepath{ff});
obj.simulationAbeta(20);
obj.calculateNumberDensityCOM;
obj.calculateNumberDensityCOM(0,1);
obj.save;
end
% 
figure; hold on;
plot(obj.ch1.results.bins,obj.ch1.results.numabeta./obj.ch1.results.volume,'Color',[0 0 1])
plot(obj.ch1.abetaSIM.results.bins,obj.ch1.abetaSIM.results.numabeta./obj.ch1.abetaSIM.results.volume,'Color',[1 0 0])

figure; hold on;
plot(obj.ch2.results.bins,obj.ch2.results.numabeta./obj.ch2.results.volume)
plot(obj.ch2.abetaSIM.results.bins,obj.ch2.abetaSIM.results.numabeta./obj.ch2.abetaSIM.results.volume)