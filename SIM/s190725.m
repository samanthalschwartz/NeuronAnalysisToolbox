close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\061019');
wb = waitbar(0);
for ff= 1:numel(filepath)
    tic
load(filepath{ff});
disp('Making calculations...');
obj.simulationAbeta(20);
obj.calculateNumberDensityCOM;
obj.calculateNumberDensityCOM(0,1);
obj.save(filepath{ff});
waitbar(ff/numel(filepath),wb);
clear obj;
toc
end
close(wb);