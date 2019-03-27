ids = [3 2 1];
channelorderingstr = {'chABeta','NR1','PSD95'}; % channel abeta, channel 1, channel2
filepath = 'G:\Hannah Dropbox SIM data\112117\NR1_PSD95_500 nM Abeta\NR1488_PSD95561_Abeta647_006_Reconstructed.nd2';

obj  = SIM();
obj.channelordering = ids;
obj.channelorderingstr = channelorderingstr;
obj.loadNDfile(filepath);

% -- this is in dipimage values! (one less than image j plane)
opts.WindowStyle='normal';
prompt = {'Bottom Plane','Top Plane'};
title = 'select the planes';
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,title,dims,definput,opts);
obj.planeBOTTOM = answer{1};
obj.planeTOP = answer{2};
obj.setimage();
obj.save('G:\Hannah Dropbox SIM data\SIM_Files\112117');

%% run all the masking, simulation etc
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
%% stupid way of saving
%% ch2
close all; clear all;
saveExcelfilename = 'G:\Hannah Dropbox SIM data\SIM_Files\071117\PSD95_112117_withNR1.xlsx';
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
wb = waitbar(0)
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
xlswrite(saveExcelfilename,matdata);

%%
%% ch1
close all; clear all;
saveExcelfilename = 'G:\Hannah Dropbox SIM data\SIM_Files\071117\NR1_112117.xlsx';
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','G:\Hannah Dropbox SIM data\SIM_Files');
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
xlswrite(saveExcelfilename,matdata);


