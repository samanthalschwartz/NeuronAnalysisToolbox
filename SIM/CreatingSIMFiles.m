ids = [2 1 3];
channelorderingstr = {'chABeta','PSD95ib','Bassoon'}; % channel abeta, channel 1, channel2
dirname = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\AB_Controls_040819';
savedir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819';
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec',dirname);
wb = waitbar(0);
for ff= 1:numel(filepath)
obj  = SIM();
obj.channelordering = ids;
obj.channelorderingstr = channelorderingstr;
obj.loadNDfile(filepath{ff});
obj.save('C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819');
obj.setimage();
waitbar(ff/numel(filepath),wb);
end
close(wb)
%%
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec',savedir);
wb = waitbar(0);
for ff= 1:numel(filepath)
    load(filepath{ff});
    % -- this is in dipimage values! (one less than image j plane)
    g = dipshow(obj.ch1.rawimage,'log')
    while(ishandle(g))
        try
                w = waitforbuttonpress;
        catch
        end
    end
    opts.WindowStyle='normal';
    prompt = {'Bottom Plane','Top Plane'};
    title = 'select the planes';
    dims = [1 35];
    definput = {'',''};
    answer = inputdlg(prompt,title,dims,definput,opts);
    obj.planeBOTTOM = answer{1};
    obj.planeTOP = answer{2};
    obj.setimage();
    obj.save();
    waitbar(ff/numel(filepath),wb);
end
close(wb)

%% run all the masking, simulation etc
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819');
wb = waitbar(0);
for ff= 1:numel(filepath)
    tic
load(filepath{ff});
obj.make_cellmask;
obj.make_maskchX;
obj.make_maskch1;
obj.make_maskch2;
disp('Finished Making Masks');
obj.make_distancemasks;
obj.measurements = [{'size'},   {'sum'}, {'Gravity'}];
obj.measure_AB;
obj.simulationAbeta(20);
obj.calculateNumberDensityCOM;
obj.calculateNumberDensityCOM(0,1);
obj.save(filepath{ff});
waitbar(ff/numel(filepath),wb);
clear obj;
toc
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


