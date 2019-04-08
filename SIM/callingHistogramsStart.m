clear all
datadir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\072617';
files = dir(fullfile(datadir,'PSD95ib488*'));

allvals = [];
for ff = 1:numel(files)
    clear obj;
    savepath = fullfile(datadir,files(ff).name);
    load(savepath);
    obj.savepath = savepath;
    currvals= obj.mindistCh1tofirstAbeta(1);
    obj.save;
    allvals = cat(2,allvals,currvals);
end
allvals = allvals(allvals<10);
figure; histogram(allvals.*obj.XYpxsize,50,'Normalization','Probability');
xlabel('Distance in microns');
title(['Histogram of Min distance from ABeta to ' obj.channelorderingstr{2}]);

allvals = [];
for ff = 1:numel(files)
    clear obj;
    savepath = fullfile(datadir,files(ff).name);
    load(savepath);
    obj.savepath = savepath;
    currvals= obj.mindistCh2tofirstAbeta(1);
    obj.save;
    allvals = cat(2,allvals,currvals);
end
allvals = allvals(allvals<10);
figure; histogram(allvals.*obj.XYpxsize,50,'Normalization','Probability');
xlabel('Distance in microns');
title(['Histogram of Min distance from ABeta to ' obj.channelorderingstr{3}]);

%%
clear all
datadir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071917';
files = dir(fullfile(datadir,'PSD95488*'));

allvals = [];
for ff = 1:numel(files)
    clear obj;
    savepath = fullfile(datadir,files(ff).name);
    load(savepath);
    obj.savepath = savepath;
    currvals= obj.mindistCh1tofirstAbeta;
    obj.save;
    allvals = cat(2,allvals,currvals);
end
allvals = allvals(allvals<10);
figure; histogram(allvals.*obj.XYpxsize,30,'Normalization','Probability');
xlabel('Distance in microns');
title(['Histogram of Min distance from ABeta to ' obj.channelorderingstr{2}]);

allvals = [];
for ff = 1:numel(files)
    clear obj;
    savepath = fullfile(datadir,files(ff).name);
    load(savepath);
    obj.savepath = savepath;
    currvals= obj.mindistCh2tofirstAbeta;
    obj.save;
    allvals = cat(2,allvals,currvals);
end
allvals = allvals(allvals<10);
figure; histogram(allvals.*obj.XYpxsize,30,'Normalization','Probability');
xlabel('Distance in microns');
title(['Histogram of Min distance from ABeta to ' obj.channelorderingstr{3}]);