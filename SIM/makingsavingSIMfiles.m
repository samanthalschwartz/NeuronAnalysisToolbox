% select a save directory
savedir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files';
if ~exist(savedir)
    mkdir(savedir)
end
% select the directory to choose files from and the corresponding
% parameters for it
topdir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\080317\500 nM Abeta_PSD95_gephyrin';
ids = [3 1 2];
channelorderingstr = {'chABeta','PSD95i','Gephyrin'}; % channel abeta, channel 1, channel2

% pick the files and then make and save SIM files
files = uipickfiles('FilterSpec',topdir,'Prompt',['Pick Files for ' channelorderingstr{1} ,...
                                            '-' channelorderingstr{2} '-' channelorderingstr{3}]);
for ff = 1:numel(files)
filepath = files{ff};
s  = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
s.dothething();
s.save(savedir);
end

%%
datafolder = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\071117\Gephyrin_Abeta_Bassoon';
savedir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071117';
if ~exist(savedir)
    mkdir(savedir)
end
files1 = dir(fullfile(datafolder,'Geph488_Abeta561_Bassoon647*.nd2'));
for ff= 1:numel(files1)
ids = [2 1 3];
channelorderingstr = {'chABeta','Gephyrin','Bassoon'}; % channel abeta, channel 1, channel2
filepath = fullfile(datafolder,files1(ff).name);
s = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
s.make_masks();
s.save(savedir);
disp(num2str(ff));
end
clear all
datafolder = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\071117\GluA1_Abeta_bassoon';
savedir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071117';
if ~exist(savedir)
    mkdir(savedir)
end
files2 = dir(fullfile(datafolder,'GluA1488_Abeta561_Bassoon647*Reconstructed.nd2'));
for ff= 1:numel(files2)
ids = [2 1 3];
channelorderingstr = {'chABeta','GluA1','Bassoon'}; % channel abeta, channel 1, channel2
filepath = fullfile(datafolder,files2(ff).name);
s = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
s.make_masks();
s.save(savedir);
disp(num2str(ff));
end
clear all
datafolder = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\071117\PSD95_Abeta_Bassoon';
savedir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071117';
if ~exist(savedir)
    mkdir(savedir)
end
files3 = dir(fullfile(datafolder,'PSD95488_Abeta561_Bassoon647*'));
for ff= 1:numel(files3)
ids = [2 1 3];
channelorderingstr = {'chABeta','PSD95','Bassoon'}; % channel abeta, channel 1, channel2
filepath = fullfile(datafolder,files3(ff).name);
s = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
obj.make_masks();
s.save(savedir);
end

%%
% datafolder = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\to fix from ndfiles\071917';
savedir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071917';
% files3 = dir(fullfile(datafolder,'Gephyrin488_Abeta561_Bassoon647*'));
orginaldir = 'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\071917\Gephyrin_Bassoon';
files3 = uigetfile('C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\071917\Gephyrin_Bassoon','MultiSelect', 'on')
for ff= 1:numel(files3)
ids = [2 1 3];
channelorderingstr = {'chABeta','Gephyrin','Bassoon'}; % channel abeta, channel 1, channel2
filepath = fullfile(orginaldir,files3{ff});
s = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
s.make_masks();
s.save(savedir);
end






