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