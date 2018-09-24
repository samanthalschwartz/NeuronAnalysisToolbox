
% this script is for setting the top and bottom slice
close all; clear all;
%% this script is for setting the top and bottom slice
% set the file name
>>>>>>> c223574a64519aa83c2ba0d4c537d33d73918a53
filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\PSD95ib488_gephyrin561_Abeta647_003_Reconstructed_SIM.mat';
load(filename);
% show the file to look through
dipshow(s.ch1.image)
GeneralAnalysis.viewMaskOverlay(s.ch1.image,s.ch1.mask);
s.ch1.mask
%
obj.planeTOP = 3; % -- set the top plane here
obj.planeBOTTOM = 11; % -- set the bottom plane here

% set the planes
obj.planeTOP = []; % -- set the top plane here
obj.planeBOTTOM = []; % -- set the bottom plane here

% save the file
save(filename,'obj');


%% this is if you want to set the same top/bottom planes for >1 file at a time
close all; clear all;
topdir = pwd;
planeTOP = 3; % -- set the top plane here
planeBOTTOM = 11; % -- set the bottom plane here

files = uipickfiles('FilterSpec',topdir,'Prompt',['Pick Files for ' channelorderingstr{1} ,...
    '-' channelorderingstr{2} '-' channelorderingstr{3}]);

for ff = 1:numel(files)
    load(files{ff});
    obj.planeTOP = 3; % -- set the top plane here
    obj.planeBOTTOM = 11; % -- set the bottom plane here
    save(filename,'obj');
end

