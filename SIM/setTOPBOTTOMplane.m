%% this script is for setting the top and bottom slice
% set the file name
filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\PSD95ib488_gephyrin561_Abeta647_003_Reconstructed_SIM.mat';
load(filename);
% show the file to look through
dipshow(s.ch1.image)
GeneralAnalysis.viewMaskOverlay(s.ch1.image,s.ch1.mask);

% set the planes
obj.planeTOP = []; % -- set the top plane here
obj.planeBOTTOM = []; % -- set the bottom plane here

% save the file
save(filename,'obj');





