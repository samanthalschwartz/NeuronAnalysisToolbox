% this script is for setting the top and bottom slice
filename = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\PSD95ib488_gephyrin561_Abeta647_003_Reconstructed_SIM.mat';
load(filename);
dipshow(s.ch1.image)
GeneralAnalysis.viewMaskOverlay(s.ch1.image,s.ch1.mask);

%%
obj.planeTOP = []; % -- set the top plane here
obj.planeBOTTOM = []; % -- set the bottom plane here
save(filename,'obj');





