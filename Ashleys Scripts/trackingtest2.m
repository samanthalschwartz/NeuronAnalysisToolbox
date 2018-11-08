savedir = 'C:\Users\schwsama\Documents\Data\zapERtrap\for SPT';
load('C:\Users\schwsama\Documents\Data\zapERtrap\all_NL1_global\040318_globalrelease_NL1_cell1_AshleyFile.mat');
sequence = aa.surfaceCargo.image;
save(fullfile(savedir,'040318_globalrelease_NL1_cell1_surfaceCargo'),'sequence');
sequence = aa.TfR.image;
save(fullfile(savedir,'040318_globalrelease_NL1_cell1_NL1channel'),'sequence');
