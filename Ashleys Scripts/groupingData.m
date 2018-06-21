savedir = 'Z:\Sam\MJK_zapERtrap_for_sam';
savename = 'Ashley062018_all_GluA1_Global';

files = uipickfiles('Prompt','select all the AshleyFile(s)','FilterSpec',savedir);

%% 
save(fullfile(savedir,savename),'files');
xlswrite(fullfile(savedir,'tempnames.xls'),files')