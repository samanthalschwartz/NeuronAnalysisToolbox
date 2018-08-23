startdir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\021318_NL1_cells2_3\';
[filename, pathname] = uigetfile(fullfile(startdir,'*.mat'),...
    'Pick the AshleyFiles to Calculate Distances','MultiSelect', 'off');

load(fullfile(pathname,filename));
aa.cleanSurfaceCargoMask_Manual(1);
% aa.cellFill.selectFullSoma();

save(fullfile(pathname,filename),'aa');
dipshow(aa.cleanedcargomask)