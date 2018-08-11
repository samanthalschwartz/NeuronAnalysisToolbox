startdir = 'Z:\Sam\Data\MJK_zapERtrap_for_sam\';
[filename, pathname] = uigetfile(fullfile(startdir,'*.mat'),...
    'Pick the AshleyFiles to Calculate Distances','MultiSelect', 'off');

load(fullfile(pathname,filename));
aa.cleanSurfaceCargoMask_Manual();
save(fullfile(pathname,filename),'aa');