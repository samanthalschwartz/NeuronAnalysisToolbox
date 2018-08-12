startdir = 'Z:\Sam\MJK_zapERtrap_for_sam\';
[filename, pathname] = uigetfile(fullfile(startdir,'*.mat'),...
    'Pick the AshleyFiles to Calculate Distances','MultiSelect', 'off');

load(fullfile(pathname,filename));
aa.cellFill.selectFullSoma();
save(fullfile(pathname,filename),'aa');
