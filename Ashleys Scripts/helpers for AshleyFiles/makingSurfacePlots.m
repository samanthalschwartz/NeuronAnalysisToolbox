pathname = 'C:\Users\schwsama\Documents\Data\zapERtrap\all_NL1_global';
filename = '021318_NL1_cells2_3_cell2_global_AshleyFile';
% startdir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\021318_NL1_cells2_3\';
% [filename, pathname] = uigetfile(fullfile(startdir,'*.mat'),...
%     'Pick the AshleyFiles to Calculate Distances','MultiSelect', 'off');

load(fullfile(pathname,filename));
if ~isempty(aa.cleanedcargomask)
dipisosurface_samchange(aa.cleanedcargomask)
else
    dipisosurface_samchange(aa.cellFill.mask*aa.surfaceCargo.mask)
end

set(gca,'FontSize',16)
%%--- to change colormaps
colormap(hsv);
colormap(flip(jet(60)));

%%--- to change aspect ratio of plot
xsize =1;
ysize = 1;
zsize = .2;
pbaspect([xsize ysize zsize])
set(gca,'clim',[ 7 33])

% %%-- to add labels
% xlabel('x-axis')
% ylabel('y-axis')
zlabel('Time (frames')
