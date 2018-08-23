%% call this to calculate and plot distances of an Ashley File
% this will prompt to select multiple files to run. Currently need to only
% select files from the same directory 

startdir = 'C:\Users\schwsama\Documents\Data\zapERtrap';
% [filename, pathname] = uigetfile(fullfile(startdir,'*.mat'),...
%     'Pick the AshleyFiles to Calculate Distances (can pick more than 1)','MultiSelect', 'on');
filename = uipickfiles('Prompt','Pick all the AshleyFiles to Calculate Distances (can pick more than 1)');
if ~iscell(filename)
    filesize = 1;
    filename = {filename};
else
    filesize = numel(filename);
end

plotflag = 1; % use 1 to keep .tif file of plots each frame, 0 to not
for ii = 1:filesize
    [pathname,currfilename,EXT] = fileparts(filename{ii});
    plotsavename = currfilename;
    plotsavedir = fullfile(pathname,plotsavename);
    if ~exist(plotsavedir,'dir')
        mkdir(fullfile(plotsavedir));
    end
    clear aa;
    disp('Loading File.....');
    load(fullfile(filename{ii}));
    disp(['File Loaded: ' fullfile(filename{ii})]);
    aa.calculateSurfaceCargoDistances(plotflag,plotsavedir);
    save(fullfile(filename{ii}),'aa');
    disp('File Saved');
end
