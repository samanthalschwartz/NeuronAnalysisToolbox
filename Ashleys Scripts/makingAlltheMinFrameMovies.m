%% this is set up to look inside all folders within toptopdir for AshleyFiles
toptopdir = 'Z:\Sam\MJK_zapERtrap_for_sam\AMB_local\101117_TfR-smFP_HA-GFP-DHFR_local';
% first go into all subdirectories in toptopdir folder
temp1 = dir(toptopdir);
temp2 = temp1(arrayfun(@(x) x.isdir,temp1));
folds = temp2(3:end); % get rid of . and ..

% now go through each folder in folds and 
% 1) find all the Ashley Files
% 2) load them in and create movie files
% 3) save the movie with the correct name

framelag = 4;
for dd = 1:numel(folds)
    currdir = fullfile(toptopdir,folds(dd).name);
   files = dir(fullfile(currdir, '*_AshleyFile.mat'));
   if ~isempty(files)
   for ff = 1:numel(files)
       filename = files(ff).name;
      load(fullfile(currdir,filename));
      nameinfo = strsplit(filename,'_');
      savename = fullfile(currdir,[nameinfo{1} '_minframemovie']);
      [h,~] = aa.plot_cargo_minFrameMovie(framelag,savename);
      close(h);
   end
   end
end
%% this is set up to just look inside topdir for AshleyFiles
topdir = {'Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\051318 GluA1 insertion',...
    '',...
    '',...
    '',...
    };
% now go through each folder in folds and 
% 1) find all the Ashley Files
% 2) load them in and create movie files
% 3) save the movie with the correct name
for tt = 1:numel(topdir)
framelag = 4;
files = dir(fullfile(topdir{tt}, '*_AshleyFile.mat'));
if ~isempty(files)
for ff = 1:numel(files)
    filename = files(ff).name;
    load(fullfile(topdir{tt},filename));
    posinfo = strfind(filename,'_');
    name = filename(1:(posinfo(2)-1));
    savename = fullfile(topdir{tt},[name '_minframemovie']);
    [h,~] = aa.plot_cargo_minFrameMovie(framelag,savename);
    close(h);
end
end
end
