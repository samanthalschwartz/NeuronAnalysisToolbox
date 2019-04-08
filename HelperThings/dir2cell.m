function files_out = dir2cell(topdir,filestr)
% function that mirrors dir call, but returns files as cell array
files = dir(fullfile(topdir,filestr));
files_unsorted = cell(1,numel(files));
for ff = 1:numel(files)
    files_unsorted{ff} = fullfile(topdir,files(ff).name);
end
files_out = natsortfiles(files_unsorted);
end
