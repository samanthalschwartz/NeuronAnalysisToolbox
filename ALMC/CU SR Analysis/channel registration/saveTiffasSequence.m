% load images and save as sequence files
[filename, pathname] = uigetfile('D:\WDKennedyLabHDDBackup\Projects\Project Cry2Olig-Gephyrin\SR\191031 controls and samples\bead images for registration\*.tif','Multiselect','on');
for ff = 1:numel(filename)
   clear sequence;
   sequence = loadtiff(fullfile(pathname,filename{ff}));
%    sequence(:,512:end) = sequence(:,512:end)/4;
   savename = strrep(filename{ff},'.','-');
   save(fullfile(pathname,savename),'sequence');
end