% copy ashley files over to computer:
savetopdir = 'C:\Users\schwsama\Documents\Data\zapERtrap';
startdir = 'Z:\Sam\MJK_zapERtrap_for_sam';
[folders2make, pathname] = uigetfile(startdir,'Select filename list file to transfer','Multiselect','on');

for ff = 1:numel(folders2make)
   file2load = fullfile(pathname,folders2make{ff});
   savedir = fullfile(savetopdir,folders2make{ff}(1:end-4));
   if ~exist(savedir)
       mkdir(savedir);
   end
   filenames = load(file2load);
   wb = waitbar(0,'copying files...');
   for currfile = 1:numel(filenames.files)
       out = strsplit(filenames.files{currfile},filesep);
       newfilename = strrep([out{end-1},'_',out{end}],' ','_');
       newfile = fullfile(savedir,newfilename);
       copyfile(filenames.files{currfile},newfile);      
       waitbar(currfile/numel(filenames.files),wb);
   end  
   close(wb)
end