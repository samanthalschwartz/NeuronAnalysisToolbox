% re-run all the SIM file masking in these folders 
folder1 = 'G:\Hannah Dropbox SIM data\SIM_Files\071917';
folder2 = 'G:\Hannah Dropbox SIM data\SIM_Files\071117';

folders = {folder1,folder2};

for ff = 1:numel(folders)
   inputsavedir = folders{ff};
   files = dir(fullfile(folders{ff},'*.mat'));
   wb = waitbar(0,folders{ff});
   for ii = 1:numel(files)
    clear obj;
    load(fullfile(folders{ff},files(ii).name));
    obj.make_maskchAB;
    obj.make_maskch1;
    obj.make_maskch2;
    obj.save(inputsavedir);
    waitbar(ii/numel(files),wb);
   end
end

