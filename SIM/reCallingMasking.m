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



%%
close all; clear all;
folder1 = 'G:\Hannah Dropbox SIM data\SIM_Files\112117';
folder2 = 'G:\Hannah Dropbox SIM data\SIM_Files\111617';
folder3 = 'G:\Hannah Dropbox SIM data\SIM_Files\080317';
folder4 = 'G:\Hannah Dropbox SIM data\SIM_Files\072617';
folders = {folder1,folder2,folder3,folder4};

for ff = 1:numel(folders)
    inputsavedir = folders{ff};
   files = dir(fullfile(folders{ff},'*.mat'));
   w1b = waitbar(0,folders{ff});
   for ii = 1:numel(files)
    clear obj;
    load(fullfile(folders{ff},files(ii).name));
    obj.make_maskchAB;
    obj.make_maskch1;
    obj.make_maskch2;
    obj.make_distancemasks;
    obj.measurements = [{'size'},   {'sum'}, {'Gravity'}];
    obj.measure_AB;
    obj.calculateNumberDensityCOM;
    obj.save(inputsavedir);
    waitbar(ii/numel(files),w1b);
   end
   close(w1b);
end