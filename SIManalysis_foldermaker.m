%% you can input the files/folders manually, otherwise if they are empty you will be prompted:
synfolder = [];
squashpath = [];
csvpath = [];
% --- examples for how to set things manually ---
% synfolder = {'G:\Sam\Data\SIM\180524 testing\gephryn488_cellfillmch_basoon647_002_synROIs',...
%     'G:\Sam\Data\SIM\180524 testing\gephryn488_cellfillmch_basoon647_002_synROIs'...
%     'G:\Sam\Data\SIM\180524 testing\gephryn488_cellfillmch_basoon647_002_synROIs'};
% squashpath = 'G:\Sam\Data\SIM\180524 testing\squashtest';
% csvpath = 'G:\Sam\Data\SIM\180524 testing\Squash\000_objFilterList.csv';

%% promts to select the correct folders and files etc
if isempty(synfolder) || isempty(squashpath) || isempty(csvpath)
%---- select the '*synROIs' folders to analyze
answer1 = questdlg('Select all the *synROIs folders you want to make analysis folders for.',...
    'Options', 'OK', 'Oops!','OK');
synfolder = uipickfiles('Prompt','Select all the synROIs folders you want to analyze');
%---- select the 'SquashMe' path to save the resulting rois
answer2 = questdlg('Select the Squash Folder to save your results, or make a new one',...
    'Options', 'OK', 'Oops!','OK');
squashpath = uigetdir('.','Select Your Squash Path or Create One');
%---- select the full path to the '*.csv' file to copy into each roi result folder
answer3 = questdlg('Select a .csv file to use a template inside each folder. Hint: normally something like **_objFilterList.csv',...
    'Options', 'OK', 'Oops!','OK');
[csvfile, csvfolder] = uigetfile('*.*','Select a .csv file to use a template inside each folder');
csvpath = fullfile(csvfolder,csvfile);
end

%% determine what the last roi analysis file number was inside the squash path and increment by 1
temp = dir(fullfile(squashpath));
if numel(temp)>3
    temp = temp(arrayfun(@(x) x.isdir,temp)); %find all directories
    temp = temp(3:end); %remove . and .. (current and parent directories from list)
    startnum = str2double(temp(end).name(1:3))+1; clear temp; %convert string to number and increment
else
    startnum = 1;
end
num = startnum; %create index to increment

%% Loop through all the synROIs folders
for ss = 1:numel(synfolder)
    currsynpath = synfolder{ss};
    [filepath,base_name,ext] = fileparts(currsynpath);
    % count how many files have 'roi' in them, divide by 4 to get how many roi
    % there are
    roifiles = dir(fullfile(synfolder{ss}, '*roi*'));
    n_folder = numel(roifiles)/4;  % how many individual rois were chosen for this image
    % Loop through and create a directory for each roi inside the current synROIs folder
    % note: each folder will have a unique number
    for ii = 1:n_folder
        % make the new folder for the file/roi inside the current synROIs folder
        newfolder_name = (sprintf(['%03d_', base_name, '_roi%d'],num,ii));
        currfold = fullfile(squashpath,newfolder_name);
        mkdir(currfold);
        % copy the csv file into the new folder, with the unique number
        % added to the name 
        csvnew_name = (sprintf('%03d_objFilterList.csv',num));
        copyfile(csvpath,fullfile(squashpath,newfolder_name,csvnew_name));
        % find all the roi associate files for the give file/roi/base_name
        n_roi = num2str(ii, '*roi%d_*');
        files = dir(fullfile(currsynpath, n_roi));
        % loop through found roi files and move them
        for j = 1:size(files)
            copyfile(fullfile(currsynpath,files(j).name), currfold);
        end
        num = num+1;
    end
end

