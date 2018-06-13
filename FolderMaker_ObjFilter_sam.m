
% add your Squash analysis path or leave blank to be promted to create a new one.
% add the csvfile name and path that will be copied into each roi folder

squashpath = '';%make a SquasshMe folder
csvfilebase = '000_objFilterList.csv';
csvpath = squashpath;
%% Get The Paths if Not Already Set Above
if isempty(squashpath)
    squashpath = uigetdir('.','Select Your Squash Path or Create One');
end
if isempty(csvpath)
    [csvfile, csvpath] = uigetfile(squashpath,['Find your ',csvfilebase,' file'],'*.*');
    copyfile(fullfile(csvpath,csvfile),fullfile(squashpath,csvfile));
else 
    searchfiles = dir(fullfile(squashpath,csvfilebase));
    if isempty(searchfiles)
        msgbox('Need to have a csvfilebase in your csvpath!','Oops');
    end
end
%% Have User Select Folders to Analyze
move_from = uigetdir(fullfile(squashpath,'..'),'Select Folder to Analyze ROIs'); %'C:\Users\Sara\Desktop\Test\con10cell1006_glu_synROIs\';
base_name = strsplit(move_from, '\'); % base_name = 'con10cell1006';
base_name = base_name{end};
% count how many files have 'roi' in them, divide by 4 to get how many roi
% there are
roifiles = dir(fullfile(move_from, '*roi*'));
n_folder = numel(roifiles)/4;

%% determine what the last roi analysis file number was and increment by 1
temp = dir(fullfile(squashpath));
if numel(temp)>3
    temp = temp(arrayfun(@(x) x.isdir,temp)); %find all directories
    temp = temp(3:end); %remove . and .. (current and parent directories from list)
    startnum = str2double(temp(end).name(1:3))+1; clear temp; %convert string to number and increment
else
    startnum = 1;
end
num = startnum; %create index to increment
%% Loop through and create the directories
for i = 1:n_folder  
    folder_name = (sprintf(['%03d_', base_name, '_roi%d'],num,i));
    mkdir(fullfile(squashpath,folder_name));
    new_name = (sprintf('%03d_objFilterList.csv',num));
    copyfile(fullfile(squashpath,csvfilebase),fullfile(squashpath,folder_name,new_name));
    % find all the rois for a given base_name
    n_roi = num2str(i, '*roi%d_*');
    files = dir(fullfile(move_from, n_roi));
    % loop through found roi files and move them
    for j = 1:size(files)
        copyfile(fullfile(move_from,files(j).name), fullfile(squashpath,folder_name));
    end
    num = num+1;
end

