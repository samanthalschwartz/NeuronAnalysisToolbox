function filepaths = selectMultipleFiles(startdir, ext, message)
%% this function calls uigetfile for multiple directories and returns a cell array of full file paths
% stardir: input string for starting directory. '' means use
% default matlab
% ext: optional file extensions to look for (string) 
% message: optional message to include
% filepaths: cellarray of fullfile file paths
%%
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile(fullfile(startdir, ['*' ext]), message,'Multiselect','on');
if ~iscell(FILENAME)
    FILENAME = {FILENAME};
end
filepaths = cellmap(@(x) fullfile(PATHNAME,x),FILENAME)';

stp = 1;
while stp == 1
    choice = questdlg('Choose more files?', ...
        'File Picker', ...
        'Yes','Nope. All done','Nope. All done');
    switch choice
        case 'Yes'
            [FILENAME, PATHNAME, FILTERINDEX] = uigetfile(fullfile(startdir, ['*' ext]), message,'Multiselect','on');
            if ~iscell(FILENAME)
                FILENAME = {FILENAME};
            end
            newfilepaths = cellmap(@(x) fullfile(PATHNAME,x),FILENAME)';
            filepaths = [filepaths;newfilepaths];
        case 'Nope. All done'
            stp = 0;
    end
end
end