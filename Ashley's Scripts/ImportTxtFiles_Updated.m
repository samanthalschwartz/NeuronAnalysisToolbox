%% Get File Path
[f,p] = uigetfile('*.txt')
file_path = fullfile(p,f);
 
%%
fid = fopen(file_path);
tline = fgetl(fid);
 
idx = 1
frappa = zeros(11,3);
while ischar(tline)
    
    % search for _FRAPPA lines
    re = regexp(tline,'(?<=_FRAPPA\tPoint\t NumberOfPoints\( \d\) : \( )\d+, \d+','match')
    
    if ~isempty(re) % if a match
        lb = regexp(tline,'(?<=Label\( )\d','match'); % get label number
        
        frappa(idx,:) = [str2num(re{:}) str2num(lb{:})];
        idx = idx+1;
    end
    
    tline = fgetl(fid)    
end
fclose(fid);

% Cleanup FRAPPA

max_label = max(frappa(:,3));
frappa(max_label+1:end,:) = []; 
frappa(:,3)=[];


%Code for a GUI that will ask me how many FRAP points
%and then I put in my code (lines 11 and 21) the variable that the GUI
%gives to me instead of 4 in this case.

% file_no=inputdlg('How many files?');

% file_no=str2double(file_no{1});
