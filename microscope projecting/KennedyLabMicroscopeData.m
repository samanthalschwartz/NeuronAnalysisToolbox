classdef KennedyLabMicroscopeData < handle
    properties
        test = [];
    end
    
    methods (Static)
        function moveMetaMorphThumbs(startdir,idstring)
            if nargin<2
                idstring = 'thumbs';
            end
            if nargin<1
                startdir = uigetdir(pwd,'Select a directory to move all ''thumb'' files to a new folder');
            end
            thumbdir = fullfile(startdir,idstring);
            if ~exist(thumbdir,'dir')
                mkdir(thumbdir)
            end
            allfiles = dir(fullfile(startdir));
            bools = cell2mat(arrayfun(@(x) (~isempty(strfind(x.name,idstring)) & ~x.isdir),allfiles,'UniformOutput',false));
            movefiles = allfiles(bools);
            for mm = 1:size(movefiles,1)
                movefile(fullfile(startdir,movefiles(mm).name),thumbdir);
            end
        end
        
        function boolarray = findfiles(allfiles,inputstr)
           boolarray = cell2mat(arrayfun(@(x) (~isempty(strfind(x.name,inputstr)) & ~x.isdir),allfiles,'UniformOutput',false));   
        end
    end
    
    
end
%%
% loading files
% first have using click a file to do something to.

% have user identify the filename
startdir = 'G:\FromMicroscopeComputer\Sam\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_10uL-1';
teststr = 'NL1-GFP_oldzap500nM_1_+488_ABbaseline_w3640_s1_t1';


% put all times together as time series

% first test the first variable
out = strfind(teststr,'_');
lett=teststr(out(end-1)+1);
switch lett
    case 't'
        % then there is a time series
        timestr = out(end+1);        
        % put together time series
        allfiles = dir(fullfile(startdir));
        boolarray = findfiles(teststr(1:timestr));
    case 's'
        % then there are multiple positions
        % no time series
        pos = out(end);
    case 'w'
        % multiple channels no positions or time
        ch = out(end)

%identify _w with a number 
% identify _s with a number
% identify _t with a number 

% find _w followed by number