topdir = 'Z:\Sam\Data!\181115 HaloPdisplay';
% topdir = 'Z:\Sam\Data!\Other Datas\180516_GCaMP6';
namestr = 'TfR-Halo_mChCellFill_6nMhalo660_longExp_timecourse_w1561_s4_t6';
% namestr = 'GCaMP_postBayKcell1_t1.TIF';

chstr = '_w';
timestr = '_t';
stagestr = '_s';

% find if there are any channels, stage positions, timeseries
[matches] = regexp(namestr,'_w\d|_t\d|_s\d','match');
[out] = regexp(namestr,'_w\d|_t\d|_s\d');

savefoldername = namestr(1:min(out)-1);

% if there are multiple times, combine all into one image and move to
% different folder
timeid = cellfun(@(x) ~isempty(strfind(x,timestr)),matches);
timepos = out(timeid);
currname = namestr(1:timepos+1);
files = dir2cell(topdir,[currname '*']);







%% --- old to delete
test=regexp(namestr,'_w*\d*(_s|_t)');

infostr = namestr(test:end);

test2=regexp(infostr,'_s\d');
infostr(test2:end)


% first test if there are defined channels
test1=regexp(namestr,'_w*\d*(_s|_t|[?.tif])','split');


% if there are not multiple channels then check if it is a time series
if isempty(test)
    % test if there are multiple timepoint
    % 1) no channels, multiple time points
    test = regexp(namestr,'_t\d*[?.tif]');
else
    % 2) no channels, no time points
end

%  if there are defined channels, test if there are multiple stage positions
test=regexp(namestr,'_w*\d*_s');
if isempty(test)
   
     test = regexp(namestr,'_t\d*[?.tif]');
      % 3) multiple channels, not multiple stage positions, test if there are multiple time points
else
    
end

% 5) if there are defined channels and multiple stage positions, test if there are multiple time points
