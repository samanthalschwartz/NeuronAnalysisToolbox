topdir = 'Z:\Sam\Data!\181121 Halo-pDisplay Kinetics';
topdir = 'Z:\Sam\Data!\Other Datas\180516_GCaMP6';
namestr = 'GFPcellfill_Rabl11mCh_TfR-Halo_timeseries_dyefr2_nowash2_w1488_t4';
namestr = 'GCaMP_postBayKcell1_t1.TIF';



test=regexp(namestr,'_w*\d*(_s|_t)');

infostr = namestr(test:end);

test2=regexp(infostr,'_s\d');
infostr(test2:end)


% first test if there are defined channels
test1=regexp(namestr,'_w*\d*(_s|_t|[?.tif])');

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
