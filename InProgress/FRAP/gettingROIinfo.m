
%load datafile
dataim = GeneralAnalysis.loadtiff_1ch('G:\FromMicroscopeComputer\190117 mScGeph FRAP\mScarGeph_GephIB_cell3_FRAPtimeseries\561allfiles.tif');
%load text file with meta info
metafile = 'G:\FromMicroscopeComputer\190117 mScGeph FRAP\mScarGeph_GephIB_cell3_FRAPtimeseries\fullseries\mScarGeph_GephIB_cell3_20190117_100729 PM\mScarGeph_GephIB_cell3.txt';

% use Ashley's script to identify regions
fid = fopen(metafile);
tline = fgetl(fid);
idx = 1;
frappa = zeros(11,3);
while ischar(tline)
    
    % search for _FRAPPA lines
    re = regexp(tline,'(?<=_FRAPPA\tPoint\t NumberOfPoints\( \d\) : \( )\d+, \d+','match');
    
    if ~isempty(re) % if a match
        lb = regexp(tline,'(?<=Label\( )\d','match'); % get label number
        
        frappa(idx,:) = [str2num(re{:}) str2num(lb{:})];
        idx = idx+1;
    end
    
    tline = fgetl(fid);  
end
fclose(fid);

% Cleanup FRAPPA
max_label = max(frappa(:,3));
frappa(max_label+1:end,:) = []; 
frappa(:,3)=[];

% test to see if locations are correct
dataim
xsz = size(dataim,2);
ysz = size(dataim,1);
for ii = 1:size(frappa,1)
 rect = rectangle('Position',[(frappa(ii,1))-6,(frappa(ii,2))-6,10,10],...
                    'EdgeColor','red',...
                    'LineWidth',1);
end


