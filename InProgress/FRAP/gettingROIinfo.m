
%load datafile
pre_shift = GeneralAnalysis.loadtiff_1ch('G:\FromMicroscopeComputer\190303 mScarGeph FRAP_FingR\GephIB\cell1highexp_mScarGeph_GephIBGFP_20190303_70622 PM\cell1_mScarGeph_GephIBGFP_w0001_z.tiff');
% pre_shift = GeneralAnalysis.loadtiff_1ch('G:\FromMicroscopeComputer\190117 mScGeph FRAP\mScarGeph_GephIBOlig_cell1_FRAPtimeseries\561allseries.tif');
[dataim,sv_arr] = GeneralAnalysis.timedriftCorrect(pre_shift);
dataim = pre_shift;
%load text file with meta info
metafile = 'G:\FromMicroscopeComputer\190303 mScarGeph FRAP_FingR\GephIB\cell1highexp_mScarGeph_GephIBGFP_20190303_70622 PM\cell1_mScarGeph_GephIBGFP.txt';
% metafile = 'G:\FromMicroscopeComputer\190117 mScGeph FRAP\mScarGeph_GephIBOlig_cell1_FRAPtimeseries\fullseries\mScarGeph_GephIBOLIG_cell1_20190117_100140 PM\mScarGeph_GephIBOLIG_cell1.txt';
% use Ashley's script to identify regions
fid = fopen(metafile);
tline = fgetl(fid);
idx = 1;
frappa = [];
while ischar(tline)
    
    % search for _FRAPPA lines
    re = regexp(tline,'(?<=_FRAPPA\tPoint\t NumberOfPoints\( \d\) : \( )\d+, \d+','match');
    
    if ~isempty(re) % if a match
        lb = regexp(tline,'(?<=Label\( )\d+','match'); % get label number
        
        frappa = [frappa ; str2num(re{:}) str2num(lb{:})];
        idx = idx+1;
    end
    
    tline = fgetl(fid);  
end
fclose(fid);

% Cleanup FRAPPA
[~,testid]=sort(frappa(:,3),'ascend');
frappa = frappa(testid,:);
% max_label = max(frappa(:,3));
% frappa(max_label+1:end,:) = []; 
% frappa(:,3)=[];

% test to see if locations are correct
%% get raw sum ROI values
regionsums = [];
boxsize = 10;
for ii = 1:size(frappa,1)
    regions{ii} = dataim(frappa(ii,1)-ceil(boxsize/2):frappa(ii,1)+ceil(boxsize/2),frappa(ii,2)-ceil(boxsize/2):frappa(ii,2)+ceil(boxsize/2),:);
    regionsums(ii,:)  = squeeze(single(sum(regions{ii},[],[1 2])))';
end
%% normalize 0-1
h = dipshow(dataim);
for ii = 1:size(frappa,1)
 rect = rectangle('Position',[(frappa(ii,1))-ceil(boxsize/2),(frappa(ii,2))-ceil(boxsize/2),boxsize,boxsize],...
                    'EdgeColor','red',...
                    'LineWidth',1);
end
startvals = [];
startvals = mean(regionsums(:,1:6),2);


uiwait(msgbox('Select a representative unbleached region','Title','modal'));
[B2,C2] = dipcrop(h);
% show bleached region
rect_unbleached = rectangle('Position',[C2(1,1),C2(1,2),C2(2,1),C2(2,2)],...
                    'EdgeColor','blue',...
                    'LineWidth',1);
uiwait(msgbox('Select a representative background region','Title','modal'));
[B1,C1] = dipcrop(h);
rect_background = rectangle('Position',[C1(1,1),C1(1,2),C1(2,1),C1(2,2)],...
                    'EdgeColor','green',...
                    'LineWidth',1);
          


background = squeeze(single(sum(B1,[],[1 2])));
unbleachROI = squeeze(single(sum(B2,[],[1 2])));
tonormvals = unbleachROI./mean(unbleachROI(1:3));
test = (regionsums'./tonormvals);
startvals = mean(test(1:3,:));

normvals = (test./startvals);
figure; plot(normvals)


