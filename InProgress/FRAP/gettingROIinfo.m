clear all; close all;
%% load datafile and drift correct
datafolder = 'G:\FromMicroscopeComputer\190410 pHujiGaba\GephIB\gephIB_pHujiGabaFRAP_cell2_20190410_120921 PM';
fileinfo = fullfile(datafolder,'gephIB_pHujiGabaFRAP_cell2_w0001_z.tiff');
metafile = fullfile(datafolder,'gephIB_pHujiGabaFRAP_cell2.txt');
pre_shift = GeneralAnalysis.loadtiff_1ch(fileinfo);
[dataim,sv_arr] = GeneralAnalysis.timedriftCorrect(pre_shift);
GeneralAnalysis.LibTiff(dataim,[fileinfo(1:end-5) 'Shift' fileinfo(end-4:end)]);
bleachfr = 5;
boxsize = 6;
%% use Ashley's script to identify regions
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

%% select the appropriate background and unbleached regions for normalization
h = dipshow(dataim,'log');
for ii = 1:size(frappa,1)
 rect = rectangle('Position',[(frappa(ii,1))-ceil(boxsize/2),(frappa(ii,2))-ceil(boxsize/2),boxsize,boxsize],...
                    'EdgeColor','red',...
                    'LineWidth',1);
end
uiwait(msgbox('Select representative unbleached regions, close window when finished','Title','modal'));
% loop through and make as many unbleached regions as you want
ub_boxes = {};
ub_rois = {};
while(ishandle(h))
    try
        [B2,C2] = dipcrop(h);
    catch
        break;
    end
    rect_unbleached = rectangle('Position',[C2(1,1),C2(1,2),C2(2,1),C2(2,2)],...
        'EdgeColor','blue',...
        'LineWidth',1);
    if isempty(ub_boxes)
        ub_boxes = {B2};
        ub_rois = {C2};
    else
        ub_boxes = [ub_boxes, {B2}];
        ub_rois = [ub_rois,{C2}];
    end
end
% replot everything
h = dipshow(dataim,'percentile');
for ii = 1:size(frappa,1)
    rect = rectangle('Position',[(frappa(ii,1))-ceil(boxsize/2),(frappa(ii,2))-ceil(boxsize/2),boxsize,boxsize],...
        'EdgeColor','red',...
        'LineWidth',1);
end
for ii = 1:numel(ub_rois)
   rect_unbleached = rectangle('Position',[ub_rois{ii}(1,1),ub_rois{ii}(1,2),ub_rois{ii}(2,1),ub_rois{ii}(2,2)],...
        'EdgeColor','blue',...
        'LineWidth',1); 
end

% show bleached region
uiwait(msgbox('Select a representative background region','Title','modal'));
[B1,C1] = dipcrop(h);
rect_background = rectangle('Position',[C1(1,1),C1(1,2),C1(2,1),C1(2,2)],...
                    'EdgeColor','green',...
                    'LineWidth',1);
%% save all regions
%-- frap regions
for ii = 1:size(frappa,1)
    x(1) = frappa(ii,1)-ceil(boxsize/2);
    y(1) = frappa(ii,2)+ceil(boxsize/2);
    
    x(2) = frappa(ii,1)-ceil(boxsize/2);
    y(2) = frappa(ii,2)-ceil(boxsize/2);
    
    x(3) = frappa(ii,1)+ceil(boxsize/2);
    y(3) = frappa(ii,2)-ceil(boxsize/2);
    
    x(4) = frappa(ii,1)+ceil(boxsize/2);
    y(4) = frappa(ii,2)+ceil(boxsize/2);
    
    roidir = fullfile(datafolder,'BleachingROIs');
    if ~exist(roidir)
        mkdir(roidir);
    end
roifilename = fullfile(roidir,['BleachingROI_' num2str(ii) '.csv']);
csvwrite(roifilename,[x',y']);
end
%-- non bleach regions
for ii = 1:numel(ub_rois)
   x(1) = ub_rois{ii}(1,1);%   ceil(ub_rois{ii}(2,1)/2);
   y(1) = ub_rois{ii}(1,2)+ub_rois{ii}(2,2);
    
    x(2) = ub_rois{ii}(1,1);%-ceil(ub_rois{ii}(2,1)/2);
    y(2) = ub_rois{ii}(1,2)%-ceil(ub_rois{ii}(2,2)/2);
    
    x(3) = ub_rois{ii}(1,1)+ub_rois{ii}(2,1);
    y(3) = ub_rois{ii}(1,2);%-ceil(ub_rois{ii}(2,2)/2);
    
    x(4) = ub_rois{ii}(1,1)+ub_rois{ii}(2,1);
    y(4) = ub_rois{ii}(1,2)+ub_rois{ii}(2,2);
    
    roidir = fullfile(datafolder,'NonbleachingROIs');
    if ~exist(roidir)
        mkdir(roidir);
    end
    roifilename = fullfile(roidir,['NonbleachingROIs_' num2str(ii) '.csv']);
csvwrite(roifilename,[x',y']);
end
%-- background region
x(1) = C1(1,1);%   ceil(ub_rois{ii}(2,1)/2);
y(1) = C1(1,2)+C1(2,2);

x(2) = C1(1,1);%-ceil(ub_rois{ii}(2,1)/2);
y(2) = C1(1,2);%-ceil(ub_rois{ii}(2,2)/2);

x(3) = C1(1,1)+C1(2,1);
y(3) = C1(1,2);%-ceil(ub_rois{ii}(2,2)/2);

x(4) = C1(1,1)+C1(2,1);
y(4) = C1(1,2)+C1(2,2);

roidir = fullfile(datafolder,'Background');
if ~exist(roidir)
    mkdir(roidir);
end
roifilename = fullfile(roidir,'Background.csv');
csvwrite(roifilename,[x',y']);
%% subtract background from ROI values and unbleached regions
regions = {};
regionsums = [];
sumbackbox = squeeze(single(sum(B1,[],[1 2])));
ppx_background = mode(sumbackbox)/(C1(2,1)*C1(2,2));
for ii = 1:size(frappa,1)
    tempbox = dataim(frappa(ii,1)-ceil(boxsize/2):frappa(ii,1)+ceil(boxsize/2),frappa(ii,2)-ceil(boxsize/2):frappa(ii,2)+ceil(boxsize/2),:);
    regions{ii} = tempbox - ppx_background;
    regionsums(ii,:)  = squeeze(single(sum(regions{ii},[],[1 2])))';
end
unbleachROIsums = [];
ub_regions = {};
for ii = 1:numel(ub_boxes)
    ub_regions{ii} = ub_boxes{ii} - ppx_background;
    unbleachROIsums(:,ii) = squeeze(single(sum(ub_regions{ii},[],[1 2])));
end
%% normalize and plot
% --- get info about unbleached regions
unbleachROIn = unbleachROIsums./unbleachROIsums(1,:);
figure; plot(unbleachROIn); hold on;
MunbleachROIn = mean(unbleachROIn,2);
sMunbleachROIn = movmean(MunbleachROIn,4);
plot(sMunbleachROIn,'--k','LineWidth',2);
xlabel('Frames');
ylabel('Smoothed Average Intensity of non-bleached Regions');
% -- get initial values to normalize to 1
unbleachROIsums = sMunbleachROIn;
tonormvals = unbleachROIsums./mean(unbleachROIsums(1:(bleachfr-1)));
blcorrected_regionsum = (regionsums'./tonormvals);

% --- create 2 different recovery traces: 
% not normalized to bleaching amount
test_1only = blcorrected_regionsum;
startvals_1only = mean(test_1only(1:(bleachfr-1),:));
normvals_1only = (test_1only./startvals_1only);
figure; plot(normvals_1only);
xlswrite(fullfile(datafolder,'Results'),normvals_1only,'1onlyNormalized')

% for data trace subtract off bleached value to set to 0.
test = blcorrected_regionsum - blcorrected_regionsum(bleachfr,:);
startvals = mean(test(1:(bleachfr-1),:));
normvals = (test./startvals);
figure; plot(normvals);
saveas(gcf,fullfile(datafolder,'Results_1NormOnly'));
close all;
% ----
xlswrite(fullfile(datafolder,'Results'),normvals,'01Normalized');
saveas(gcf,fullfile(datafolder,'Results_01Norm'));
close all;
save(fullfile(datafolder,'Results'));

