%% do this correctly
% make some colormaps
blackjet = flip(jet(255));
blackjet(1,:) = [0 0 0]; blackjet(end,:) = [0 0 0];


ch_dhfr = 'F:\Sam\050118 NL1 insertion\cell4-C3.tif';
ch_cellfill = 'F:\Sam\050118 NL1 insertion\cell4-C2.tif';


% cell1 = readtimeseries('G:\Sam\050118 NL1 insertion\cell3-C3.tif','',[],1,0);
% cell1before = bfopen('G:\Sam\050118 NL1 insertion\C3-cell3.tif');
cell1before = bfopen(ch_dhfr);
cell1 = dip_image(zeros([size(cell1before{1,1}{1}),size(cell1before{1,1},1)]));
for ii = 1:size(cell1,3)
 cell1(:,:,ii-1) = cell1before{1,1}{ii};
end

cropx = 15:size(cell1,1)-15;
cropy = 15:size(cell1,2)-15;
cropt = 0:size(cell1,3)-1;

gcell1 = gaussf(cell1(cropx,cropy,cropt),[1 1 0]);
ddcell1 = dxx(gcell1) + dyy(gcell1);
ddcell2 = -ddcell1;
ddcell2(ddcell2<0) = 0;
threshval = multithresh(single(ddcell2),2);
[out,thres] = threshold(ddcell2,'fixed',threshval(1));

% loop through each frame and value label by frame number
labeled_im = dip_image(zeros(size(out)));
for ii = 1:size(out,3)
    temp = label(out(:,:,ii-1),1);
    temp(temp>0) = ii;
    labeled_im(:,:,ii-1) = temp;
end

% cellmaskbefore = bfopen('G:\Sam\050118 NL1 insertion\C2-cell3.tif');
cellmaskbefore = bfopen(ch_cellfill);
cellmask = dip_image(zeros([size(cellmaskbefore{1,1}{1}),size(cellmaskbefore{1,1},1)]));
for ii = 1:size(cellmask,3)
 cellmask(:,:,ii-1) = cellmaskbefore{1,1}{ii};
end

cellmask = cellmask(cropx,cropy,cropt);
cellmask2 = gaussf(cellmask,[1 1 0]);
cellmask2 = dxx(cellmask2) + dyy(cellmask2);
cellmask2 = -cellmask2;
cellmask2(cellmask2<0) = 0;

threshval1 = multithresh(single(cellmask2),2);
[cm1,thres1] = threshold(cellmask2,'fixed',threshval1(1));
threshval2 = multithresh(single(cellmask),2);
[cm2,thres2] = threshold(cellmask,'fixed',threshval2(1));
cellmask2use = cm1|cm2;
testval = closing(closing(cellmask2use));

incellstuff = labeled_im*testval;
test = min(incellstuff,incellstuff>0,3);
test(test>(size(out,3)+1)) = 0;

% now use cell mask to only include objects that are within the cell
% cellmask = readtimeseries('G:\Sam\050118 NL1 insertion\C2-cell3.tif','',[],1,0);


h = dipshow(test(15:end-15,15:end-15,:),blackjet);
dipmapping(h,[0 size(out,3)]);
diptruesize(h,100);
% get colorbar tick info
colorunit = size(out,3)/255;
numofcolbarval = 4;
colbarplace =[0:numofcolbarval]*255/4;
colbarval = floor(colbarplace * colorunit);
c = colorbar;
c.Location = 'WestOutside';
c.Ticks = colbarplace;
c.TickLabels = colbarval;
c.FontSize = 16;
c.Label.String = 'First Frame with Intensity';
c.Label.FontSize = 16;
h.OuterPosition = h.OuterPosition + [0 0 400 50];

cellmasksum = sum(cellmask2use,[],3);
cellmasksum(cellmasksum>0)= 1;
B = bwperim(logical(cellmasksum));
dipshow(B)
B = bwboundaries(logical(cellmasksum));
roipoly(B)
saveas(h,fullfile([ch_dhfr(1:end-4) '_minframeImage']),'png');
saveas(h,fullfile([ch_dhfr(1:end-4) '_minframeImage']),'fig');
close(h);

