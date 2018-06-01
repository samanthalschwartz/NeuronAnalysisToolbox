%% do this correctly
% make some colormaps
blackjet = jet;
blackjet = flip(blackjet);
blackjet(1,:) = [0 0 0]; blackjet(end,:) = [0 0 0];
hsvblack = hsv;
hsvblack(end,:) = [0 0 0];
hotblack = hot;
hotblack(end,:) = [0 0 0];


% cell1 = readtimeseries('G:\Sam\050118 NL1 insertion\C3-cell3.tif','',[],1,0);
% cell1before = bfopen('G:\Sam\050118 NL1 insertion\C3-cell3.tif');
cell1before = bfopen('G:\Sam\050118 NL1 insertion\cell4-C3.tif');
cell1 = dip_image(zeros([size(cell1before{1,1}{1}),size(cell1before{1,1},1)]));
for ii = 1:size(cell1,3)
 cell1(:,:,ii-1) = cell1before{1,1}{ii};
end


movierange = [0:37,39:size(cell1,3)-1];
gcell1 = unif(cell1(:,:,movierange),[1 1 0]);
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
cellmaskbefore = bfopen('G:\Sam\050118 NL1 insertion\cell4-C2.tif');
cellmask = dip_image(zeros([size(cellmaskbefore{1,1}{1}),size(cellmaskbefore{1,1},1)]));
for ii = 1:size(cellmask,3)
 cellmask(:,:,ii-1) = cellmaskbefore{1,1}{ii};
end

cellmask = cellmask(:,:,movierange);
% cellmask2 = gaussf(cellmask,[2 2 1])
cellmask2 = dxx(cellmask) + dyy(cellmask);
cellmask2 = -cellmask2;
cellmask3 = gaussf(cellmask2);
cellmask3(cellmask3<0) = 0;

seeds = minima(cellmask2(:,:,10),1,0);
      image_out = waterseed(seeds,cellmask2(:,:,10),1);
      
image_out = watershed(cellmask2(:,:,10),1,0);

testval = dip_image(logical(cellmask3>20));
testval = closing(closing(testval));

incellstuff = labeled_im*testval;
test = min(incellstuff,incellstuff>0,3);
test(test>(size(out,3)+1)) = 0;





% now use cell mask to only include objects that are within the cell
% cellmask = readtimeseries('G:\Sam\050118 NL1 insertion\C2-cell3.tif','',[],1,0);



h = dipshow(test(15:end-15,15:end-15,:),blackjet);
dipmapping(h,[0 size(out,3)]);
diptruesize(h,100)
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

%% 
cell1 = bfopen('G:\Sam\050118 NL1 insertion\C3-cell3.tif');
cell1di = dip_image(zeros([size(cell1{1,1}{1}),size(cell1{1,1},1)]));
for ii = 1:size(cell1di,3)
 cell1di(:,:,ii-1) = cell1{1,1}{ii};
end


%%
cellmask = readtimeseries('G:\Sam\050118 NL1 insertion\C2-cell1.tif','',[],1,0)
% cellmask2 = gaussf(cellmask,[2 2 1])
cellmask2 = dxx(cellmask) + dyy(cellmask);
cellmask2 = -cellmask2;
cellmask3 = gaussf(cellmask2)
cellmask3(cellmask3<0) = 0;
testval = dip_image(logical(cellmask3>20));
overlay(cellmask,testval,hot)


bwim1=adaptivethreshold(single(cellmask(:,:,10)),11,0.03,1);
dipshow(bwim1)

threshval2 = multithresh(single(cellmask3),3);
[out,threshval2] = threshold(single(cellmask),'fixed',threshval2(end));
overlay(cellmask,out)

image_out = watershed(cellmask2(:,:,20),1,0,0)
%%
h = dipshow(test*testval,blackjet);
dipmapping(h,[0 size(out,3)])

