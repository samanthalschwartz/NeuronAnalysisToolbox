%% script to calculate distance for each object
ch_dhfr = 'F:\Sam\050118 NL1 insertion\cell3-C3.tif';
ch_cellfill = 'F:\Sam\050118 NL1 insertion\cell3-C2.tif';

% first create cell mask based on cell fill channel
cellmaskbefore = bfopen(ch_cellfill);
cellmask = dip_image(zeros([size(cellmaskbefore{1,1}{1}),size(cellmaskbefore{1,1},1)]));
for ii = 1:size(cellmask,3)
 cellmask(:,:,ii-1) = cellmaskbefore{1,1}{ii};
end
% remove any unwanted pixels 
cropx = 15:size(cellmask,1)-15;
cropy = 15:size(cellmask,2)-15;
cropt = 0:size(cellmask,3)-1;
cellmask = cellmask(cropx,cropy,cropt);
cellmask0 = gaussf(cellmask,[1 1 1]);
cellmask1 = dxx(cellmask0) + dyy(cellmask0);
cellmask2 = -cellmask1;
cellmask2(cellmask2<0) = 0;
threshval2 = multithresh(single(cellmask2),2);
threshval1= multithresh(single(cellmask0),2);
[out1,thres1] = threshold(cellmask,'fixed',threshval1(1));
[out,thres] = threshold(cellmask2,'fixed',threshval2(1));
cm = out1 | out; %-- this is the cellmask now
%now have user select soma mask to compare distance from:
h = dipshow(cellmask1(:,:,10));
diptruesize(h,200);
[soma, v] = diproi(h);
close(h)

%% now loop through and label aggregates in each frame and then get distances from selected soma
dhfrbefore = bfopen(ch_dhfr);
cell1 = dip_image(zeros([size(dhfrbefore{1,1}{1}),size(dhfrbefore{1,1},1)]));
for ii = 1:size(cell1,3)
 cell1(:,:,ii-1) = dhfrbefore{1,1}{ii};
end
cell1 = cell1(cropx,cropy,cropt);



gcell1 = gaussf(cell1,[1 1 0]);
ddcell1 = dxx(gcell1) + dyy(gcell1);
ddcell2 = -ddcell1;
ddcell2(ddcell2<0) = 0;
threshval = multithresh(single(ddcell2),2);
[dhfrout,dhfrthres] = threshold(ddcell2,'fixed',threshval(1));

distVtimemat =[];
wb = waitbar(0,'Analyzing Distances...');
bw2 = soma;
for ii = 0:(size(dhfrout,3)-1)
   
    cmframe = closing(cm(:,:,ii));
    temp = dhfrout(:,:,ii)*cmframe;
     dhfr_labeled = label(temp,1);
    D2 = bwdistgeodesic(logical(cmframe), logical(bw2), 'quasi-euclidean');
    if ~any(dhfr_labeled>0)
        continue;
    else
        labeledmat = zeros(max(dhfr_labeled),2);
        labeledmat(:,1) = ii;
        for ll = 1:max(dhfr_labeled)
            bw1 = squeeze((dhfr_labeled == ll));
            D1 = bwdistgeodesic(logical(cmframe), logical(bw1), 'quasi-euclidean');
            D = D1+D2;
            D = round(D * 8) / 8;
            D(isnan(D)) = inf;
            paths = imregionalmin(D);
            paths_thinned_many = bwmorph(paths, 'thin', inf);
% %             %-- to visualize path
            figure;
            P = false(size(logical(D)));
            P = imoverlay(P, ~logical(cmframe), [1 0 0]);
            P = imoverlay(P, paths, [.5 .5 .5]);
            P = imoverlay(P, paths_thinned_many, [1 1 0]);
            P = imoverlay(P, logical(dhfr_labeled), [1 1 1]);
            imshow(P, 'InitialMagnification', 'fit')
% %             %---
            dist = size(find(paths_thinned_many),1);
            labeledmat(ll,2) = dist;
        end
        distVtimemat = [distVtimemat;labeledmat];
    end
    waitbar(ii/(size(dhfrout,3)-1),wb);
end
close(wb);
save(fullfile([ch_dhfr(1:end-4) '_distances']),'distVtimemat');

%% plot results
ff= figure;
hist3(distVtimemat,'CDataMode','auto','Nbins',[max(distVtimemat(:,1))-2 30],'FaceColor','interp');
ff= figure;
toplot = distVtimemat;
toplot(toplot(:,2)>10,:) = [];

cnt = hist3(distVtimemat,'CDataMode','auto','FaceColor','interp');
hist3(distVtimemat,'CDataMode','auto','FaceColor','interp')
figure
hist3(toplot,'CDataMode','auto','FaceColor','interp')
xlabel('Time (Frames)');
ylabel('Distance from Soma (pixels)');
view(2)
figure; surf(cnt)