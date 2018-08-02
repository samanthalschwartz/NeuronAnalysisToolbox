uiopen('E:\Sam\Data\180516_GCaMP6\-BayK\cell1pre.tif',1);
uiopen('C:\Users\schwsama\Documents\Data\Hannah Calcium\Abeta1Post.tif',1)
% uiopen('E:\Sam\Data\180516_GCaMP6\+BayK\+BayK_cell1.tif',1)
gimage = gaussf(image,[2 2 2]);
diffimage = dzz(gimage);
adiff = abs(diffimage);
aa = GeneralAnalysis.imgThreshold_fixedUserInput(adiff)

tt = threshold(diffimage,'otsu');
tt = threshold(adiff,'otsu');
sadiff = sum(tt,[],3);
mask = sadiff>0;
lbl = label(mask,1,20,);

%%
msr = measure(lb,adiff,'P2A');
szs = msr.size;
bins 
h = histogram(szs,300)

% lblmov = repmat(lbl,1,1,size(gimage,3));
%% for each label, get intensity over the entire movie
trace = zeros(max(lbl),size(lblmov,3));
wb = waitbar(0,'Calculating...');
for ll = 1:max(lbl)%[1:8,10:max(lbl)]
%     tic
currmask = lbl==ll;
bcurrmask = bdilation(currmask);
mask2use = repmat(bcurrmask,1,1,size(image,3));
imgINmask = image*mask2use;
sumval = sum(imgINmask,[],[1 2]);
sizeval = sum(bcurrmask,[],[1,2]);
trace(ll,:) = sumval./sizeval;
% toc
waitbar(ll/max(lbl),wb);
end
close(wb)

allsum = sum(trace,2);
[~, ordx] = sort(allsum, 'ascend');
ord_trace = trace(ordx,:);
%% now plot
times = 1:size(ord_trace,2);
figure;
% wb = waitbar(0,'Plotting some things...');
% cnt = 0;
for ii = 1:size(ord_trace,1)
    %     waitbar(cnt/ii*jj,wb);
    for jj = 1:size(ord_trace,2)
        %         if vals{ii}(jj)
        %             FaceColor = 'k';
        %         else
        %             FaceColor = 'w';
        %         end
        p = patch([times(jj),times(jj)+1,times(jj)+1,times(jj)],[ii-1, ii-1, ii, ii],ord_trace(ii,jj));
        set(p,'FaceColor','flat','EdgeColor','none','CData',);
        %
        %         cnt = cnt+1;
    end
end
title('GCamp Intensity (AU)','FontSize',16)
xlabel('Frame','FontSize',16);
ylabel('Puncta #','FontSize',16)

%% old to delete
ss = threshold(sadiff,'otsu');

ww = watershed(-squeeze(adiff),2);
sadiff(ww) = 0;
mask = sadiff>0;
lbl = label(mask,1)

thresh = GeneralAnalysis.imgThreshold_fixedUserInput(adiff, sadiff);
sthresh = GeneralAnalysis.imgThreshold_fixedUserInput(sum(thresh,[],3));

I = single(sadiff);
se = strel('disk', 20);
Ie = imdilate(I, se);

Iobr = imreconstruct(Ie, I);
dipshow(Iobr)

fgm = imregionalmax(Iobrcbr);


D = dt(mask)
DL = watershed(squeeze(D));
bgm = DL == 0;
figure
imshow(bgm), title('Watershed ridge lines (bgm)')