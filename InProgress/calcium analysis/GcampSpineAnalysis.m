%% load in files and all that (cell fill and align etc)
ga = GeneralAnalysis(); %make just an empty version of the class to be able to use ga as a shortcut
% -- get pre file
[FILENAME_pre, PATHNAME_pre] = uigetfile(fullfile(pwd,'*.*'),'Select a pre condition calcium image');
prompt_pre = 'Select the Unique File Identifier: use * for wildcard';
name = 'Calcium Data Selection';
defaultanswer = {FILENAME_pre};
numlines = 1;
flstring_pre=inputdlg(prompt_pre,name,numlines,defaultanswer);
flstring_pre = flstring_pre{1};
% -- get post file
[FILENAME_post, PATHNAME_post] = uigetfile(fullfile(PATHNAME_pre,'*.*'),'Select a post condition calcium image');
prompt_post = 'Select the Unique File Identifier: use * for wildcard';
name = 'Calcium Data Selection';
defaultanswer = {FILENAME_post};
numlines = 1;
flstring_post=inputdlg(prompt_post,name,numlines,defaultanswer);
flstring_post = flstring_post{1};
% -- get cell fill file
[FILENAME_cf, PATHNAME_cf] = uigetfile(fullfile(PATHNAME_post,'*.*'),'Select the corresponding cell fill image');

img_pre = dip_image(ga.loadtiffseries(PATHNAME_pre,flstring_pre));
img_post = dip_image(ga.loadtiffseries(PATHNAME_post,flstring_post));
cellfill_pre = ga.loadtiff_1ch(fullfile(PATHNAME_cf,FILENAME_cf));

% align post image with last frame of pre-image
pre_lastframe = img_pre(:,:,end);
post_firstframe = img_post(:,:,end);
img_post_foralign = cat(3,pre_lastframe,post_firstframe);
[out,sv_arr] = ga.timedriftCorrect(img_post_foralign);
newimg_post = ga.applyshift2series(img_post,sv_arr);

fullseq = cat(3,img_pre,newimg_post);
 % align cellfill to img
imgsum = sum(fullseq,[],3);
corrimg = cat(3,sum(fullseq(:,:,end-5:end),[],3),cellfill_pre);
dc_corrimg = ga.timedriftCorrect(corrimg);
cellfill = dc_corrimg(:,:,end);
%% identify calcium
maxsz = 20;
tic;
ca = fullseq;
thresh_param = 1;
h = waitbar(0,'Identifying Calcium Transients (Finding Intensity Changes) ....');
% h = msgbox('Thresholding Image (Identifying Calcium Transients)....','Thanks');
imgg = gaussf(ca,[1 1 0]);
waitbar(.2,h);
diffimage = dzz(imgg);
waitbar(.3,h,'Identifying Calcium Transients (Smoothing Data)....');
adiff = abs(diffimage);
lp = ga.imgLaplaceCutoff(adiff,[1 1 0],[1 1 0]); %this helps smooth image (especially enhancing round objects)
waitbar(.4,h);
lp_ = lp.^thresh_param;
waitbar(.5,h,'Identifying Calcium Transients (Thresholding)....');
tt = threshold(lp_,'otsu');
lb = label(tt);
waitbar(.7,h,'Identifying Calcium Transients (Measuring)....');
msr = measure(tt,ca,{'size','DimensionsCube'});
spanval = max(msr.DimensionsCube([1 2],:));
% remove labels that are too big
largemasksbool = spanval>maxsz;
largemasksID = find(largemasksbool);
slb = single(lb);
for ii = 1:numel(largemasksID)
slb(slb==largemasksID(ii)) = 0;
end
cleanlb = dip_image(slb);
waitbar(1,h);
close(h)
toc;
% underimgin = repmat(cellfill,1, 1, size(fullseq(:,:,500:end),3));
% newmask =ga.cleanUpMask_manual_square(underimgin,cleanlb(:,:,500:end)>0);
% newcleanlb = cat(3,cleanlb(:,:,0:949),newmask);
% mask4overlay = bdilation(tt,2);
% [h,overlayim] = GeneralAnalysis.viewMaskOverlayPerimStatic(ca,mask4overlay);
%%

%% for each label, get intensity over the entire movie

stt = sum(cleanlb,[],3);
stt = stt>0;
slb = label(stt);
msrsum = measure(slb,cellfill,{'size','P2A'});
longids = find(msrsum.P2A>0.7);
sslb = single(slb);
for ii = 1:numel(longids)
    sslb(sslb==longids(ii)) = 0;
end
cleaner = dip_image(sslb);

mask2clean = cleaner>0;
cleanslb =ga.cleanUpMask_manual_square(squeeze(cellfill),cleaner>0,200);

tracelb= label(cleanslb);
trace = zeros(max(tracelb),size(fullseq,3));
traceraw = zeros(max(tracelb),size(fullseq,3));

wb = waitbar(0,'Quantifying Calcium Change in ROIs...');
for ll = 1:max(tracelb)%[1:8,10:max(lbl)]
%     tic
currmask = tracelb==ll;
bcurrmask = bdilation(currmask);
mask2use = repmat(bcurrmask,1,1,size(fullseq,3));
sumval = sum(fullseq,mask2use,[1 2]);
% imgINmask = image*mask2use;
% sumval = sum(imgINmask,[],[1 2]);
% sizeval = sum(bcurrmask,[],[1,2]);
sumtrace = single(sumval);
tracemed = median(sumtrace);
trace(ll,:) = sumtrace./tracemed;
traceraw(ll,:) = sumtrace;
% toc
waitbar(ll/max(tracelb),wb);
end
close(wb)

allsum = sum(traceraw,2);
[~, ordx] = sort(allsum, 'ascend');
ord_trace = traceraw(ordx,:);
% now plot
h = msgbox('Plotting HeatMap....');
times = 1:size(ord_trace,2);
figure;
% wb = waitbar(0,'Plotting some things...');
% cnt = 0;
for ii = 1:size(ord_trace,1)
    for jj = 1:size(ord_trace,2)
        p = patch([times(jj),times(jj)+1,times(jj)+1,times(jj)],[ii-1, ii-1, ii, ii],ord_trace(ii,jj));
        set(p,'FaceColor','flat','EdgeColor','none');
        
    end
end
title('GCamp Intensity (AU)','FontSize',16)
xlabel('Frame','FontSize',16);
ylabel('Puncta #','FontSize',16)
xlim([0 size(ord_trace,2)]);
ylim([0 size(ord_trace,1)]);
c = colorbar;
close(h);

values = cat(2,ordx,ord_trace);
mat2write = num2cell(values);
addrow = nan(1,size(mat2write,2));
addrow = num2cell(addrow);
addrow{1} = 'Label ID';
addrow{2} = 'Absolute Calcium Intensity per frame';
mat2write = [addrow;mat2write];
 savefile = 'C:\Users\sammy\Documents\MATLAB\calciumtest.xlsx';
 xlswrite(savefile,mat2write);