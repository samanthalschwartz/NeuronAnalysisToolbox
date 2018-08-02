ga = GeneralAnalysis();
[FILENAME, PATHNAME] = uigetfile(fullfile(pwd,'*.*'),'Select an example calcium image');
prompt = 'Select the Unique File Identifier: use * for wildcard';
name = 'Calcium Data Selection';
defaultanswer = {FILENAME};
numlines = 1;
flstring=inputdlg(prompt,name,numlines,defaultanswer);
flstring = flstring{1};

[FILENAME_cf, PATHNAME_cf] = uigetfile(fullfile(PATHNAME,'*.*'),'Select the corresponding cell fill image');

img = dip_image(ga.loadtiffseries(PATHNAME,flstring));
cellfill_pre = ga.loadtiff_1ch(fullfile(PATHNAME_cf,FILENAME_cf));

% align cellfill to img
imgsum = sum(img,[],3);
corrimg = cat(3,imgsum,cellfill_pre);
dc_corrimg = ga.timedriftCorrect(corrimg);
cellfill = dc_corrimg(:,:,end);


imgg = gaussf(img,[1 1 0]);
diffimage = dzz(imgg,[2 2 0]);
adiff = abs(diffimage);
lp = ga.imgLaplaceCutoff(adiff,[1 1 0],[1 1 0]);
tt = threshold(lp,'otsu');
lb = label(tt,inf,30,600);
% --- optional thresholding
tic
msr = measure(lb,adiff,'P2A','size','sum','feret');
toc

%% bill analysis
[f,x] = ecdf(msr.size);


%%


ids = find(msr.P2A>1.2 | msr.P2A<0.8);
cleanlb = lb;
tic
for ii = ids
cleanlb(lb==ii) = 0;
end
toc
% -------


tot = sum(cleanlb,[],3);
tot = tot>0;
% ga.viewMaskOverlayPerimStatic(cellfill,bdilation(tot,1));

totseries = repmat(bdilation(tot,1),1,1,size(img,3));
% ga.viewMaskOverlayPerimStatic(cellfill,bdilation(tot,1))
% ga.viewMaskOverlayPerimStatic(img,totseries)
llb = label(tot);

mskin = bdilation(tot,1);
newmask = ga.cleanUpMask_manual(cellfill,bdilation(tot,1));
llb = label(newmask);

%% for each label, get intensity over the entire movie
trace = zeros(max(llb),size(img,3));
wb = waitbar(0,'Quantifying Calcium Change in ROIs...');
for ll = 1:max(llb)%[1:8,10:max(lbl)]
%     tic
currmask = llb==ll;
bcurrmask = bdilation(currmask);
mask2use = repmat(bcurrmask,1,1,size(img,3));
sumval = sum(img,mask2use,[1 2]);
% imgINmask = image*mask2use;
% sumval = sum(imgINmask,[],[1 2]);
% sizeval = sum(bcurrmask,[],[1,2]);
tracemed = median(sumval);
trace(ll,:) = sumval./tracemed;
% toc
waitbar(ll/max(llb),wb);
end
close(wb)

allsum = sum(trace,2);
[~, ordx] = sort(allsum, 'ascend');
ord_trace = trace(ordx,:);
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
