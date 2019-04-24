clear all
filename = 'G:\FromMicroscopeComputer\190415 Cry2Olig-GephIB-HaloTag-CAGregulated\clone2-bluelight_timeseries\clone2-bluelight_timeseries_w2640_s1_t.tiff';
uiopen(filename,1);
timage = GeneralAnalysis.timedriftCorrect(image);
GeneralAnalysis.LibTiff(timage,[filename(1:end-5) '_Shift' filename(end-4:end)]);

% maxval = prctile(image(:),99.99)*2;
newim = timage;
newim(newim<0) = 0;
firstfr = newim(:,:,0);
% newim(newim>maxval) = 0;
dipshow(newim);
img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(firstfr,1,1);
normim = img_laplcutoff./max(img_laplcutoff);
mask = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff);
% mask = bdilation(mask,1);
GeneralAnalysis.viewMaskOverlay(firstfr,mask);

fullmask = repmat(mask,[1 1 size(newim,3)]);
GeneralAnalysis.viewMaskOverlay(newim,fullmask);

lb = label(fullmask,1);
vals = zeros(max(lb),size(newim,3));
wb = waitbar(0);
for ll = 1:max(lb)
currob = newim.*(lb==ll);

vals(ll,:) = single(squeeze(sum(newim,lb==ll,[1 2])));
waitbar(ll/max(lb),wb);
end
close(wb);
forplot = vals./vals(:,1);
test = sort(forplot,1,'descend');

figure; plot(test')

figure; h = heatmap(test(1:end-20,:));
set(h,'GridVisible','off')
set(h,'Colormap',jet)
set(h,'ColorLimits',[1 1.8])
%%
val = 20;
test1 = mean(test(1:val,:));
testsd = std(test(1:val,:))./sqrt(val); 
xval = [0:0.5:(size(test1,2)/2 )]; xval = xval(1:end-1);
assert(size(xval,2) == size(test1,2))
figure;e = errorbar(xval,[test1],[testsd]);
e.CapSize = 2;
e.LineWidth = 0.3;
e.Color = 'k';
hold on; 
p = plot(xval,test1);
p.LineWidth = 1.5;
p.Color = 'k';
xlim([0 xval(end)]); ylim([1 1.7]);
set(gca,'FontSize',20)
ylabel('

