%% load and extract images
% uiopen('F:\Hiester et al., 2017 Frontiers Submission\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
image = loadtiff('F:\Hiester et al., 2017 Frontiers Submission\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_3.tif');% uiopen('E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
sep = image(:,:,:,3);
cellfill = image(:,:,:,1);
ga = GeneralAnalysis;
sepim = dip_image(sep);
%% mask cell fille
msk_cellfill = threshold(gaussf(cellfill),'otsu');
max_msk_cellfill = max(msk_cellfill,[],3);
%% mask spines in time movie
% dggsep = ga.imgDggCutoff(gaussf(sep));
% varsep = varif(dggsep,7,'elliptic');
% spines = threshold(varsep,'otsu');
% spines_manual = ga.imgThreshold_fixedUserInput(varsep)
%% try white top hat
image_out = opening(sep,5,'elliptic');
wth = sep - image_out;
dggsep = ga.imgDggCutoff(wth);
varsep = varif(dggsep,5,'elliptic');
spines = threshold(varsep.^.8,'otsu');
h = ga.viewMaskOverlayPerimStatic(sepim,spines)
dipmapping(h,'colormap',jet)
dipmapping(h,[0 2500]);
diptruesize(h,100);

%% mask spines with max projection to make time course
max_sep = max(sep,[],3);
max_spines = max(spines,[],3);
wshed = watershed(-gaussf(max_sep),2);
max_spines(wshed) = 0;
lb = label(max_spines,1,30); %label with minimum size 30 pixels
spines_fixed = repmat(lb>0,1,1,size(sep,3));
%% view mask
h = ga.viewMaskOverlayPerimStatic(dip_image(sep),spines_fixed)
dipmapping(h,'colormap',jet)
dipmapping(h,[0 2500]);
diptruesize(h,100);
%% sum each ROI over time
spinelabel = label(spines_fixed);
sums = zeros(max(spinelabel),size(sep,3));
wb = waitbar(0);
sepim = dip_image(sep);
for ll = 1:max(spinelabel)
sums(ll,:) = single(sum(sepim,spinelabel == ll,[1 2]));
waitbar(ll/max(spinelabel),wb);
end
close(wb)


msrarray = cell(1,size(sep,3));
wb = waitbar(0);
for ll = 0:(size(sep,3)-1)
    msr = measure(spinelabel(:,:,ll),sepim(:,:,ll),{'size','sum'});
    msrarray{ll+1} = msr;
    waitbar((ll+1)/size(sep,3),wb);
end
close(wb)

sums = zeros(max(spinelabel),size(sep,3));
sizes = [];
%reshape array
for ll = 1:max(spinelabel)
    for tt = 1:size(sep,3)
        
        sums(ll,tt) = msrarray{tt}(ll).sum;
        
    end
end



f = figure; 
trace = sums./sums(:,1);
allsum = sum(trace,2);
[~, ordx] = sort(allsum, 'descend');
ord_trace = trace(ordx,:);

hm = heatmap(ord_trace)
hm.GridVisible = 'off';
hm.Colormap = jet(50);
hm.ColorLimits = [0 4];

%% get boundaries
% logspines = logical(max_spines);
% b = bwboundaries(logspines);
% figure; hold on;
% for ii =1:numel(b)-2
%     boundary = b{ii};
%     plot(boundary(:,2),...
%         boundary(:,1),'g','LineWidth',2);
%     hold on;
% end
% xlim([1 size(sep,2)]);
% ylim([1 size(sep,1)]);
%% view rois from file
strFilename = 'E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_3_ROI1.zip';
[sROI] = ReadImageJROI(strFilename);
% fh = dipshow(sep); hold on;
for ii = 1:numel(sROI)
ysize = sROI{ii}.vnRectBounds(3)-sROI{ii}.vnRectBounds(1);
xsize = sROI{ii}.vnRectBounds(4)-sROI{ii}.vnRectBounds(2);
rectbounds = [sROI{ii}.vnRectBounds(2),sROI{ii}.vnRectBounds(1),xsize,ysize];
rectangle('Position',rectbounds,'EdgeColor','g')
end

%%
spinesum = sum(spines,[],3);
spinesummask = bdilation(spinesum>0,2);
spinesummask3 = repmat(spinesummask,1,1,size(sep,3));
sumshaft = sum(sepim,max_msk_cellfill-spinesummask3,[1 2]);
sumshaft./sumshaft(1)

sumspine = sum(sepim,spinelabel>0,[1 2]);
sumspine./sumspine(1)