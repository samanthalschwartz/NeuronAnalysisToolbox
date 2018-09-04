%% load and extract images
uiopen('F:\Hiester et al., 2017 Frontiers Submission\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
% uiopen('E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1);
sep = image(:,:,:,3);
cellfill = image(:,:,:,1);
ga = GeneralAnalysis;
%% mask cell fille
msk_cellfill = threshold(gaussf(cellfill),'otsu');
%% mask spines in time movie
dggsep = ga.imgDggCutoff(gaussf(sep));
varsep = varif(dggsep,7,'elliptic');
spines = threshold(varsep,'otsu');
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

f = figure; 
hm = heatmap(sums./sums(:,1))
hm.GridVisible = 'off';
hm.Colormap = jet(50);
% hm.ColorLimits = [0 3]

%% get boundaries
logspines = logical(max_spines);
b = bwboundaries(logspines);
figure; hold on;
for ii =1:numel(b)-2
    boundary = b{ii};
    plot(boundary(:,2),...
        boundary(:,1),'g','LineWidth',2);
    hold on;
end
xlim([1 size(sep,2)]);
ylim([1 size(sep,1)]);


