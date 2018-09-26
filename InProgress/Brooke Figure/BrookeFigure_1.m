% --- first load in pre image and filter
ca = loadtiff('C:\Users\schwsama\Documents\Data\For Brooke\269dFP+281_light_pre_2.tif');
threshold_parameter = 1; % smaller number is more sensitive, larger number is less sensitive
ca_g = gaussf(ca,[1 1 0]);
ca_gnorm = ca_g;%- medif(ca_g);
ca_dzz = dzz(ca_gnorm,[2 2 3]);
ab_ca_dzz = abs(ca_dzz);
gab_ca_dzz = gaussf(ab_ca_dzz,[1 1 0]);
tt = threshold(ab_ca_dzz^threshold_parameter,'otsu');
lb = label(tt,inf,0,100);
msr = measure(lb,ca,{'DimensionsCube'});
maxsizes = max(msr.DimensionsCube);
badis = maxsizes>20;
newlables = GeneralAnalysis.removeLabels(lb,msr.ID(badis));
sum_tt = sum(newlables,[],3)>0;
sumlb = label(sum_tt);
msr_sum = measure(sumlb,sumlb,{'DimensionsCube','Size'});
maxsizes_sum = max(msr_sum.DimensionsCube);
badis = maxsizes_sum>15 | msr_sum.Size<3;
newSUMlables = GeneralAnalysis.removeLabels(sumlb,msr_sum.ID(badis));
% newSUMlables = bdilation(newSUMlables>0,1);
newSUMlables = label(newSUMlables>0,1);

%-- now load in post image. find shift. shift mask and do another round of
%mask size thresholding
ca_post = loadtiff('C:\Users\schwsama\Documents\Data\For Brooke\269dFP+281_light_post_2.tif');
sv1 = findshift(ca_post(:,:,1),ca(:,:,end));
padwidth = max(sv1)*10;
postSUMlabels = shift_withpadding(squeeze(newSUMlables),sv1,padwidth);
postSUMlabels = round(real(postSUMlabels));
%-first ask: are all the same labels in both pre and post
p1 = single(unique(postSUMlabels));
p0 = single(unique(newSUMlables));
msr_post = measure(postSUMlabels,postSUMlabels,{'Size'});
badids1 = msr_post.ID(msr_post.Size<3);
badids2 = setdiff(p1,p0);
badidspost = unique([badids1,badids2]);
if ~isempty(badids1|badids2)
newSUMlables = GeneralAnalysis.removeLabels(newSUMlables,badidspost);
postSUMlabels = GeneralAnalysis.removeLabels(postSUMlabels,badidspost);
end
finalLabels_pre = repmat(newSUMlables,1,1,size(ca,3));
finalLabels_post = repmat(postSUMlabels,1,1,size(ca,3));

%--check if it's OK
GeneralAnalysis.viewMaskOverlayPerimStatic(dip_image(ca),finalLabels_pre>0);


%-- now quantify inside the mask
measure_structure = GeneralAnalysis.measureMaskTimeSeries(finalLabels_pre(:,:,0:200),ca(:,:,1:201));
raw_intensities = measure_structure.Sum./measure_structure.Size;
intensities_t0 = mean(measure_structure.Sum(:,1:5)./measure_structure.Size(:,1:5),2);
intensities = raw_intensities./intensities_t0;
figure; hmap = heatmap(intensities);
hmap.GridVisible = 'off';
hmap.Colormap = jet(100);

%%
measure_structure_post = GeneralAnalysis.measureMaskTimeSeries(finalLabels_post(:,:,0:200),ca_post(:,:,1:201));
raw_intensities_post = measure_structure_post.Sum./measure_structure_post.Size;
% intensities_t0 = mean(measure_structure_post.Sum(:,1:5)./measure_structure_post.Size(:,1:5),2);
intensities_post = raw_intensities_post./intensities_t0;
figure; hmap = heatmap(intensities_post);
hmap.GridVisible = 'off'
hmap.Colormap = jet(100);
%%
% pks = cell(1,size(intensities,1));
peak_numbers = zeros(size(intensities,1),1);
for ii = 1:size(intensities,1)
    allpeaks = findpeaks(intensities(ii,:));
    peaks = findpeaks(intensities(ii,:),'MinPeakHeight',1.1,'MinPeakProminence',2*std(allpeaks));
    peak_numbers(ii) = size(peaks,2);
end
%%
% pks = cell(1,size(intensities,1));
peak_numbers_post = zeros(size(intensities_post,1),1);
for ii = 1:size(intensities_post,1)
    allpeaks_post = findpeaks(intensities_post(ii,:));
    peaks = findpeaks(intensities_post(ii,:),'MinPeakHeight',1.1,'MinPeakProminence',2*std(allpeaks_post));
    peak_numbers_post(ii) = size(peaks,2);
end
%%
tic
wow_pre = changem(finalLabels_pre(:,:,1),[-1;peak_numbers],[single(unique(finalLabels_pre))']);
toc
tic
wow_post = changem(finalLabels_post(:,:,1),[-1;peak_numbers_post],[single(unique(finalLabels_post))']);
toc
%%
dipshow(wow_pre,[-1 10]); dipmapping(gcf,'jetblack')
dipshow(wow_post,[-1 10]); dipmapping(gcf,'jetblack')

%%
peaknum_forplot = peak_numbers - peak_numbers_post;
tic
wow_forplot = changem(finalLabels_post(:,:,1),[-1;peaknum_forplot],[single(unique(finalLabels_post))']);
toc
%% test values
lbids = [single(unique(finalLabels_post(:,:,1)))']; 
clmp = flip(red(max(peaknum_forplot)*2));
clmpneg = flip(blue(20));

clmP = hot(max(peaknum_forplot)-min(peaknum_forplot));
dipshow(ca)
for ll = 1:size(lbids(2:end),1)
currlb = lbids(ll+1);
col = peaknum_forplot(ll);
% if col<0
% %     facecolor = clmpneg(abs(col),:);
%     facecolor = [0 0 1];
% elseif col==0
%     facecolor = [0 0 .1];
% else
%     facecolor = clmp(col*2,:);
% end

facecolor = clmP(

vts = bwboundaries(imbinarize(single(finalLabels_post(:,:,1) == ll)));
x = vts{1}(:,2);
y = vts{1}(:,1);
patch(x,y,facecolor,'FaceAlpha',.8)
end
