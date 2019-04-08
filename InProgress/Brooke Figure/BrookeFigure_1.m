% --- first load in pre image and filter
savedir = 'C:\Users\schwsama\Documents\Data\For Brooke\';
ca = loadtiff('C:\Users\schwsama\Documents\Data\For Brooke\269dFP+281_light_pre_2.tif');
threshold_parameter = 1; % smaller number is more sensitive, larger number is less sensitive
ca_g = gaussf(ca,[2 2 0]);
ca_gnorm = ca_g;%- medif(ca_g);
ca_dzz = dzz(ca_gnorm,[1 1 3]);
ab_ca_dzz = abs(ca_dzz);
gab_ca_dzz = gaussf(ab_ca_dzz,[1 1 0]);
tt = threshold(ab_ca_dzz^threshold_parameter,'otsu');
lb = label(tt,inf,0,100);
msr = measure(lb,ca,{'DimensionsCube'});
maxsizes = max(msr.DimensionsCube);
badis = maxsizes>30;
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
GeneralAnalysis.viewMaskOverlayPerimStatic(dip_image(ca_post),finalLabels_post>0);


%-- now quantify inside the mask
measure_structure = GeneralAnalysis.measureMaskTimeSeries(finalLabels_pre(:,:,:),ca(:,:,:));
raw_intensities = measure_structure.Sum./measure_structure.Size;
intensities_t0 = mean(measure_structure.Sum(:,1:5)./measure_structure.Size(:,1:5),2);
intensities = raw_intensities./intensities_t0;
figure; hmap = heatmap(intensities);
hmap.GridVisible = 'off';
hmap.Colormap = jet(100);

%%
measure_structure_post = GeneralAnalysis.measureMaskTimeSeries(finalLabels_post(:,:,:),ca_post(:,:,:));
raw_intensities_post = measure_structure_post.Sum./measure_structure_post.Size;
% intensities_t0 = mean(measure_structure_post.Sum(:,1:5)./measure_structure_post.Size(:,1:5),2);
intensities_post = raw_intensities_post./intensities_t0;
figure; hmap = heatmap(intensities_post);
hmap.GridVisible = 'off';
hmap.Colormap = jet(100);
%%
% pks = cell(1,size(intensities,1));
peak_numbers = zeros(size(intensities,1),1);
for ii = 1:size(intensities,1)
    
    trace = movmean(intensities(ii,:),3);
    allpeaks = findpeaks(trace);
    srtpks = sort(allpeaks);
    mean_pre = mean(srtpks(1:50));
    std_pre = std(srtpks(1:50));
    [pk,lk] = findpeaks(trace,'MinPeakHeight',1.2*mean_pre,'MinPeakProminence',std_pre*3);
    peak_numbers(ii) = size(pk,2);
    f = figure; plot(1:size(intensities(ii,:),2),trace,lk,pk,'o');
    title(['Intensity Pre Light - Trace # ' num2str(ii)]);
    xlabel('Frame');
    ylabel('F/F0');
    saveas(f,fullfile(savedir,['PreTrace_' num2str(ii)]),'png');
    close(f);
end
%%
% pks = cell(1,size(intensities,1));
peak_numbers_post = zeros(size(intensities_post,1),1);
for ii = 1:size(intensities_post,1)
    trace = movmean(intensities_post(ii,:),3);
%     [allpeaks_post,locs] = findpeaks(trace,'MinPeakProminence',0.2);
    [allpeaks_post,locs] = findpeaks(trace);
    srpks = sort(allpeaks_post);
    mean_post = mean(srpks(1:100));
    std_post = std(srpks(1:100));
    [pk,lk]  = findpeaks(trace,'MinPeakHeight',1.2*mean_post,'MinPeakProminence',3*std_post);
    peak_numbers_post(ii) = size(pk,2);
    f = figure; plot(1:size(intensities_post(ii,:),2),trace,lk,pk,'o');
    title(['Intensity Post Light - Trace # ' num2str(ii)]);
    xlabel('Frame');
    ylabel('F/F0');
    saveas(f,fullfile(savedir,['PostTrace_' num2str(ii)]),'png');
    close(f);
end
%%
% zeroIDs=  peak_numbers==0;
% peak_numbers(zeroIDs) = [];
% peak_numbers_post(zeroIDs)= [];
% 
% prelabs = single(unique(finalLabels_pre))';
% prelabs(zeroIDs) = [];
% postlabs = single(unique(finalLabels_post))';
% postlabs(zeroIDs) = [];


% tic
% wow_pre = changem(finalLabels_pre(:,:,1),[-1;peak_numbers],[prelabs]);
% toc
% tic
% wow_post = changem(finalLabels_post(:,:,1),[-1;peak_numbers_post],[postlabs]);
% toc
%%
% dipshow(wow_pre,[-1 10]); 
% dipshow(wow_post,[-1 10]); 

%%
peaknum_forplot_pre = peak_numbers - peak_numbers_post;
peaknum_forplot = peaknum_forplot_pre;
% peaknum_forplot = abs(min(peaknum_forplot_pre)) + peaknum_forplot_pre
% tic
% wow_forplot = changem(finalLabels_post(:,:,1),[-1;peaknum_forplot],[single(unique(finalLabels_post))']);
% toc
%% test values
lbids = [single(unique(finalLabels_post(:,:,1)))']; 
clmp = hot(ceil(max(peaknum_forplot).*1)); 
clmp = clmp(end-max(peaknum_forplot):end,:)
clmpneg = winter(abs(min(peaknum_forplot)));

clmP = hot(max(peaknum_forplot)-min(peaknum_forplot));
dipshow(ca)
for ll = 1:size(lbids(2:end),1)
currlb = lbids(ll+1);
col = peaknum_forplot(ll);
if col<0
    facecolor = clmpneg(abs(col),:);
%     facecolor = [0 0 1];
elseif col==0
    facecolor = [0 .2 1];
else
    facecolor = clmp(ceil(col),:);
end

% facecolor = clmP(
% 
vts = bwboundaries(imbinarize(single(bdilation((finalLabels_pre(:,:,1) == ll),2))));
x = vts{1}(:,2);
y = vts{1}(:,1);
patch(x,y,facecolor,'FaceAlpha',.9)
end
dipmapping(gcf,'colormap',flip(gray))

%
figure;
overallcolor  = [clmpneg;[0 .2 0];clmp];
for yy = 1:size(overallcolor,1)
    facecolor = overallcolor(yy,:)
    patch([1 3 3 1],[yy yy yy+1 yy+1],facecolor,'EdgeColor','none')
    
end
sortvals = sort(unique(peaknum_forplot));
yticks(1:size(peaknum_forplot,1))

yticklabels(sortvals(1):sortvals(end));
