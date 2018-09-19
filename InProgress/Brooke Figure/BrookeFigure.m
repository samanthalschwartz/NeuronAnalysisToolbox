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
badis = maxsizes>10;
newlables = GeneralAnalysis.removeLabels(lb,msr.ID(badis));
sum_tt = sum(newlables,[],3)>0;
sumlb = label(sum_tt);
msr_sum = measure(sumlb,sumlb,{'DimensionsCube'});
maxsizes_sum = max(msr_sum.DimensionsCube);
badis = maxsizes_sum>15;
newSUMlables = GeneralAnalysis.removeLabels(sumlb,msr_sum.ID(badis));


%--check if it's OK
sumtt3 = repmat(newSUMlables>0,1,1,size(ca,3));
GeneralAnalysis.viewMaskOverlayPerimStatic(dip_image(ca),sumtt3);

%-- now quantify inside the mask
measure_structure = GeneralAnalysis.measureMaskTimeSeries(sumtt3(:,:,0:200),ca(:,:,1:201));
raw_intensities = measure_structure.Sum./measure_structure.Size;
intensities_t0 = mean(measure_structure.Sum(:,1:5)./measure_structure.Size(:,1:5),2);
intensities = raw_intensities./intensities_t0;
figure; hmap = heatmap(intensities);
hmap.GridVisible = 'off'
hmap.Colormap = jet(100);

% pks = cell(1,size(intensities,1));
peak_numbers = zeros(size(intensities,1),1);
for ii = 1:size(intensities,1)
    peaks = findpeaks(intensities(ii,:),'MinPeakHeight',1.2);
    peak_numbers(ii) = size(peaks,2);
end
%%
oldim = label(sumtt3(:,:,1));
tic
wow = changem(oldim,[-1;peak_numbers],[0,1:max(oldim)]);
toc
%% now analyze post light
ca_post = loadtiff('C:\Users\schwsama\Documents\Data\For Brooke\269dFP+281_light_post_2.tif');
% sv1 = findshift(a,b)
% sb = shift(b,sv2);
premask2 = squeeze(sumtt3(:,:,0));
sv1 = findshift(ca_post(:,:,1),ca(:,:,end))
newmask = shift(premask,sv1,1);
newmask = real(newmask);
newmask = newmask>0.1;
test = label(newmask);% ---- BAD!!! FT wraps around -- need to figure out what to do here!!
test(test == max(test)) = 0;
newmask = test;
% msrpost = measure(newmask,single(newmask),{'size'});
% msrpre = measure(premask,single(premask),{'size'});
% lbpost = label(newmask>0)
sum_post3 = repmat(newmask,1,1,size(ca_post,3));
%%
measure_structure_post = GeneralAnalysis.measureMaskTimeSeries(sum_post3(:,:,0:200),ca_post(:,:,1:201));
raw_intensities_post = measure_structure_post.Sum./measure_structure_post.Size;
% intensities_t0 = mean(measure_structure_post.Sum(:,1:5)./measure_structure_post.Size(:,1:5),2);
intensities_post = raw_intensities_post./intensities_t0;
figure; hmap = heatmap(intensities_post);
hmap.GridVisible = 'off'
hmap.Colormap = jet(100);

% pks = cell(1,size(intensities,1));
peak_numbers_post = zeros(size(intensities_post,1),1);
for ii = 1:size(intensities_post,1)
    peaks = findpeaks(intensities_post(ii,:),'MinPeakHeight',1.2);
    peak_numbers_post(ii) = size(peaks,2);
end
%%
tic
wow_post = changem(newmask,[-1;peak_numbers_post],[0,1:max(newmask)]);
toc
