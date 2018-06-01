close all
clear all
%% addfilepaths
addpath(genpath('G:\Sam\Kennedy MATLAB folder\Sam'));
filepath488 = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\maxproj_488light\homerGFP_+488_maxproj_cell6.tif';
filepath561 = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\maxproj_488light\GluA1_+488_maxproj_cell6.tif';
filepath561pre = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\maxproj_pre488\GluA1_dark_maxproj_cell6.tif';
filepath561post = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\maxproj_post488\GluA1_post488_maxproj_cell6.tif';
savedir = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\results';
savename = 'cell6';
%% load files

img488 = GeneralAnalysis.loadtiff_1ch(filepath488);
img561 = GeneralAnalysis.loadtiff_1ch(filepath561);
img561pre = GeneralAnalysis.loadtiff_1ch(filepath561pre);
img561post = GeneralAnalysis.loadtiff_1ch(filepath561post);

%% correct for drift over time from 561 channel and apply the same correction to 488 channel 
[shifted_img561,sv_arr] = GeneralAnalysis.timedriftCorrect(img561);
shifted_img488 = GeneralAnalysis.applydriftCorrect(img488(:,:,1:end),sv_arr);
shifted_img488 = cat(3,img488(:,:,0),shifted_img488);
% shift pre488 files onto 488 channel
svpre = findshift(shifted_img561(:,:,0),img561pre(:,:,0));
img561pre(:,:,0) = shift(img561pre(:,:,0),svpre);
[shifted_img561pre,sv_arrpre] = GeneralAnalysis.timedriftCorrect(img561pre);
% shift post488 files onto 488 channel
svpost = findshift(shifted_img561(:,:,end),img561post(:,:,0));
img561post(:,:,0) = shift(img561post(:,:,0),svpost);
[shifted_img561post,sv_arrpost] = GeneralAnalysis.timedriftCorrect(img561post);
% create concatenated movie
fullimg561 = cat(3,shifted_img561pre,shifted_img561,shifted_img561post);
h = dipshow(fullimg561,[0 1000]);
%% make some masks using Homer1c channel (488) -- need to check and make sure that this masking is working for each file
gsig = [1 1 1];
lsig = [1 1 1];
img488_1 = medif(shifted_img488,3);
img488_2 = gaussf(img488_1,gsig);
img4lcutoff = img488_2;
img488_3 = GeneralAnalysis.imgLaplaceCutoff(img4lcutoff,gsig,lsig);
img488_4 = img488_3;%.^1.2;
img4mask = img488_4;
maskpre1 = mask_img(img4mask);
%  make labels bigger
maskpre = GeneralAnalysis.bwmorph_timeseries(maskpre1,'thicken',1);

ll = slice_op('watershed',-img4lcutoff,2);
maskws = maskpre;
maskws(ll) = 0;
cm4overlay = hot(256);
cm4overlay(end,:) = [1 1 0];
overl = overlay(img488_4,maskws);
h = dipshow(overl,cm4overlay);
dipmapping(h,[100 600])
mask = logical(maskws);

%%  now loop through and label the image
lbl = label(mask,2);
goodid = [];
goodlbl = lbl;
allgoodmask = mask;
for ii = 1:max(lbl)
    % check that label is in every frame
    out = find(lbl==ii);
    [x, y, z] = ind2sub(size(lbl),out);
    Lia = ismember(1:size(lbl,3),z);
    if sum(Lia)~=size(lbl,3)
        goodlbl(lbl==ii) = false;
        allgoodmask(lbl==ii) = 0;
        continue;
    end
    goodid = [goodid;ii];
end
%%

pre_lbl = repmat(goodlbl(:,:,0),1,1,size(shifted_img561pre,3));
post_lbl = repmat(goodlbl(:,:,end),1,1,size(shifted_img561post,3));
full_lbl = cat(3,pre_lbl,goodlbl,post_lbl);
rawIntensity_traces = zeros(size(full_lbl,3),numel(goodid));
meanIntensity_traces = zeros(size(full_lbl,3),numel(goodid));
full_mask = full_lbl;
full_mask(full_mask>0) = 1;
for jj=1:numel(goodid)
    gi = goodid(jj);
    vals = squeeze(sum(fullimg561,full_lbl==gi,[1 2]));
    sums = squeeze(sum(full_mask,full_lbl==gi,[1 2]));
    means = vals./sums;
    rawIntensity_traces(:,jj) = single(vals);
    meanIntensity_traces(:,jj) = single(means);
end

normraw = rawIntensity_traces./mean(rawIntensity_traces(1:size(pre_lbl,3),:));
normmean = meanIntensity_traces./mean(meanIntensity_traces(1:size(pre_lbl,3),:));
figure; plot(normraw); hold on; plot(mean(normraw,2),'--','LineWidth',2,'Color','k');
figure; plot(normmean);hold on; plot(mean(normmean,2),'--','LineWidth',2,'Color','k');

%% -- can start from here after loading saved file
minval = 1.5;
nr_ids = sum(normraw>minval,1);
nm_ids = sum(normmean>minval,1);
figure; plot(normraw(:,nr_ids>0)); hold on; plot(mean(normraw(:,nr_ids>0),2),'--','LineWidth',2,'Color','k');
figure; plot(normmean(:,nm_ids>0));hold on; plot(mean(normmean(:,nm_ids>0),2),'--','LineWidth',2,'Color','k');
figure;  plot(mean(normmean(:,nm_ids>0),2),'--','LineWidth',2,'Color','k');
[f,x] = ecdf(max(normmean,[],1));
figure; plot(x,f);
%%
%% now find 'good' rois and plot ontop of 488 image the perimeter of the mask
gdtest = sum(normmean>minval,1);
gdids = goodid(gdtest>0);
% 
% pickednorms = norms(gdtest>0,:);
% figure; plot(pickednorms');
% figure; plot(mean(pickednorms)');

%%
% plot max value distribution
% [f,x] = ecdf(max(norms,[],2));
% figure;plot(x,f)
% xlabel('Fold Increase in GluA1 Fluroscence at Homer1 Puncta');
% ylabel(['Cumulative Probability']);

% make perimeter for responders
perim_mask = dip_image(zeros(size(full_lbl)));
for ll = 1:size(gdids,1)
    perim = maskperim(full_lbl==gdids(ll));
    perim_mask = perim_mask | perim;
end
% make perimeter for non-responders
nonresp_perim = maskperim(logical(full_mask));


ovl488 = (medif(fullimg561,2));
% ovl561(nonresp_perim) = max(ovl561)*1.5;
ovl488(nonresp_perim) = 0;
ovl488(perim_mask) = max(ovl488)*20;
cm4overlay = bone(256);
cm4overlay(1,:) = [0 0 1]; %set all perims
cm4overlay(end,:) = [1 0 0]; %set all responder perims
h1 = dipshow(ovl488,cm4overlay);
dipmapping(h1,'global');
dipmapping(h1,'lin');
dipmapping(h1,[100 1000]);
diptruesize(h1,150)
% % perim_mask
ovlglua1 = overlay(medif(shifted_img488,3),nonresp_perim(:,:,(size(shifted_img561pre,3)):(end - size(shifted_img561post,3))));
h2 = dipshow(ovlglua1,cm4overlay);
dipmapping(h2,'global');
dipmapping(h2,'lin');
dipmapping(h2,[100 600]);
diptruesize(h2,150)

% %% plot cdf of % responding at different threshold levels
% [f,x] = ecdf(max(norms,[],1));
% figure; plot(x,f);
%% save some info
clear results
results.fullimg561 = single(fullimg561);
results.shifted_img488 = single(shifted_img488);
results.fullmask = single(full_mask);
results.full_lbl = single(full_lbl);
results.rawIntensity_traces = rawIntensity_traces;
results.meanIntensity_traces = meanIntensity_traces;
results.normraw = normraw;
results.normmean = normmean;
results.presize = size(shifted_img561pre,3);
results.postsize = size(shifted_img561post,3);
save(fullfile(savedir,savename),'results');
%% make movies
addpath('G:\Sam\OlderMatlabFiles\sam\2016');
%-- GluA1 results with cuttoff
GluA1filename = fullfile(savedir,[savename '_GluA1cutoff=' num2str(minval)]);
writeDipImageMovie(h1,GluA1filename)
%-- homer with masks movie
homer1filename = fullfile(savedir,[savename '_Homer1c']);
writeDipImageMovie(h2,homer1filename)

