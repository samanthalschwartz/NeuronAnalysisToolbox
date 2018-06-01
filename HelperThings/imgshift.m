close all
clear all
addpath(genpath('G:\Sam\Kennedy MATLAB folder\Sam'));
filepath488 = 'Z:\Sam\Data!\180516_GluA1mCh-HomerGFP\maxproj_488light\GluA1_+488_maxproj_cell3.tif';
filepath561 = 'Z:\Sam\Data!\180516_GluA1mCh-HomerGFP\maxproj_488light\homerGFP_+488_maxproj_cell3.tif';
img488 = GeneralAnalysis.loadtiff_1ch(filepath488);
img561 = GeneralAnalysis.loadtiff_1ch(filepath561);

shifted_im488 = dip_image(zeros(size(img488)));
shifted_im561_b = dip_image(zeros(size(img561)));
imgref488 = squeeze(img488(:,:,0));
imgref561 = squeeze(img561(:,:,0));
shifted_im488(:,:,0) = imgref488;
shifted_im561_b(:,:,0) = imgref561;

for ii = 1:(size(img488,3)-1)
    imgcurr488= squeeze(img488(:,:,ii));
    imgcurr561= squeeze(img561(:,:,ii));
    sv1 = findshift(imgref488,imgcurr488,'iter',0);
    sv2 = findshift(imgref561,imgcurr561,'iter',0);
    shift_img488 = shift(imgcurr488,sv1,1);
    shift_img561 = shift(imgcurr561,sv2,1);
    shifted_im488(:,:,ii) = shift_img488;
    shifted_im561_b(:,:,ii) = shift_img561;
end


shifted_im488(shifted_im488 > (2^14 - 2*255)) = 0; %remove weird saturated pixels.
shifted_im561_b(shifted_im561_b > (2^14 - 2*255)) = 0; %remove weird saturated pixels.
shifted_im488(shifted_im488<0) = 0;
shifted_im561_b(shifted_im561_b<0) = 0;
out = joinchannels('rgb',shifted_im488,shifted_im561_b);
% dipshow(out,[-100 2250]);
%%
medfilt = medif(shifted_im561_b);
gsig = [2 2 1];
lsig = [2 2 0];
img_out = LoG_threshold(shifted_im561_b,gsig,lsig);
mask = mask_img(img_out);
cm4overlay = hot(256);
cm4overlay(end,:) = [1 1 0];
overl = overlay(shifted_im561_b,mask);
h = dipshow(overl,cm4overlay);
dipmapping(h,[100 600])
%%  now loop through and label the image
lbl = label(mask,2);
goodid = [];
for ii = 1:max(lbl)
    % check that label is in every frame
    out = find(lbl==ii);
    [x, y, z] = ind2sub(size(lbl),out);
    Lia = ismember(1:size(lbl,3),z);
    if sum(Lia)~=size(lbl,3)
        continue;
    else
    goodid = [goodid;ii];
    end
end

recruit_traces = zeros(numel(goodid),size(shifted_im488,3));
for jj = 1:numel(goodid)
    gi = goodid(jj);
%     Fo = squeeze(sum(shifted_im488(:,:,0:2),lbl(:,:,0:2)==gi,[1 2])/sum(mask(:,:,0:2),lbl(:,:,0:2)==gi,[1 2]));
    F = squeeze(sum(shifted_im488,lbl==gi,[1 2])/sum(mask,lbl==gi,[1 2]));
    recruit_traces(jj,:) = single(F);
%     Fo = squeeze(sum(shifted_im488(:,:,0),lbl(:,:,0)==gi,[1 2]));
%     F = squeeze(sum(shifted_im488,lbl==gi,[1 2]));
%     recruit_traces(jj,:) = single(F./Fo);
end
figure; plot(recruit_traces');
norms = recruit_traces./recruit_traces(:,1);
figure; plot(mean(norms)');
figure; plot(norms');
%% now find 'good' rois and plot ontop of 488 image the perimeter of the mask
gdtest = sum(norms>1.5,2);
gdids = goodid(gdtest>0);

pickednorms = norms(gdtest>0,:);
figure; plot(pickednorms');
figure; plot(mean(pickednorms)');
% make perimeter for responders
perim_mask = dip_image(zeros(size(lbl)));

%%
% plot max value distribution
[f,x] = ecdf(max(norms,[],1));
figure;plot(x,f)
for ll = 1:size(gdids,1)
    perim = maskperim(lbl==gdids(ll));
    perim_mask = perim_mask | perim;
end
% make perimeter for non-responders
nonresp_perim = maskperim(mask);


ovl488 = shifted_im488;
% ovl561(nonresp_perim) = max(ovl561)*1.5;
ovl488(nonresp_perim) = 0;
ovl488(perim_mask) = max(ovl488)*1.5;;
cm4overlay = bone(256);
cm4overlay(1,:) = [1 1 0]; %set all perims
cm4overlay(end,:) = [0 0 1]; %set all responder perims

h = dipshow(ovl488,cm4overlay);
dipmapping(h,'global')
dipmapping(h,'lin')
dipmapping(h,[100 1000]);
% perim_mask
overlay(shifted_im561_b,nonresp_perim)
%% plot cdf of % responding at different threshold levels
[f,x] = ecdf(max(norms,[],1));
figure; plot(x,f);