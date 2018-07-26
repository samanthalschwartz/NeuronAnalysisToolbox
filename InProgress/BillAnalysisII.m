%% load data and prepare for analysis
datafile = 'E:\Sam\Data\Bill\Composite_pos875.tif';
ga = GeneralAnalysis();
full = ga.loadtiff_2ch(datafile);
cellfill= image(:,:,:,1);
geph = image(:,:,:,3);
[cellfill,sv] = ga.timedriftCorrect(dip_image(cellfill));
geph = ga.applydriftCorrect(dip_image(geph),sv);
dipshow(cellfill)
dipshow(geph)
gsig = 1;
lsig = 1;
gephl = ga.imgLaplaceCutoff(geph,gsig,lsig);
cellfilll = ga.imgLaplaceCutoff(cellfill,gsig,lsig);
[geph_mask1,threshval] = ga.imgThreshold_fixedUserInput(gephl);
[cellfill_mask,threshval] = ga.imgThreshold_fixedUserInput(cellfilll);
cellfillmask = logical(ga.thicken(cellfill_mask,2));
geph_mask2 = geph_mask1*cellfillmask;
gephmask = bdilation(geph_mask2);
[h,overlayim] = ga.overlay(geph,gephmask);
[h,overlayim] = ga.viewMaskOverlayPerimStatic(geph,gephmask)
gephmask_cleaned = ga.cleanUpMask_manual(geph,gephmask);
gephwtshed = ga.watershed_timeseries(-gaussf(geph),1);
gephclean2 = gephmask_cleaned;
gephclean2(gephwtshed==1) = 0;
%%
c = stretch(cellfill);
g= stretch(geph);
ovl = joinchannels('rgb',c,g*15,c)
%%
lbl_pre = label(logical(gephclean2),2);
lbl = GeneralAnalysis.findLabelsInMask(lbl_pre,dip_image(cellfillmask));
goodid = [];
wb = waitbar(0,'Labeling...');
goodlbl = lbl_out;
for ii = 1:max(lbl)
    % check that label is in every frame
    out = find(lbl==ii);
    [x, y, z] = ind2sub(size(lbl),out);
    Lia = ismember(1:size(lbl,3),z);
    if sum(Lia)~=size(lbl,3)
        goodlbl(lbl==ii) = false;
        continue;
    else
    goodid = [goodid;ii];
    end
     waitbar(ii/max(lbl),wb);
end
close(wb);
%%
recruit_traces = zeros(numel(goodid),size(geph,3));
wb = waitbar(0,'Calculating...');
for jj = 1:numel(goodid)
    gi = goodid(jj);
    F = squeeze(sum(geph,lbl_out==gi,[1 2])./sum(geph(:,:,0),lbl(:,:,0)==gi,[1 2]));
    recruit_traces(jj,:) = single(F);
    waitbar(jj/numel(goodid));
end
close(wb);
figure; plot(recruit_traces');
% norms = recruit_traces./recruit_traces(:,1);
% figure; plot(mean(norms)');
% figure; plot(norms');
%%
minval = 1.5;
gdtest = sum(recruit_traces>minval,2);
gdids = goodid(gdtest>0);
perim_mask = dip_image(zeros(size(lbl_out)));
wb = waitbar(0,'making perimeters...');
for ll = 1:size(gdids,1)
    incmask = lbl_out==gdids(ll);
    incmask = bdilation(incmask,1);
    perim =ga.maskperim(incmask);
    perim_mask = perim_mask | perim;
    waitbar(ll/size(gdids,1),wb);
end
close(wb);
% make perimeter for non-responders
noincmask = logical(goodlbl);
noincmask = bdilation(noincmask);
nonresp_perim = ga.maskperim(noincmask);
%%
cm4overlay = bone(256);
cm4overlay(1,:) = [0 0 1]; %set all perims
cm4overlay(end,:) = [1 0 0]; %set all responder perims
im = geph;
im(nonresp_perim) = 0;
im(perim_mask) = max(geph)*2;
h1 = dipshow(im,cm4overlay);
dipmapping(h1,'global');
dipmapping(h1,'lin');
dipmapping(h1,[100 1000]);
diptruesize(h1,150)

%%
clear gdtest
gdtest = max(recruit_traces,[],2);
[f,x] = ecdf(gdtest(gdtest<3));
figure; plot(x,f);
figure; h = histogram(gdtest(gdtest<3),70);
xlabel('Max Fold Intensity Increase (AU)','FontSize',14);
ylabel('Number of ROIs','FontSize',14);

