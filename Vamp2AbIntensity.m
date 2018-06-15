% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_5uL-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_2ndOnly647_5uL-1';
% % filenameFL = '180614FLtoxin_2ndOnly647_5uL-2_w1561_s2.TIF'
% filenameFL = '180614FLtoxin_2ndOnly647_5uL-1_w1561_s2.TIF'
% % filenameVamp = '180614FLtoxin_2ndOnly647_5uL-2_w2640_s2.TIF'
% filenameVamp = '180614FLtoxin_2ndOnly647_5uL-1_w2640_s2.TIF'
% fltoxin = GeneralAnalysis.loadtiff_1ch(fullfile(datadir, filenameFL));
% vamp2 = GeneralAnalysis.loadtiff_1ch(fullfile(datadir, filenameVamp));
% joinchannels('rgb',stretch(fltoxin)*.99,vamp2*0,stretch(vamp2))
%%
% load in vamp2 channel data
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_5uL-1';
% files = dir(fullfile(datadir,'*w2640*.TIF'));
% ii=1;
% im = GeneralAnalysis.loadtiff_1ch(fullfile(datadir,files(ii).name));
% imsum = sum(im,[],3); imsumg = gaussf(imsum);
% [mask,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(imsumg);
% mask_closing = bclosing(mask,2,-1)
% mask_out = GeneralAnalysis.bwmorph_timeseries(mask,'close',inf);
% bwmframe = bwmorph(single(mask),'fill',inf);
% bwmframe = bwmorph(single(bwmframe),'close',inf);
% bwmframe = bwmorph(single(bwmframe),'bridge',inf);
% bwmframe = bwmorph(single(bwmframe),'clean',inf);
% se = strel('disk',10);
% mask_closing = imclose(single(mask),se);
%%
% load in vamp2 channel data and select rois to analyze
datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614GFPcontrol_vamp2AB647_3ug-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_5uL-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_10uL-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_2ndOnly647_5uL-1';
files = dir(fullfile(datadir,'*w2640*.TIF'));
for ff=1:numel(files)
im = GeneralAnalysis.loadtiff_1ch(fullfile(datadir,files(ff).name));
imsum = sum(im,[],3); imsumg = gaussf(imsum);
h = dipshow(imsumg,'log');
filerois={};
fileroicords={};
ii=1; %index for user selected regions
while ishandle(h)
    try
        [filerois{ii}, fileroicords{ii}] = diproi(h);
        patch(fileroicords{ii}(:,1),fileroicords{ii}(:,2),'r');hold on;
        %                         roi = filerois{ii}; roicords = fileroicords{ii};
        %                         save(fullfile(savedirlist{ff},[NAME '_ROI_#' num2str(ii)]),'roicords','roi');
    catch ME % catch this issue from dipimage so that script can continue
        if isequal(ME.message,'You closed the window! That wasn''t the deal!')
            continue;
        end
    end
    ii=ii+1;
end
save(fullfile(datadir,[files(ff).name(1:end-4) '_ROIs']),'filerois','fileroicords');
end
%% load each mask in and give it a new number
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_5uL-1';
% files = dir(fullfile(datadir,'*.mat'));
% for ff = 3:numel(files)
%   curr = load(fullfile(datadir,files(ff).name));
% labeledim = curr.filerois{1}; 
% for rr = 2:numel(curr.filerois)
%    temp = curr.filerois{rr}*rr;
%    labeledim = labeledim + temp;
% end
% save(fullfile(datadir,[files(ff).name(1:end-4) '_labeled']),'labeledim');
% end
%% 
% load in each mask
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_5uL-1';
% roifiles = dir(fullfile(datadir,'*ROIs_labeled.mat'));
% imfiles = dir(fullfile(datadir,'*w2640*.TIF'));
% ii = 2;
% im = GeneralAnalysis.loadtiff_1ch(fullfile(datadir,imfiles(ii).name));
% sumim = sum(im,[],3);
% currrois = load(fullfile(datadir,roifiles(ii).name));
% labeledim = logical(labeledim);
% labeledim = label(logical(currrois.labeledim));
% msr = measure(currrois.labeledim,sumim,{'mean','mass','size'})
%  measure mean intensity
% measure mean intensity in other channel
% plot scatter plot
% or categorize based on min mean intensity 
%%
datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614GFPcontrol_vamp2AB647_3ug-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_2ndOnly647_5uL-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_5uL-1';
% datadir = 'G:\Sam\Data\180614 RBL vamp2Ab\180614FLtoxin_vamp2AB647_10uL-1';
roifiles = dir(fullfile(datadir,'*ROIs.mat'));
vamp_imfiles = dir(fullfile(datadir,'*w2640*.TIF'));
% FLtox_imfiles = dir(fullfile(datadir,'*w1561*.TIF'));
FLtox_imfiles = dir(fullfile(datadir,'*488LOW*.TIF'));

FLtoxPos_vampvals = [];
FLtoxNeg_vampvals = [];
figure; hold on; fa = gca;
cols = hsv(numel(roifiles));
for ii = 1:numel(roifiles)
    vampvals = [];
    FLtoxvals = [];
    vamp_imfile = GeneralAnalysis.loadtiff_1ch(fullfile(datadir,vamp_imfiles(ii).name));
    vamp_im = squeeze(sum(vamp_imfile,[],3));
    FLtox_imfile = GeneralAnalysis.loadtiff_1ch(fullfile(datadir,FLtox_imfiles(ii).name));
    FLtox_im = squeeze(sum(FLtox_imfile,[],3));
    
    curroi = load(fullfile(datadir,roifiles(ii).name));
    
    for mm = 1:numel(curroi.filerois)
        if isempty(curroi.filerois{mm})
            continue
        end
        labeledim = label(curroi.filerois{mm});
        vamp_msr = measure(labeledim,vamp_im,{'mean','mass','size'});
        FLtox_msr = measure(labeledim,FLtox_im,{'mean','mass','size'});
        vampvals = cat(1,vampvals,vamp_msr.mean);
        FLtoxvals = cat(1,FLtoxvals,FLtox_msr.mean);
        if FLtox_msr.mean>20000
            out = drawpolygon(FLtox_im,curroi.fileroicords{mm},max(FLtox_im)*2);
            dipshow(out,'lin');
            FLtoxPos_vampvals = cat(1,FLtoxPos_vampvals,vamp_msr.mean);
        else
            FLtoxNeg_vampvals = cat(1,FLtoxNeg_vampvals,vamp_msr.mean);
        end
    end
    plot(fa,vampvals,FLtoxvals,'.','MarkerSize',14,'Color',cols(ii,:)); hold on;
end
xlabel('Full Length Vamp2 Intensity'); ylabel('Toxin Expression');
% figure; plot(vampvals,FLtoxvals,'.','MarkerSize',14)

figure; histogram(FLtoxPos_vampvals,10); hold on;
histogram(FLtoxNeg_vampvals,10)