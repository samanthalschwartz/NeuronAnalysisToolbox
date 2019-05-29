uiopen('E:\Sample PV SyPhy mRuby\2.21.19\Slice 2 BS_XY1553295460_Z00_T0_C0.tif',1)
c0 = image;
uiopen('E:\Sample PV SyPhy mRuby\2.21.19\Slice 2 BS_XY1553295460_Z00_T0_C1.tif',1)
c1 = image;
uiopen('E:\Sample PV SyPhy mRuby\2.21.19\Slice 2 BS_XY1553295460_Z00_T0_C2.tif',1)
c2 = image;
% set images
emma.c0.image = c0;
emma.c1.image = c1;
emma.c2.image = c2;
% set mask cell fill
% [~,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(gaussf(c2(:,:,8)));
% maskc3= gaussf(c2,[1 1 0])>threshval;
emma.c2.mask = maskc3;
emma.c2.thresh = threshval;
% cleanmask = GeneralAnalysis.cleanUpMask_manual_square(dip_image(emma.c2.image),emma.c2.mask);
% emma.c2.mask = cleanmask;

% set ch1 mask  synaptophysin (red)
img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(c1);
% [mask1,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff(:,:,5))
threshval = 2.6900e+03;
emma.c1.mask = img_laplcutoff>threshval;
emma.c1.thresh = threshval;
% GeneralAnalysis.viewMaskOverlay(c1,emma.c1.mask )


% set ch0 mask PV cell fill (green)
img_laplcutoff0 = GeneralAnalysis.imgLaplaceCutoff(c0);
% [mask0,threshval0] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff0(:,:,5))
threshval0 = 425.7595;
emma.c0.mask = img_laplcutoff0>threshval0;
emma.c0.thresh = threshval0;

%%
GeneralAnalysis.viewMaskOverlay(c0,emma.c0.mask)
joinchannels('rgb',stretch(emma.c1.mask.*emma.c2.mask), stretch(emma.c0.mask.*emma.c2.mask), stretch(emma.c2.mask))


vd = viewer3D();
vd.numchannels = 3;
vd.image{1}  = dip_image(emma.c2.mask);
vd.image{2}  = dip_image(emma.c1.mask.*emma.c2.mask);
vd.image{3}  = dip_image(emma.c0.mask.*emma.c2.mask);
% vd.image{2}  = maskch1;
% vd.image{3}  = maskch0;
vd.scaleval{1} = 10;
vd.scaleval{2} = 10;
vd.scaleval{3} = 10;

vd.colors(2,:) = [1 0 0];
vd.colors(3,:) = [0 1 0];
vd.dimension = [1 1 0.3];
vd.initialize3Dimage

set(gca,'XLim',[210 320])
set(gca,'ZLim',[5 50])
%% making just the cell fill

%%
axis vis3d
clear F G
firstview = [100 12]; view(firstview);[az,el] = view;
lastview = [99 48];
loops1 = 50;
steps = (lastview-firstview)/loops1
F(loops1) = struct('cdata',[],'colormap',[]);
for j = 1:loops1
    view(az+((j-1).*steps(1)),el+((j-1).*steps(2)));
    drawnow
    F(j) = getframe(gcf);
end
G(loops1) = struct('cdata',[],'colormap',[]);
steps = (lastview-firstview)/loops1;
[az,el] = view;
for j = 1:loops1
    view(az-((j-1).*steps(1)),el-((j-1).*steps(2)));
    drawnow
    G(j) = getframe(gcf);
end
filesavename = 'E:\Sample PV SyPhy mRuby\Slice 2 BS_XY1553295460_Z00_T0_c2only'
myVideo = VideoWriter(filesavename,'MPEG-4');
myVideo.FrameRate = 15;  % Default 30
myVideo.Quality = 100;
open(myVideo);
test= cat(2,F,G);
writeVideo(myVideo, test);
close(myVideo);

%%
vd1 = viewer3D();
vd1.numchannels = 1;
vd1.image{1}  = dip_image(emma.c2.mask);

% vd.image{2}  = maskch1;
% vd.image{3}  = maskch0;
vd1.scaleval{1} = 10;
vd1.dimension = [1 1 0.3];

vd1.colors(1,:) = [1 1 1];
vd1.initialize3Dimage

%%
% img_laplcutoff2 = GeneralAnalysis.imgLaplaceCutoff(gaussf(c2,[1 1 0]));
% [mask,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff2(:,:,9));
% [mask,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(gaussf(c2(:,:,9)));
% dog = abs(gaussf(c2,[4 4 0]) - gaussf(c2,[2 2 0]));
%  dog_cutoff = -dog;
%         dog_cutoff(dog_cutoff<0) = 0;



% 
%  img_dcccutoff = imgDccCutoff(img_in,gsig)
% 
% img_dggcutoff = GeneralAnalysis.imgDccCutoff(c2,[1 1 0])
% 
% 
% mask = gaussf(c2,[1 1 0])>threshval
% mask = img_laplcutoff2>threshval
% cleanmask = GeneralAnalysis.cleanUpMask_byframe_square(dip_image(c2),mask);
% mask2 = logical(cleanmask);
% 