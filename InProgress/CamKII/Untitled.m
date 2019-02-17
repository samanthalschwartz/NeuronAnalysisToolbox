uiopen('G:\FromMicroscopeComputer\190118 CamKII virus test\300nL_I206KTD_+488\300nL_I206KTD_+488_w1488_s5_t.tiff',1);
% remove saturated pixels
maxval = prctile(image(:),99.99)*2;
newim = image;
newim(newim>maxval) = 0;
dipshow(newim);
% mask
img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(newim,1,1);
normim = img_laplcutoff./max(img_laplcutoff);
mask = GeneralAnalysis.imgThreshold_fixedUserInput(normim*newim);
GeneralAnalysis.viewMaskOverlay(newim,mask);

sumim = sum(mask,[],3)>0;
lbsum  = label(sumim);

sumim3d = repmat(sumim,1,1,size(image,3));
GeneralAnalysis.viewMaskOverlay(newim,sumim3d);
im561 = dip_image(image);
mask = repmat(sumim,1,1,size(im561,3));
intens_time = sum(im561.*mask,[],[1 2]);
vals = single(squeeze(intens_time./intens_time(1)));
figure; plot(vals)
