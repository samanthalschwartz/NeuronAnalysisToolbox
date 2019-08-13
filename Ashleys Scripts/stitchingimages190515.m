close all; clear all
filename = '\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\ACTIVITY DEPENDENCE\062519_Matt_Activity_NL1_local\slip2_bic\1_merge.tif';
uiopen(filename);
soma = image;
close all;
uiopen('\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\ACTIVITY DEPENDENCE\062519_Matt_Activity_NL1_local\slip2_bic\1_merge_2.tif');
dendrites = image;


close all;


 [dend_timedrift_ch1,svarr] = GeneralAnalysis.timedriftCorrect(dendrites(:,:,:,1));
 [dend_timedrift_ch3] = GeneralAnalysis.applydriftCorrect(dip_image(dendrites(:,:,:,3)),svarr);
% 
 [soma_timedrift_ch1,svarr1] = GeneralAnalysis.timedriftCorrect(soma(:,:,:,1));
 [soma_timedrift_ch3] = GeneralAnalysis.applydriftCorrect(dip_image(soma(:,:,:,3)),svarr1);
% 
 [ch1stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dend_timedrift_ch1,soma_timedrift_ch1);
 dipshow(ch1stitchimage,'log')



[ch1stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dendrites(:,:,:,1),soma(:,:,:,1),ccpeak);
dipshow(ch1stitchimage,'log')

[ch3stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dendrites(:,:,:,3),soma(:,:,:,3),ccpeak);
figure
dipshow(ch3stitchimage,'log')

[ch2stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dendrites(:,:,:,2),soma(:,:,:,2),ccpeak);
% stitchimage = cat(4,ch1stitchimage,ch2stitchimage,ch3stitchimage);
stitchimage = cat(4,ch1stitchimage,ch2stitchimage,ch3stitchimage);
save([filename(1:end-4) '_stitch'],'stitchimage');

[ch1stitchimage_timedrift,svarr] = GeneralAnalysis.timedriftCorrect(ch1stitchimage);
[ch2stitchimage_timedrift] = GeneralAnalysis.applydriftCorrect(dip_image(ch2stitchimage),svarr);
[ch3stitchimage_timedrift] = GeneralAnalysis.applydriftCorrect(dip_image(ch3stitchimage),svarr);
stitchimage = cat(4,ch1stitchimage_timedrift,ch3stitchimage_timedrift);
save([filename(1:end-4) '_stitch'],'stitchimage');
LibTiff(ch1stitchimage_timedrift,[filename(1:end-4) '_ch1_stitch']);
LibTiff(ch2stitchimage_timedrift,[filename(1:end-4) '_ch2_stitch']);
LibTiff(ch3stitchimage_timedrift,[filename(1:end-4) '_ch3_stitch']);

aa = AshleyAnalysis();
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
% aa.TfR = channelTfR();
aa.cellFill.setimage(ch1stitchimage_timedrift);
aa.surfaceCargo.setimage(ch3stitchimage_timedrift);
% aa.TfR.setimage(ch2stitchimage);
aa.path_channel_cellfill = [filename(1:end-4) '_stitch'];
save([filename(1:end-4) '_stitch_AshleyFile'],'aa');







