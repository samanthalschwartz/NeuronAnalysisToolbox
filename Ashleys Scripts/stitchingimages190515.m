close all; clear all
filename = 'G:\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\060118\TIFF files\merges\slip1_2_merge.tif';
uiopen(filename);
soma = image;
close all;
uiopen('G:\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\060118\TIFF files\merges\slip1_2_dendrites_merge.tif');
dendrites = image;


close all;
[ch1stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dendrites(:,:,:,2),soma(:,:,:,2));
dipshow(ch1stitchimage,'log')

[ch2stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dendrites(:,:,:,1),soma(:,:,:,1),ccpeak);

[ch3stitchimage, ccpeak] = GeneralAnalysis.stitch2movies(dendrites(:,:,:,3),soma(:,:,:,3),ccpeak);

stitchimage = cat(4,ch1stitchimage,ch2stitchimage,ch3stitchimage);
save([filename(1:end-4) '_stitch'],'stitchimage');
LibTiff(stitchimage,[filename(1:end-4) '_stitch']);

aa = AshleyAnalysis();
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
% aa.TfR = channelTfR();
aa.cellFill.setimage(ch1stitchimage);
aa.surfaceCargo.setimage(ch3stitchimage);
% aa.TfR.setimage(ch2stitchimage);
aa.path_channel_cellfill = [filename(1:end-4) '_stitch'];
save([filename(1:end-4) '_stitch_AshleyFile'],'aa');







