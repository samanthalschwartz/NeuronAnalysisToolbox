close all; clear all
datadir = 'G:\FromMicroscopeComputer\190316 pHujiIntensityTesting\FRAP\488FRAP\cell1_488FRAP_pHujiGaba_Cry2OligGephIB__20190316_84852 PM';
savedir = datadir;
%%
ch488_filebase = 'cell1_488FRAP_pHujiGaba_Cry2OligGephIB__w0000_z*';
ch488files = dir2cell(datadir,ch488_filebase);
clear im_488;
im_488 = [];
for ii = 1:numel(ch488files)
image = loadtiff(ch488files{ii});
im_488 = cat(4,im_488,image);
end
image_488 = max(dip_image(im_488),[],4);
GeneralAnalysis.LibTiff(image_488,fullfile(savedir,ch488_filebase(1:end-1)));
%%
ch561_filebase = 'cell1_488FRAP_pHujiGaba_Cry2OligGephIB__w0001_z*';
ch561files = dir2cell(datadir,ch561_filebase);
clear im_561;
im_561 = [];
for ii = 1:numel(ch561files)
image = loadtiff(ch561files{ii});
im_561 = cat(4,im_561,image);
end
image_561 = max(dip_image(im_561),[],4);

GeneralAnalysis.LibTiff(image_561,fullfile(savedir,ch561_filebase(1:end-1)));

