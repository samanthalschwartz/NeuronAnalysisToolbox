close all; clear all
datadir = 'F:\FromMicroscopeComputer\190318 HaloGabaTesting\488bleach\cell5_HaloGaba647_GephIB_488FRAP_20190319_64825 PM';
savedir = datadir;
% cd([datadir '/..']);
% savedir = pwd;
%%
ch561_filebase = 'cell5_HaloGaba647_GephIB_488FRAP_w0001_z*';
ch561files = dir2cell(datadir,ch561_filebase);
im_561 = [];
for ii = 1:numel(ch561files)
    image = loadtiff(ch561files{ii});
    im_561(:,:,:,ii) = image;
end
image_561 = max(dip_image(im_561),[],4);
GeneralAnalysis.LibTiff(image_561,fullfile(savedir,ch561_filebase(1:end-1)));
%%
ch488_filebase = 'cell5_HaloGaba647_GephIB_488FRAP_w0000_z*';
ch488files = dir2cell(datadir,ch488_filebase);
im_488 = [];
for ii = 1:numel(ch488files)
    image = loadtiff(ch488files{ii});
    im_488(:,:,:,ii) = image;
end
image_488 = max(dip_image(im_488),[],4);
GeneralAnalysis.LibTiff(image_488,fullfile(savedir,ch488_filebase(1:end-1)));