close all; clear all
num_base = 4;
num_repeats = 6;
repeatlen_561only = 8;
repeatlen_561_488 = 1;
repeatlength_tot = repeatlen_561only+repeatlen_561_488;
totlen = (num_base+num_repeats)*repeatlength_tot;
real_totlen = num_base+repeatlength_tot*num_repeats;
datadir = 'G:\FromMicroscopeComputer\190303 mScarGeph FRAP_FingR\Cry2Olig_GephIB\cell1_mScarGeph_Cry2OligGephIBGFP_20190303_92949 PM';
savedir = datadir;
% cd([datadir '/..']);
% savedir = pwd;
%%
ch561_filebase = 'cell1_mScarGeph_Cry2OligGephIBGFP_w0001_z*';
ch561files = dir2cell(datadir,ch561_filebase);
clear im_561;
for ii = 1:numel(ch561files)
    clear im_561frame im_base561 im_post561 image dimage
    image = loadtiff(ch561files{ii});
    dimage_561 = dip_image(image);
    if ii == 1
        im_561 = zeros(size(dimage_561,1),size(dimage_561,2),real_totlen,numel(ch561files));
    end
    assert(size(image,3) == totlen);
    im_base561 = dimage_561(:,:,0:repeatlength_tot:(num_base-1)*repeatlength_tot);
    im_post561 = dimage_561(:,:,(num_base)*repeatlength_tot:end);
    im_561frame = cat(3,im_base561,im_post561);
    im_561(:,:,:,ii) = im_561frame;
end
image_561 = max(dip_image(im_561),[],4);
GeneralAnalysis.LibTiff(image_561,fullfile(savedir,ch561_filebase(1:end-1)));
%% this section is wrong!!!
ch488_filebase = 'cell1_mScarGeph_Cry2OligGephIBGFP_w0000_z*';
ch488files = dir2cell(datadir,ch488_filebase);
clear im_488;
% im_post488_b = dimage_488(:,:,(num_base+1)*repeatlength_tot:(repeatlength_tot):end); % this is weird crap image over and over
for ii = 1:numel(ch488files)
clear im_488frame im_base488 im_post488 image dimage
    image = loadtiff(ch488files{ii});
    dimage_488 = dip_image(image);
    if ii == 1
        im_488 = zeros(size(dimage_488,1),size(dimage_488,2),real_totlen,numel(ch488files));
    end
assert(size(image,3) == totlen)
im_base488 = dimage_488(:,:,0:repeatlength_tot:(num_base-1)*repeatlength_tot);
im_post488 = dimage_488(:,:,(num_base+1)*repeatlength_tot-1:(repeatlength_tot):end);
im_post488rep = repmat(im_post488,[1 1 1 repeatlength_tot]);
im_post488rep=reshape(im_post488rep,[size(im_post488,1),size(im_post488,2),(size(im_post488,3)*repeatlength_tot)]);
im_488frame = cat(3,im_base488,im_post488rep);
im_488(:,:,:,ii) = im_488frame;
end
image_488 = max(dip_image(im_488),[],4);
GeneralAnalysis.LibTiff(image_488,fullfile(savedir,ch488_filebase(1:end-1)));



    
    





