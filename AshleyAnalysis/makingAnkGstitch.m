filename = 'Y:\Lab Projects\zapERtrap\Raw Data\AISexpts\071618\NL1\TIFF files\071618_DHFR-GFP-NL1_AnkG-mCherry_antiGFP640_global_5_stitched.tif';
uiopen(filename); close all;

ankG = image(:,:,:,1); ankG = permute(ankG,[2 1 3]);
cellfill = image(:,:,:,2); cellfill = permute(cellfill,[2 1 3]);
cargo = image(:,:,:,3); cargo = permute(cargo,[2 1 3]);

aa = AshleyAnalysis();
aa.path_3ch_datafile = filename;
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
aa.TfR = channelTfR();
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
% aa.TfR.setimage(im_array(:,:,:,2));
aa.cellFill.trim_rawimage;
aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;

aa.maskImages;

aa.surfaceCargo.removeBoundaryMaskArtifact;

aa.cleanCellFillMask_Manual;
aa.cleanSurfaceCargoMask;
aa.cleanSurfaceCargoMask_Manual;
joinchannels('rgb',aa.cellFill.mask,aa.cleanedcargomask)
