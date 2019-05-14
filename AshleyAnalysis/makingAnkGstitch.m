filename = 'F:\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\051418\TIFF files\merges\4_stitch_merge_AshleyFile.tif';
uiopen(filename); close all;

cellfill = image(:,:,:,1); cellfill = permute(cellfill,[2 1 3]);
tfr = image(:,:,:,2); tfr = permute(tfr,[2 1 3]);
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
aa.surfaceCargo.trim_rawimage;
aa.maskImages;
% aa.surfaceCargo.removeBoundaryMaskArtifact;
% aa.cleanCellFillMask_Manual;
% aa.cleanSurfaceCargoMask;
aa.cleanSurfaceCargoMask_Manual;
aa.calc_cargo_minFrame
aa.plotCargoHeatMap
%%
filename_cellfill = 'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_stitch_cellfill.tif.tiff';
filename_cargo = 'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_stitch_cargo.tif.tiff';
filename_tfr = 'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_stitch_tfr.tif.tiff';

uiopen(filename_cellfill); close all;
cellfill = squeeze(image); %cellfill = permute(cellfill,[2 1 3]);
uiopen(filename_tfr); close all;
tfr = squeeze(image); %tfr = permute(tfr,[2 1 3]);
uiopen(filename_cargo); close all;
cargo = squeeze(image); %cargo = permute(cargo,[2 1 3]);

aa = AshleyAnalysis();
aa.path_channel_cellfill = filename_cellfill;
aa.path_channel_surfaceCargo = filename_cargo;
aa.path_channel_TfR = filename_tfr;
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
aa.TfR = channelTfR();
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
% aa.TfR.setimage(im_array(:,:,:,2));
aa.cellFill.trim_rawimage;
aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;
aa.surfaceCargo.trim_rawimage;
aa.maskImages;
% aa.surfaceCargo.removeBoundaryMaskArtifact;
% aa.cleanCellFillMask_Manual;
% aa.cleanSurfaceCargoMask;
aa.cleanSurfaceCargoMask_Manual;
aa.calc_cargo_minFrame
aa.plotCargoHeatMap
