dend_filename = 'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_dendrites_merge.tif';
uiopen(dend_filename); close all;
tfr_dend = image(:,:,:,1); tfr_dend = permute(tfr_dend,[2 1 3]);
cellfill_dend = image(:,:,:,2); cellfill_dend = permute(cellfill_dend,[2 1 3]);
cargo_dend = image(:,:,:,3); cargo_dend = permute(cargo_dend,[2 1 3]);

soma_filename = 'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_merge.tif';
uiopen(soma_filename); close all;
tfr_soma = image(:,:,:,1); tfr_soma = permute(tfr_soma,[2 1 3]);
cellfill_soma = image(:,:,:,2); cellfill_soma = permute(cellfill_soma,[2 1 3]);
cargo_soma = image(:,:,:,3); cargo_soma = permute(cargo_soma,[2 1 3]);

% check that stitching works for cellfill
ccpeak = {235,980};
soma = squeeze(sum(cellfill_soma,3));
dendrites = squeeze(sum(cellfill_dend,3));
[stitchimage, ccpeak] = GeneralAnalysis.stitch2images(dendrites,soma,ccpeak,1);
dipshow(stitchimage);
% now make all the stitched images
ccpeakmov = cell(size(cellfill_dend,3),1);
ccpeakmov(:) = {ccpeak};
% cellfill
[stitchmovie_cellfill, ccpeak] = GeneralAnalysis.stitch2movies(cellfill_dend,cellfill_soma,ccpeakmov);
GeneralAnalysis.LibTiff(stitchmovie_cellfill,'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_stitch_cellfill.tif')
% cargo
[stitchmovie_cargo, ccpeak] = GeneralAnalysis.stitch2movies(cargo_dend,cargo_soma,ccpeakmov);
GeneralAnalysis.LibTiff(stitchmovie_cargo,'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_stitch_cargo.tif');
% TfR
[stitchmovie_tfr, ccpeak] = GeneralAnalysis.stitch2movies(tfr_dend,tfr_soma,ccpeakmov);
GeneralAnalysis.LibTiff(stitchmovie_tfr,'G:\zapERtrap\old\AMB_previous\GluA1release_REs\060118\merges\slip2\1_stitch_tfr.tif');
