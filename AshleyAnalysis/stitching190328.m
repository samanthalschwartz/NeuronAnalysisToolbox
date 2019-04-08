uiopen('G:\zapERtrap\AMB_previous\051318 NL1 insertion\3_merge.tif',1)
fillsoma = image(:,:,:,1);
dipshow(fillsoma)
othersoma = image(:,:,:,2);
dipshow(othersoma)
cargosoma = image(:,:,:,3);
dipshow(cargosoma)

uiopen('G:\zapERtrap\AMB_previous\051318 NL1 insertion\3_dendrites_merge.tif',1);
filldends = image(:,:,:,1);
dipshow(filldends)
otherdends = image(:,:,:,2);
dipshow(otherdends)
cargodends = image(:,:,:,3);
dipshow(cargodends)

[stitchsoma, ccpeak] = GeneralAnalysis.stitch2movies(filldends,fillsoma);
[stitchsoma, ccpeak] = GeneralAnalysis.stitch2images(filldends,fillsoma);
[stitchother, ~] = GeneralAnalysis.stitch2images(otherdends,othersoma,ccpeak);
[stitchcargo, ~] = GeneralAnalysis.stitch2images(cargodends,cargosoma,ccpeak);

LibTiff(stitchsoma,'G:\zapERtrap\AMB_previous\051318 NL1 insertion\3_fullcell_soma.tif')
LibTiff(stitchother,'G:\zapERtrap\AMB_previous\051318 NL1 insertion\3_fullcell_other.tif')
LibTiff(stitchcargo,'G:\zapERtrap\AMB_previous\051318 NL1 insertion\3_fullcell_cargo.tif')