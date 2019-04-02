uiopen('F:\zapERtrap\AMB_previous\051318 NL1 insertion\2test_fullcell_merge.tif',1);
cellfill = image(:,:,:,1);
tfr = image(:,:,:,2);
cargo = image(:,:,:,3);
%%

aa = AshleyAnalysis();
aa.path_3ch_datafile = 'F:\zapERtrap\AMB_previous\051318 NL1 insertion\2test_fullcell_merge.tif';
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
aa.TfR = channelTfR();
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
% aa.TfR.setimage(im_array(:,:,:,2));

aa.cellFill.mask_img;
aa.cellFill.viewMaskOverlayPerim;
aa.surfaceCargo.lsig = [1 1 0];
aa.surfaceCargo.gsig = [1 1 1];
aa.surfaceCargo.mask_img_highsens

img_m = medif(obj.image,3);
           img_mg = gaussf(aa.surfaceCargo.image,aa.surfaceCargo.gsig);
           img_laplcutoff = GeneralAnalysis.imgLaplaceCutoff(img_mg,aa.surfaceCargo.lsig,aa.surfaceCargo.gsig);
           %            glim = gaussf(img_laplcutoff);
           %            obj.mask = GeneralAnalysis.imgThreshold(img_laplcutoff);
           [maskpre,threshval] = GeneralAnalysis.imgThreshold_fixedUserInput(img_laplcutoff);
           %            maskpre = GeneralAnalysis.imgThreshold(glim.^1.2);
                     
           mask = GeneralAnalysis.bwmorph_timeseries(maskpre,'thicken',1);
           ll = slice_op('watershed',-img_mg,2);
           maskws = mask;
           maskws(ll) = 0;
           
           aa.surfaceCargo.mask = logical(maskws);
%            obj.mask = logical(mask);
           aa.surfaceCargo.backgroundvalue = threshval;




aa.surfaceCargo.viewMaskOverlayPerim;

aa.cellFill.selectSoma();

% clean up the image:
uiwait(msgbox('Select regions in the mask to remove. Once you are satisfied, close the window.','Clean UP','modal'));
aa.cleanSurfaceCargoMask_Manual();
h = aa.plot_cargo_minFrame();
M = aa.plotDensityperTime([18,210,400]);
