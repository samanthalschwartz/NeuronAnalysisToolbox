%% Script to Identify CellFill and Recruitment Objects to make an 'AshleyFile' Object

close all; clear all;
%-- set imaging parameters:

%%
filename = '\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\GLOBAL RELEASE\GluA1\051718\TIFF files\merges\slip2_2_merge_stitch.tiff';
releaseframe = 6;
[FILEPATH,NAME,EXT] = fileparts(filename);
temp = strsplit(NAME,'_,_');
savename = fullfile(FILEPATH,temp{1});

uiopen(filename); close all;

TfR = image(:,:,:,1); TfR = permute(TfR,[2 1 3]);
cellfill = image(:,:,:,2); cellfill = permute(cellfill,[2 1 3]);
cargo = image(:,:,:,3); cargo = permute(cargo,[2 1 3]);

%%
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
% aa.TfR.setimage(im_array(:,:,:,2));

[img_out_cf,sv_arr] = GeneralAnalysis.timedriftCorrect(aa.cellFill.image);
aa.cellFill.setimage(img_out_cf);
img_out_sc = GeneralAnalysis.applydriftCorrect(aa.surfaceCargo.image,sv_arr);
aa.surfaceCargo.setimage(img_out_sc);



% trim the image 
aa.cellFill.ROI_trim = [];
aa.cellFill.trim_rawimage();
aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;
aa.surfaceCargo.trim_rawimage;
aa.cleanedcargomask = []; %just make sure to reset cleaned cargomask
aa.imagingparams.releaseframe = releaseframe;
aa.cellFill.mask_img;
aa.surfaceCargo.mask_img_highsens;
% aa.surfaceCargo.viewMaskOverlayPerim;

% --- run this to look at cell fill mask compared to surface cargo mask
% realmask = aa.surfaceCargo.mask;
% aa.surfaceCargo.mask = aa.surfaceCargo.mask.*aa.cellFill.mask;
% perim = dt(aa.cellFill.mask);
% bin_im = (perim==1);
% aa.surfaceCargo.mask = aa.surfaceCargo.mask+bin_im;
% aa.surfaceCargo.viewMaskOverlayPerim;
% aa.surfaceCargo.mask = realmask; %--- DON'T FORGET TO RUN THIS LINE TO RESET THE SURFACE CARGO MASK ----
% 
% select cell soma for the file 
aa.cellFill.selectSoma();

% clean up the image:
%Call this line instead if you want to clean up by frame
% aa.cleanSurfaceCargoMaskbyFrame_Manual();


% now make the min distance image
%h = aa.plot_cargo_minFrame();
savename = filename(1:end-11);
close all;

aa.cargo_heatmap = [];
h = aa.plotCargoHeatMap;
% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
save(fullfile([savename '-AshleyFile.mat']), 'aa','-v7.3'); 
=======
h = aa.plotCargoHeatMap(1); % put 1 as input to reset
% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
save(fullfile([savename '_AshleyFile.mat']), 'aa', '-v7.3'); 
>>>>>>> d8f20e53a948e5ac7407e0cc743f5dbf71613ebf

saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
saveFigure_eps(h, fullfile([savename '_timeHeatMap']),'Arial');
close(h);