%% Script to Identify CellFill and Recruitment Objects to make an 'AshleyFile' Object
% Update datafilepath variable below and then click the 'Run' button in the tool bar above. 

%postrelease(2).frame_start = 21;
%postrelease(2).frame_end = 'end';
%postrelease(2).framerate = 1;
% %%
% filename = '\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\GLOBAL RELEASE\NL1\051318\TIFF files\3_stitched.tif';
% [FILEPATH,NAME,EXT] = fileparts(filename);
% temp = strsplit(NAME,'_');
% savename = fullfile(FILEPATH,temp{1});
% 
% uiopen(filename); close all;
% 
% TfR = image(:,:,:,2); TfR = permute(TfR,[2 1 3]);
% cellfill = image(:,:,:,1); cellfill = permute(cellfill,[2 1 3]);
% cargo = image(:,:,:,3); cargo = permute(cargo,[2 1 3]);
% 
close all; clear all;
%-- set imaging parameters:
releaseframe = 12;
postrelease.framerate =2; %in min/frame
%%
filename = '\\data\dept\SOM\PHARM\All\Research\KennedyLab\Lab Projects\zapERtrap\Raw Data\GLOBAL RELEASE\NL1\050118\TIFF files\cell3_stitched.tif';
[FILEPATH,NAME,EXT] = fileparts(filename);
temp = strsplit(NAME,'_,_');
savename = fullfile(FILEPATH,temp{1});

uiopen(filename); close all;

TfR = image(:,:,:,2); TfR = permute(TfR,[2 1 3]);
cellfill = image(:,:,:,1); cellfill = permute(cellfill,[2 1 3]);
cargo = image(:,:,:,3); cargo = permute(cargo,[2 1 3]);

%%
% to view image series 'N' to go to next time slice 'P' for previous time
% slice. Use 'Mapping' tool bar at top to change the look up table.
% 'I' to zoom in or 'O' to zoom out.
aa = AshleyAnalysis();
aa.path_3ch_datafile = filename;
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
aa.TfR = channelTfR();
%%
aa.imagingparams.releaseframe = releaseframe;
aa.imagingparams.postrelease.framerate = postrelease.framerate;
%%

% use last index to represent which stack to use
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
aa.TfR.setimage(TfR);
% 
[img_out_cf,sv_arr] = GeneralAnalysis.timedriftCorrect_parfor(aa.cellFill.image);
aa.cellFill.setimage(img_out_cf);
img_out_sc = GeneralAnalysis.applydriftCorrect(aa.surfaceCargo.image,sv_arr);
aa.surfaceCargo.setimage(img_out_sc);
img_out_tf = GeneralAnalysis.applydriftCorrect(aa.TfR.image,sv_arr);
aa.TfR.setimage(img_out_tf);
%%
% trim the image 
aa.cellFill.ROI_trim = [];
aa.cellFill.trim_rawimage();
aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;
aa.surfaceCargo.trim_rawimage;
aa.cleanedcargomask = []; %just make sure to reset cleaned cargomask
% aa.TfR.ROI_trim = aa.cellFill.ROI_trim;
% aa.TfR.trim_rawimage;

% % to see transferin channel
% h = dipshow(aa.TfR.image);
% dipmapping(h,'lin')
% 
% % to see cellFill channel
% h = dipshow(aa.cellFill.image);
% dipmapping(h,'lin')
% % to see surface Cargo channel
% h = dipshow(aa.surfaceCargo.image);
% dipmapping(h,'lin')

% mask cell -- change these to modify image smoothing
% aa.cellFill.lsig = [1 1 1];
% aa.cellFill.gsig = [1 1 1];
aa.cellFill.mask_img;
% mask surface Cargo -- change these to modify image smoothing
% aa.surfaceCargo.lsig = [1 1 0];
% aa.surfaceCargo.gsig = [1 1 1];
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
aa.cleanSurfaceCargoMaskbyFrame_Manual();
%aa.cleanSurfaceCargoMask_Manual();


% now make the min distance image
%h = aa.plot_cargo_minFrame();
savename = filename(1:end-4);
close all;

aa.cargo_heatmap = [];
% h = aa.plotCargoHeatMap;
% now save the object
imgparam.maxtime = 120;
% imgparam.colormap = '';
h = aa.plotCargoHeatMap(1,imgparam); % put 1 as input to reset
% now save the object
fullsavename = fullfile([savename '_AshleyFile.mat']);
aa.save(fullsavename); 

saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
saveFigure_eps(h, fullfile([savename '_timeHeatMap']),'Arial');
close(h);