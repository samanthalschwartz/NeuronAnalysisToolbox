%% Script to Identify CellFill and Recruitment Objects to make an 'AshleyFile' Object
% Update datafilepath variable below and then click the 'Run' button in the tool bar above. 
close all; clear all;
%-- set imaging parameters:
baselineframe_start = 1; % first frame number that baseline acquisition begins
baselineframe_end = 11; % last frame number of baseline
baselineframerate = 1; % frame rate in minutes/frame 
releasetime = 2; % time in minutes, after release, that first frame of post release starts
postrelease(1).frame_start = 12; % first frame number of post release
postrelease(1).frame_end = 'end'; % last frame number of post release - or 'end' if post release goes until end of series
postrelease(1).framerate = 2; % frame rate in minutes/frame
%postrelease(2).frame_start = 21;
%postrelease(2).frame_end = 'end';
%postrelease(2).framerate = 1;
%%
filename = 'Y:\Lab Projects\zapERtrap\Raw Data\AISexpts\071618\NL1\TIFF files\071618_DHFR-GFP-NL1_AnkG-mCherry_antiGFP640_global_8_stitched.tif';
[FILEPATH,NAME,EXT] = fileparts(filename);
temp = strsplit(NAME,'_');
savename = fullfile(FILEPATH,temp{1});

uiopen(filename); close all;

TfR = image(:,:,:,1); TfR = permute(TfR,[2 1 3]);
cellfill = image(:,:,:,2); cellfill = permute(cellfill,[2 1 3]);
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
aa.imagingparams.baselineframe_start =baselineframe_start;
aa.imagingparams.baselineframe_end =baselineframe_end;
aa.imagingparams.baselineframerate =baselineframerate;
aa.imagingparams.releaseframe = releasetime;
aa.imagingparams.postrelease = postrelease;
%%

% use last index to represent which stack to use
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
% aa.TfR.setimage(im_array(:,:,:,2));

[img_out_cf,sv_arr] = GeneralAnalysis.timedriftCorrect(aa.cellFill.image);
aa.cellFill.setimage(img_out_cf);
img_out_sc = GeneralAnalysis.applydriftCorrect(aa.surfaceCargo.image,sv_arr);
aa.surfaceCargo.setimage(img_out_sc);
% img_out_tf = GeneralAnalysis.applydriftCorrect(aa.TfR.image,sv_arr);
% aa.TfR.setimage(img_out_tf);

% trim the image 
aa.cellFill.ROI_trim = [];
aa.cellFill.trim_rawimage();
aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;
aa.surfaceCargo.trim_rawimage;
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
aa.cellFill.lsig = [1 1 1];
aa.cellFill.gsig = [1 1 1];
aa.cellFill.mask_img;
aa.cellFill.viewMaskOverlayPerim;
% mask surface Cargo -- change these to modify image smoothing
aa.surfaceCargo.lsig = [1 1 0];
aa.surfaceCargo.gsig = [1 1 1];
aa.surfaceCargo.mask_img_highsens
aa.surfaceCargo.viewMaskOverlayPerim;

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
uiwait(msgbox('Select regions in the mask to remove. Once you are satisfied, close the window.','Clean UP','modal'));
aa.cleanSurfaceCargoMask_Manual();
% aa.cleanSurfaceCargoMask_Manual(1); % call this line instead if you want to start again
aa.cleanSurfaceCargoMaskbyFrame_Manual();

%To view cleanedcargomask overlayed on top of cellfill
GeneralAnalysis.viewMaskOverlay(aa.cellFill.image,aa.cleanedcargomask)

% now make the min distance image
% h = aa.plot_cargo_minFrame();
close all;
aa.cargo_heatmap = [];
imgparam.maxtime = 60;
h = aa.plotCargoHeatMap(1,imgparam);
% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
save(fullfile([savename '_AshleyFile.mat']), 'aa'); 

saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
saveas(h,fullfile([savename '_timeHeatMap']),'eps');
close(h);