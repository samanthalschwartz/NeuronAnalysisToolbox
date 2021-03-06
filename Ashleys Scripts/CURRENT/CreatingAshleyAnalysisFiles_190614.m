%% Script to Identify CellFill and Recruitment Objects to make an 'AshleyFile' Object
% Update datafilepath variable below and then click the 'Run' button in the tool bar above. 
close all; clear all;
%-- set imaging parameters:
releaseframe = 11; % first frame of post release starts
%%
filename = 'C:\Users\bourkea\Dropbox\Data from Cloud\091319_GluA1_activitydependence\TTX\3_merge_concat+washout.tif';
[FILEPATH,NAME,EXT] = fileparts(filename);
temp = strsplit(NAME,'_');
savename = fullfile(FILEPATH,temp{1});

uiopen(filename); close all;

 %aa.convertfromload() %run this line of code if there are errors with
% SelectSoma, etc.

cellfill = image(:,:,:,1); 
TfR = image(:,:,:,2); 
cargo = image(:,:,:,3);

TfR = permute(TfR,[2 1 3]);
cellfill = permute(cellfill,[2 1 3]);
cargo = permute(cargo,[2 1 3]);


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
%%

% use last index to represent which stack to use
aa.cellFill.setimage(cellfill);
aa.surfaceCargo.setimage(cargo);
% aa.TfR.setimage(im_array(:,:,:,2));

[img_out_cf,sv_arr] = timedriftCorrect_parfor(aa.cellFill.image);

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
% IF CELL FILL MASK DOES NOT LOOK GOOD, DO NOT PROCEED!!!
aa.cellFill.lsig = [1 1 1];
aa.cellFill.gsig = [1 1 1];
aa.cellFill.mask_img; % two step masking:
% first image is to mask based on intensity (select dim axons in
% background, etc.)
% second image is to mask based on edge detection (select part of image
% edge, etc.)
aa.cellFill.viewMaskOverlayPerim;
% mask surface Cargo -- change these to modify image smoothing
lsig = [1 1 0];
gsig = [1 1 1];
aa.surfaceCargo.mask_img_highsens();
dipshow(aa.surfaceCargo.mask)
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
aa.cleanSurfaceCargoMask_RemoveRegion % call this line and press T to remove regions in a spatially-restricted manner

tempmask = aa.cleanedcargomask; % safety net in case the cleaned surface cargo mask is mistakenly reset
save(fullfile([savename '_TempCleanedSurfaceCargoMask']), 'tempmask');

%Remove the surface accumulations that last for only 1-2 frames (most
%likely junk)
minNumframes = 3; %keeps only the accumulations that last for 3+ frames
aa.removeShorterCargos(minNumframes,1:size(aa.cleanedcargomask,3)-minNumframes);
%run this so you can keep all the accumulations that occur in the last few
%frames (otherwise these would be falsely removed)

% now make the min distance image
% h = aa.plot_cargo_minFrame();
close all;
h = aa.plotCargoHeatMap(1);
% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
save(fullfile([savename '_AshleyFile.mat']), 'aa', '-v7.3'); 

saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
saveas(h,fullfile([savename '_timeHeatMap']),'eps');
close(h);