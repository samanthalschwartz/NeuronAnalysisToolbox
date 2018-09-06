%% Script to Identify CellFill and Recruitment Objects to make an 'AshleyFile' Object
% Update datafilepath variable below and then click the 'Run' button in the tool bar above. 
close all; clear all;
%-- set imaging parameters:
baselineframe_start = 1;
baselineframe_end = 6;
baselineframerate = 2;
releasetime = 1;
postreleaseframe_start = 7;
postreleaseframe_end = 'end';
postreleaseframerate = 2;

%% select - in order - the pre and post files you want to concat
files = uipickfiles('prompt','select the 2 images to concatenate',...
    'FilterSpec','Z:\Sam\MJK_zapERtrap_for_sam\AMB_globalrelease\041718_DHFR-GFP-NL1_TfRmChSEP_antiHA');
im_array1 = GeneralAnalysis.loadtiff_2ch(files{1});
im_array2 = GeneralAnalysis.loadtiff_2ch(files{2});
im_array = cat(3,im_array1,im_array2);
[FILEPATH,NAME,EXT] = fileparts(files{1});
temp = strsplit(NAME,'_');
savename = fullfile(FILEPATH,temp{1});

%%
% to view image series 'N' to go to next time slice 'P' for previous time
% slice. Use 'Mapping' tool bar at top to change the look up table.
% 'I' to zoom in or 'O' to zoom out.
aa = AshleyAnalysis();
aa.path_3ch_datafile = files;
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
aa.TfR = channelTfR();
%%
aa.imagingparams.baselineframe_start =baselineframe_start;
aa.imagingparams.baselineframe_end =baselineframe_end;
aa.imagingparams.baselineframerate =baselineframerate;
aa.imagingparams.releasetime =releasetime;
aa.imagingparams.postreleaseframe_start =postreleaseframe_start;
aa.imagingparams.postreleaseframe_end =postreleaseframe_end;
aa.imagingparams.postreleaseframerate =postreleaseframerate;
%%

% use last index to represent which stack to use
aa.cellFill.setimage(im_array(:,:,:,1));
aa.surfaceCargo.setimage(im_array(:,:,:,2));
% aa.TfR.setimage(im_array(:,:,:,2));

[img_out_cf,sv_arr] = GeneralAnalysis.timedriftCorrect(aa.cellFill.image);
aa.cellFill.setimage(img_out_cf);
img_out_sc = GeneralAnalysis.applydriftCorrect(aa.surfaceCargo.image,sv_arr);
aa.surfaceCargo.setimage(img_out_sc);
% img_out_tf = GeneralAnalysis.applydriftCorrect(aa.TfR.image,sv_arr);
% aa.TfR.setimage(img_out_tf);

% trim the image 
aa.cellFill.ROI_trim = [];
aa.cellFill.trim_rawimage()
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


% now save the object
% save(fullfile(datafilepath,[savename '_AshleyFile.mat']), 'aa'); 
save(fullfile([savename '_AshleyFile.mat']), 'aa'); 

% now make the min distance image
h = aa.plot_cargo_minFrame();
saveas(h,fullfile([savename '_timeHeatMap']),'fig');
saveas(h,fullfile([savename '_timeHeatMap']),'png');
close(h);