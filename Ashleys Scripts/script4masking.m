%% Script to Identify CellFill and Recruitment Objects to make an 'AshleyFile' Object
% Update datafilepath variable below and then click the 'Run' button in the tool bar above. 
close all; clear all;
%% add the correct path to your .tiff file data
datafilepath = 'Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\GluA1release_REs\060118\merges\slip2\2_dendrites_merge.tif';
%%
% to view image series 'N' to go to next time slice 'P' for previous time
% slice. Use 'Mapping' tool bar at top to change the look up table.
% 'I' to zoom in or 'O' to zoom out.
addpath(genpath('Z:\Lab Resources\Analysis Resources\Matlab Resource\NeuronAnalysisToolBox'));
aa = AshleyAnalysis();
% aa.path_channel_cellfill = fullfile(datafilepath,'cell1_r.tif');
% aa.path_channel_surfaceCargo = fullfile(datafilepath,'cell1_fr.tif');
% aa.path_channel_TfR = fullfile(datafilepath,'cell1_g.tif');
aa.loadImages;
aa.path_3ch_datafile = datafilepath;
aa.cellFill = channelCellFill();
aa.surfaceCargo = channelSurfaceCargo();
aa.TfR = channelTfR();

% use last index to represent which stack to use
im_array = GeneralAnalysis.loadtiff_3ch(datafilepath);
aa.cellFill.setimage(im_array(:,:,:,1));
aa.surfaceCargo.setimage(im_array(:,:,:,3));
aa.TfR.setimage(im_array(:,:,:,2));
% trim the image 
aa.cellFill.ROI_trim = [];
aa.cellFill.trim_rawimage()
aa.surfaceCargo.ROI_trim = aa.cellFill.ROI_trim;
aa.surfaceCargo.trim_rawimage;
aa.TfR.ROI_trim = aa.cellFill.ROI_trim;
aa.TfR.trim_rawimage;

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
aa.cellFill.lsig = [0 0 2];
aa.cellFill.gsig = [0 0 2];
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

% now save the object
save([datafilepath(1:end-4) '_AshleyFile.mat'], 'aa'); 

% now make the min distance image
h = aa.plot_cargo_minFrame();
saveas(h,fullfile([datafilepath(1:end-4) '_timeHeatMap']),'fig');
saveas(h,fullfile([datafilepath(1:end-4) '_timeHeatMap']),'png');
close(h);