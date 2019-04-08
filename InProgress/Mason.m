close all; clear all;
%% load the files
%--- select the raw tif file path for each channel
tfrch_path = 'G:\090718_3xCADmChGluA1_GFPSPl_TfRHT_ddsthrombin_n2_3_w3640.tif';
glua1_path = 'G:\090718_3xCADmChGluA1_GFPSPl_TfRHT_ddsthrombin_n2_3_w1561.tif';
synpod_path = 'G:\090718_3xCADmChGluA1_GFPSPl_TfRHT_ddsthrombin_n2_3_w2488.tif';
%--- load the files
tfr = loadtiff(tfrch_path);
synp = loadtiff(synpod_path);
glua1 = loadtiff(glua1_path);
%--- get max projection
tfr = max(tfr,[],3);
synp = max(synp,[],3);
glua1 = max(glua1,[],3);

%% smooth the TfR channel
tfr_g = GeneralAnalysis.imgLaplaceCutoff(tfr);
%% try automatic thresholding
test = multithresh(single(tfr_g),2);
tfr_mask = tfr_g>test(2);
% --- dilation of the mask step
tfr_mask = bdilation(tfr_mask,1);
%% select the good ones
good_tfr_mask = logical(GeneralAnalysis.cleanUpMaskKeepers_manual_square(tfr,tfr_mask));
%% get the info from the mask
msrquants = {'size','sum'};
tfr_info = measure(good_tfr_mask,tfr,msrquants);
synp_info = measure(good_tfr_mask,synp,msrquants);
glua1_info = measure(good_tfr_mask,glua1,msrquants);
t = table([1:size(tfr_info)]',tfr_info.size',synp_info.size',glua1_info.size',...
    tfr_info.sum',synp_info.sum',glua1_info.sum',...
    (tfr_info.sum./tfr_info.size)',(synp_info.sum./synp_info.size)',(glua1_info.sum./glua1_info.size)');
t.Properties.VariableNames = {'ROI','TfRMaskSize','SynpMaskSize','Glua1MaskSize',...
    'TfRMaskTotalIntensity', 'synpMaskTotalIntensity', 'Glua1MaskTotalIntensity',...
     'TfRMaskMeanIntensity','synpMaskMeanIntensity','glua1MaskMeanIntensity'};
%% write results to excel file - in directory where tfr data is
[FILEPATH,NAME,EXT] = fileparts(glua1_path);
savename= fullfile(FILEPATH,[NAME '_ResultsTable.xls']);
writetable(t,savename);
%% write ROIs as tiff files to overlay onto images 

