%% --- create an 'instance' of the class
sr = SRPipeline();

%% --- allow user to set parameters that you want to give them the ability to
sr.acq_nchannels = 2;                   % number of channels acquired in experiment
sr.TMF_do_TMF = 'true';                 % should TMF be performed? ('true' or 'false')
sr.TMF_save_TMF = 'yes';                % save the temporal median filtered data? options: 'no', 'tif', 'mat'. 'tif': saves each TMF corrected channel as a tif, 'mat': saves the TMF corrected channel in the output Matlab file    
sr.Tstorm_camera_pixelsize = 100;       % camera pixelsize in nm
sr.drift_onoff = 'on';                  % turn drift correction on/off, Important if drift has been corrected prior to running the pipeline, e_g. using CAT_app
sr.drift_segpara = 500;                 % segmentation parameters, how many frames are put together
sr.image_gen_generate_image = 'true';   % automatically generated images?

%% ---- run the analysis

%----this part can all also be done inside the class (as part of the AnalyzeDataFiles call) or can be broken out like
% this for flexibility---
sr.get_raw_data_files(); %prompt user to select files
sr.selectSaveDir(); % prompt user to select save dir
sr.startMIJI(); % open MIJI

%---- run the analysis!
sr.AnalyzeDataFiles();

%---- this can also be done inside the class or broken out
MIJ.exit;
%% -- some extra things: can call these at anytime to see params
sr.viewCameraparams();
sr.viewTMFparams();
sr.viewAnalysisparams();
sr.viewdriftparams();
sr.viewImageGenparams();
