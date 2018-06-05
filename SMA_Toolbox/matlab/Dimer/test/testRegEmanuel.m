%
% Analysis of Emanuel's 2015/9/22 datasets.
%
% 

%
%home_dir = '/nfs/olah/home/mjo/LidkeLab/Data';
home_dir = '/home/mjo/LidkeLab/Data/Emanuel';

pixel_size = 0.167; %microns
fid_pattern = 'fid*.mat'; %Patter of FID files

%% Emanuel 2015/06/19 SPTRegisterChannels
sptreg_dir = fullfile(home_dir,'Emanuel-SPT','150619-SPTRegisterChannels');
beads_file = fullfile(sptreg_dir,'Beads#0001-2015-6-19-15-6-6.mat');
bg_file = fullfile(sptreg_dir,'Background#0001-2015-6-19-14-52-47.mat');
sptreg_pattern = 'ChannelRegistation-2015-6-19-14-*.mat';
spt_reg = RegistrationAnalysis.loadSPTRegisterChannels(sptreg_dir, sptreg_pattern, beads_file, bg_file);

%% Sammy ChannelRegistrationV1  2013/2015
crv1_dir = fullfile(home_dir,'Emanuel-SPT','ChannelRegistrationV1');
crv1_pattern_2015 = 'ChannelRegistration-2015-02-26-12-*.mat';
crv1_2015_reg = RegistrationAnalysis.loadChannelRegistrationV1(crv1_dir, crv1_pattern_2015);
crv1_pattern_2013 = 'ChannelRegistration-2013-05-03-17-*.mat';
crv1_2013_reg = RegistrationAnalysis.loadChannelRegistrationV1(crv1_dir, crv1_pattern_2013);

%% Emanuel 2015/9/22 HMM Fiducials
hmm_calibration_dir = fullfile(home_dir,'Emanuel-SPT','Calibration');
fid_file = fullfile(hmm_calibration_dir,'fiducial-2015-9-22-17-33.mat'); %This one has cuttoff on top row.  Cannot be accurately fit
beads_file = fullfile(hmm_calibration_dir,'Beads-2015-9-22-15-15-60.mat');
bg_file = fullfile(hmm_calibration_dir,'Background-2015-9-22-15-3-44.mat');
regHMM = RegistrationAnalysis.loadHMMFiducialData(fid_file, beads_file, bg_file, pixel_size);

%% Emanuel 2015/10/28 Grid fiducuals
data_dir = fullfile(home_dir,'Emanuel-SPT/151028_testing-fiducials');
beads_file = fullfile(data_dir,'gain_fid1.mat');
bg_file = fullfile(data_dir,'background1.mat');
reg1028 = RegistrationAnalysis(data_dir,fid_pattern,beads_file, bg_file, pixel_size);

%% Emanuel 2015/10/29 Gain 0 Grid fiducuals
data_dir = fullfile(home_dir,'Emanuel-SPT/151029_testing-fiducials/EMgain0/With-zoom');
beads_file = fullfile(data_dir,'Beads-2015-10-29-12-40-10.mat');
bg_file = fullfile(data_dir,'Background-2015-10-29-12-27-55.mat');
reg1029_gain0 = RegistrationAnalysis(data_dir,fid_pattern,beads_file, bg_file, pixel_size);

%% Emanuel 2015/10/29 Gain 2500 Grid fiducuals
data_dir = fullfile(home_dir,'Emanuel-SPT/151029_testing-fiducials/EMgain2500');
beads_file = fullfile(data_dir,'Beads-2015-10-29-12-48-9.mat');
bg_file = fullfile(data_dir,'Background-2015-10-29-12-47-3.mat');
reg1029_gain2500 = RegistrationAnalysis(data_dir,fid_pattern,beads_file, bg_file, pixel_size);
