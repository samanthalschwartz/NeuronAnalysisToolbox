%
% Example of batch processing for Emanuel's 2015/9/22 dataset.
%
% 

%% Globals
%data_dir = '/home/mjo/LidkeLab/Data/Emanuel/150922/';
data_dir = '/nfs/olah/home/mjo/LidkeLab/Data/Emanuel/150922/';
dataset_name = '2015-9-22';
pixelSize = 0.167;
roiNames = {'Channel1', 'Channel2'};
wavelengths = [585, 655];

%% Fiducial data reading
calibration_dir = fullfile(data_dir,'Calibration');
all_fid_files = {'fiducial-2015-9-22-15-55.mat','fiducial-2015-9-22-17-33.mat','fiducial-2015-9-22-18-51.mat'};
beads_file = fullfile(calibration_dir,'Beads-2015-9-22-15-15-60.mat');
bg_file = fullfile(calibration_dir,'Background-2015-9-22-15-3-44.mat');
fid_file=fullfile(calibration_dir,all_fid_files{1});
fprintf('Using HMM Fiducial Registration File: %s\n',fid_file);

reg = RegistrationAnalysis.loadHMMFiducialData(fid_file, beads_file, bg_file, pixelSize);
[chMapFunc, chMapStruct] = reg.getOptimalMapping();
fprintf('Optimal Mapping Alg:%s RMSE:%.4g nm MaxError:%.4g nm\n',chMapStruct.algorithm, chMapStruct.RMSE, chMapStruct.maxError);
reg.save();
fprintf('Registration saved as: %s\n',reg.saveFilePath);
regFile = reg.saveFilePath;


%% Batch SPData processing
% This allows the easy creation of SPData files for each dataset in the folder
example_spdata = fullfile(data_dir,'HA585655_EGF-50nM-2015-9-22-17-3-58.spdata'); %location of an 'example' spdata to make all the others datasets like
spdata_overwrite = false; % Set as false does not overwite existing .spdata files
spdata_files = SPData.batchProcess(data_dir, '*.mat', example_spdata, spdata_overwrite);
fprintf('SPData Batch Processing Complete!\n');
fprintf('->Processed %i .spdata Files \n',numel(spdata_files));
disp(spdata_files');

%% Batch RPT Tracking
% This allows the easy tracking of dataset for each ROI in an spdata.  We provide an example for each
% of the two ROIs from the SPData which correspond the two different channels.  In general both channels 
% from a single dataset should be tracked individually and saved to produce the example files.
example_rpt_left = fullfile(data_dir,'RPT','HA585655_EGF-50nM-2015-9-22-17-3-58_Channel1.rpt'); %Example RPT for the left (channel1) ROI
example_rpt_right = fullfile(data_dir,'RPT','HA585655_EGF-50nM-2015-9-22-17-3-58_Channel2.rpt'); %Example RPT for the right (channel2) ROI
example_rpt = {example_rpt_left, example_rpt_right};
ROInames = {'Channel1','Channel2'};
rpt_overwrite=false; % Set as false does not overwite existing .rpt files
RPT.batchProcess(data_dir, '*.spdata', ROInames, example_rpt, rpt_overwrite);

%% Batch RPT Registration Analysis

RegisteredPairAnalysis.batchProcess(data_dir,'*.spdata', roiNames, regFile, wavelengths);

% 
% %% Get all pairs from every dataset and save them in one big cell array
% allPairs = RegisteredPairAnalysis.batchPairExport(data_dir,'*.regpairs');
% 
% %Check a subdirectory 'Analysis' exists
% analysisDir = fullfile(data_dir,'Analysis');
% if ~exist(analysisDir,'dir');
%     mkdir(analysisDir);
% end
% savePairsFile = fullfile(analysisDir,sprintf('pairs.%s.mat',dataset_name)); % File to save pairs datastructure in
% save(savePairsFile, 'allPairs','-mat'); %Save the pairs 
