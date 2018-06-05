%
% Example of batch processing for Emanuel's 2015/9/22 dataset.
%
% 

%% Globals
data_dir = '/home/mjo/LidkeLab/Data/Emanuel/50922/';
dataset_name = '2015-9-22';
pixelSize = 0.167;

%% Fiducial data reading
calibration_dir = fullfile(data_dir,'Calibration');
fid_files = cellmap(@(f) fullfile(calibration_dir,f), ...
    {'fiducial-2015-9-22-15-55.mat','fiducial-2015-9-22-17-33.mat','fiducial-2015-9-22-18-51.mat'});
beads_file = fullfile(calibration_dir,'Beads-2015-9-22-15-15-60.mat');
bg_file = fullfile(calibration_dir,'Background-2015-9-22-15-3-44.mat');

reg = RegistrationAnalysis(calibration_dir, fid_files, beads_file, bg_file, pixelSize);


%This method reads the data captured with the old HMM-related Fiducial Registration Capture software
% we convert this into our more universal fiducialData structure which simply allows it to be used as
% a registration method to create a new RegisterdPairAnalysis object;
fiducialData = RegisteredPairAnalysis.readHMMFiducialRegistrationData(fid_file, beads_file, bg_file);
RegisteredPairAnalysis.testFiducialRegistration(fiducialData);

%% Batch SPData processing
% This allows the easy creation of SPData files for each dataset in the folder

example_spdata = fullfile(data_dir,'.spdata'); %location of an 'example' spdata to make all the others datasets like
spdata_overwrite = false; % Set as false does not overwite existing .spdata files
SPData.batchProcess(data_dir, '*.mat', example_spdata, spdata_overwrite);

%% Batch RPT Tracking
% This allows the easy tracking of dataset for each ROI in an spdata.  We provide an example for each
% of the two ROIs from the SPData which correspond the two different channels.  In general both channels 
% from a single dataset should be tracked individually and saved to produce the example files.
example_rpt_left = fullfile(data_dir,'RPT','.rpt'); %Example RPT for the left (channel1) ROI
example_rpt_right = fullfile(data_dir,'RPT','.rpt'); %Example RPT for the right (channel2) ROI
example_rpt = {example_rpt_left, example_rpt_right};
rpt_overwrite=false; % Set as false does not overwite existing .rpt files
RPT.batchProcess(data_dir, '*.spdata', example_rpt, rpt_overwrite);

%% Batch RPT Registration Analysis
example_regpairs = fullfile(data_dir,'RegisteredPairAnalysis','.regpairs');
RegisteredPairAnalysis.batchProcess(data_dir,'*.spdata',example_regpairs);

%% Get all pairs from every dataset and save them in one big cell array
allPairs = RegisteredPairAnalysis.batchPairExport(data_dir,'*.regpairs');

%Check a subdirectory 'Analysis' exists
analysisDir = fullfile(data_dir,'Analysis');
if ~exist(analysisDir,'dir');
    mkdir(analysisDir);
end
savePairsFile = fullfile(analysisDir,sprintf('pairs.%s.mat',dataset_name)); % File to save pairs datastructure in
save(savePairsFile, 'allPairs','-mat'); %Save the pairs 
