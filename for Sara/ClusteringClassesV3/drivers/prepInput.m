clear all
close all

%Here's the path for all the data
%O:\Cell Path\Lidke Lab\Arianne\SR microscopy\For Sam\no IL-3
%O:\Cell Path\Lidke Lab\Arianne\SR microscopy\Data for pair corr\Control
%clustdatadir = 'Z:\Arianne\SR microscopy\For Sam\no IL-3'
%srdatadir    = 'Z:\Arianne\SR microscopy\Data for pair corr\Control';
clustdatadir = '.';
srdatadir    = '.';

pixel2nm = 16000/150;

clusfiles = dir(fullfile(clustdatadir,'*.mat'));

SRc = SRcluster();   % constructor
SRc.PlotFigures = true;

SRc.E = 50;       % cutoff distance (nm)
SRc.LoS = 0.01;   % level of significance
SRc.A_ROI = 1;    % area of the ROI (nm^2)

for ii = 1%:numel(clusfiles) %loop through each file
   %load(fullfile(clustdatadir,clusfiles(ii).name));
   %stsplit_output = strsplit(clusfiles(ii).name,'_Cell')
   %srdemoname = stsplit_output{1};
   load('../data/Control_0017_2014_4_22_13_6_49_Cell4variables.mat');
   for rois = 1%:numel(regs) %loop through each roi within the file
      % make input to function call:
      % get coordinate matrix
      coordmatrix = regs(rois).points;   % nm
      % get corresponding sigma values
      %load(fullfile(srdatadir,[srdemoname '.mat']));
      load('../data/Control_0017_2014_4_22_13_6_49.mat');
      ids = arrayfun(@(x, y) find(SRtest.Results.X*106 == x &        ...
                                  SRtest.Results.Y*106 == y),        ...
                     regs(rois).points(:,1), regs(rois).points(:,2), ...
                     'UniformOutput', false);
      stdevmatrix = [double(SRtest.Results.X_STD(cell2mat(ids))), ...
                     double(SRtest.Results.Y_STD(cell2mat(ids)))] .* pixel2nm;
    
      % get the corresponding stdev of Drift Correction
      Sigma_Reg = std(SRtest.DriftCorrect_XYShift) .* pixel2nm; 
    
      %------------now call function and save output somehow:
      % format of output needs to be modified   -------
      SRc.clusterSR(coordmatrix, stdevmatrix, Sigma_Reg);
      SRcollapseFig = SRc.plotSRcollapse();
      %cutoffFigs = SRc.cutoffPlots();
      %SRclusterFig = SRc.plotSRclusters();
      %[results, analysisFigs] = SRc.analyzeSRclusters();
   end
end
