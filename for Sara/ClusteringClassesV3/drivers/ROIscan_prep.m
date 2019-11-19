clear all
close all

DataDirs = {...
'.', ...
%'D:\hD-FC SR data\20160916 hD-FC sext\Yst3\Results',...
%'D:\hD-FC SR data\20160916 hD-FC sext\Yst6\Results',...
%'D:\hD-FC SR data\20160923 hD-FC T35_sext\sext\Yst1\Results',...
%'D:\hD-FC SR data\20160923 hD-FC T35_sext\sext\Yst2\Results',...
}

SRc = SRcluster();
%SRc.PvalueStatistics = true;
SRc.ShrinkFactor = 0.5;
SRc.LoS = 0.01;
SRc.A_ROI = (1000*1000);
SRc.PlotFigures = false;
SRc.Printing = false;

n_ROIs = 0;
for j = 1 : numel(DataDirs)
   DataDir = DataDirs{j};
   Files = dir(fullfile(DataDir, '*_ROI.mat'));
   FileName = fullfile(DataDir, Files.name)
   load(FileName);
    
   %if isnan(Sigma_Reg)
      Sigma_Reg = [0, 0];
   %end

   for i = 1 : numel(RoI)
      n_ROIs = n_ROIs + 1;
      % Assume single label.
      X = double(RoI{i}.X{1});   % nm
      Y = double(RoI{i}.Y{1});   % nm
      X_STD = double(RoI{i}.X_STD{1});   % nm
      Y_STD = double(RoI{i}.Y_STD{1});   % nm
  
      [xy_SR, sigma_SR, combined] = ...
         SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
      XY_orig{n_ROIs} = [X, Y];
      Sigma_orig{n_ROIs} = [X_STD, Y_STD];
      XY_SR{n_ROIs} = xy_SR;
      Sigma_SR{n_ROIs} = sigma_SR;
      Combined{n_ROIs} = combined;
   end
end

save('ROI_collapsed.mat', ...
     'XY_orig', 'Sigma_orig', 'XY_SR', 'Sigma_SR', 'Combined');
