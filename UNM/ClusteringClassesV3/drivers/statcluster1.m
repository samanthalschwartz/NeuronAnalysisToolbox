clear all
close all

DataDirs = { ...
'.'
}

collapsing = true;

c = Clustering();
%c.Fig_ext = '';
H_nm = 1000;
V_nm = 1000;
particle_types = {'data'};

for j = 1:length(DataDirs)
   DataDir = DataDirs{j};
   Files = dir(fullfile(DataDir, '*_ROI.mat'));
   FileName = fullfile(DataDir, Files.name);
   BaseName = regexprep(Files.name, '_Results_ROI.mat', '');
   load(FileName);

   %if isnan(Sigma_Reg)
      Sigma_Reg = [0, 0];
   %end

   for i = 1:length(RoI)
      X = double(RoI{i}.X);   % nm
      Y = double(RoI{i}.Y);   % nm
      X_STD = double(RoI{i}.X_STD);   % nm
      Y_STD = double(RoI{i}.Y_STD);   % nm

      if collapsing
         SRc = SRcluster();
         SRc.ShrinkFactor = 0.5;
         [xy_SR, sigma_SR, combined] = ...
            SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
         XY = xy_SR;
      else
         XY = [X, Y];
      end

      c.Results = 'results';
      P{1} = XY;
      base_name = sprintf('%s-ROI%d', BaseName, i);
      c.cluster_stats('Hopkins', P, base_name, particle_types, H_nm, V_nm);
      c.cluster_stats('Ripley',  P, base_name, particle_types, H_nm, V_nm);
   end
end
