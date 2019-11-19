clear all
close all

DataDirs = {...
'F:\Reanalyzed data\MG20150521-YstAB\Calb\Yst6\Results3',...
'F:\Reanalyzed data\MG20150521-YstAB\Calb\Yst8\Results3',...
'F:\Reanalyzed data\MG20150521-YstAB\Calb\Yst9\Results3',...
'F:\Reanalyzed data\MG20150521-YstAB\Calb\Yst10\Results3',...
'F:\Reanalyzed data\MG20150521-YstAB\Calb\Yst11\Results3',...
    }

n_objs_collapsed = [];
singles_per_total_clusters = [];
min_c2c_dist_total = [];
radius_per_2orMore_cluster = [];
numclust_123 = [];  numclust_all = [];
numdensity_3 = [];

for j = 1:length(DataDirs)
    DataDir = DataDirs{j};
    Files = dir([DataDir, '\*_ROI.mat']);
    FileName = [DataDir, '\', Files.name]
    load(FileName);
    
      if isnan(Sigma_Reg)
         Sigma_Reg = [10, 10];
      end

    for i = 1:length(RoI)
        X = double(RoI{i}.X);   % nm
        Y = double(RoI{i}.Y);   % nm
        X_STD = double(RoI{i}.X_STD);   % nm
        Y_STD = double(RoI{i}.Y_STD);   % nm

        SRc = SRcluster();
        SRc.Cutoff = 10 : 10 : 1000
        SRc.LoS
        SRc.A_ROI = (1500*1500);
        SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
        %cutoffFigs = SRc.cutoffPlots();
        SRc.PlotFigures=false;
        SRc.Algorithm = 'DBSCAN_Daszykowski';   SRc.minPts = 3;   SRc.E = 50;

         [results, analysisFigs] = SRc.analyzeSRclusters();
         n_objs_collapsed = [n_objs_collapsed, results.n_objs_collapsed_pass1];
         singles_per_total_clusters = [singles_per_total_clusters, results.singles_per_total_clusters];
         min_c2c_dist_total = [min_c2c_dist_total, results.min_c2c_dists];
         radius_per_2orMore_cluster = [radius_per_2orMore_cluster, results.radius_per_2orMore_cluster];
         numclust_123 = [numclust_123; results.numclust]; 
         numclust_all = [numclust_all, sum(results.numclust)];
         numdensity_3 = [numdensity_3, results.numdensity(3)];
    end
end

Dir = 'C:\Users\MGraus\Documents\MATLAB\SRcluster\';
save([Dir, 'AYst.mat'], 'n_objs_collapsed', 'singles_per_total_clusters', 'min_c2c_dist_total',...
   'radius_per_2orMore_cluster', 'numclust_123', 'numclust_all', 'numdensity_3');

