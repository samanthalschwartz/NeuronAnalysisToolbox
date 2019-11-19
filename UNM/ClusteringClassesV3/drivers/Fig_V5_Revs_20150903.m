% SR yeast beta-glucan exposure
%revision figures 5 and Ripley's for 3rd submission , Fig_V5_20150520
% started by Jia.L, 2015.08.28
% using DBSCAN for post-collapsed clustering, 2015.09.07
% add simulation data, 2015.09.24

%% Simulation 1x10^(-5) domain density
close all
clear all
clc

cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Fig_Simulation')

Rho_d=1.0e-05;%domain density(nm^-2)
nSim=100;% how many time to run the simulation
for ii=1:nSim  
    [results, results_o, results_c, results_t,results_simulated] = testSingleLabel(Rho_d);
    N_Collapsed{ii}=results.numobjs_collapsed_pass1;
    N_Col2True{ii}=results_simulated.N_observations;%observations per true localization
    N_near2truedist{ii}=results_simulated.nearest_true_dists;
    N_true{ii}=results_t.n_points;
    N_collapsed{ii}=results_c.n_points;
    N_observed{ii}=results_o.n_points;
end

clear results results_c results_o results_simulated results_t ii
close all
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Fig_Simulation')
save 10e-5

%over lap the 2 figure
figure(1)
ah1=axes;
hist(cell2mat(N_Collapsed),40);% cyan dots vs. blue dots
hold on
hist(cell2mat(N_Col2True),20);% cyan dots vs. blue dots
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');%cyan vs blue, red bar
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k');%cyan vs red, blue bar
hold off
box off
title({'Simulations'})
ylabel({'Frequency'})
%xlabel({'Collapsed observations per localization'})
legend('localizations','simulated localizations')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(ah1,'LineWidth',3,'FontSize',18)

figure(2)
ah2=axes;
figurematrix=[cell2mat(N_true),...
             cell2mat(N_collapsed),...
             cell2mat(N_observed)];

figuregroup=[(ones(size(cell2mat(N_true)))),...
             (ones(size(cell2mat(N_collapsed))))+1,...
             (ones(size(cell2mat(N_observed))))+2]; 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Simulations'})
ylabel({'Counts'})
xlabel({''})
set(gcf,'color',[1 1 1])
set(ah2,'XTickLabel',{'Simulated localization ','Localizations','Observations'},...
     'XTick',[1 2 3])
set(ah2,'LineWidth',3,'FontSize',18)

figure(3)
ah3=axes;
hist(cell2mat(N_near2truedist),40);% cyan dots vs. blue dots
box off
title({'Simulations'})
ylabel({'Frequency'})
xlabel({'nearest to true localization'})
set(gcf,'color',[1 1 1])
set(ah3,'LineWidth',3,'FontSize',18)

saveas(figure(1),'Hist_Simulation_Distribution.fig')
saveas(figure(1),'Hist_Simulation_Distribution.tif')

saveas(figure(2),'Sim_STATS_Counting.fig')
saveas(figure(2),'Sim_STATS_Counting.tif')

saveas(figure(3),'Hist_Dist.fig')
saveas(figure(3),'Hist_Dist.tif')

%% Simulation 5x10^(-5) domain density
close all
clear all
clc

cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Fig_Simulation')

Rho_d=5.0e-05;%domain density(nm^-2)
nSim=100;% how many time to run the simulation
for ii=1:nSim  
    [results, results_o, results_c, results_t,results_simulated] = testSingleLabel(Rho_d);
    N_Collapsed{ii}=results.numobjs_collapsed_pass1;
    N_Col2True{ii}=results_simulated.N_observations;%observations per true localization
    N_near2truedist{ii}=results_simulated.nearest_true_dists;
    N_true{ii}=results_t.n_points;
    N_collapsed{ii}=results_c.n_points;
    N_observed{ii}=results_o.n_points;
end

clear results results_c results_o results_simulated results_t ii
close all
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Fig_Simulation')
save 50e-5

%over lap the 2 figure
figure(1)
ah1=axes;
hist(cell2mat(N_Collapsed),40);% cyan dots vs. blue dots
hold on
hist(cell2mat(N_Col2True),20);% cyan dots vs. blue dots
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');%cyan vs blue, red bar
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k');%cyan vs red, blue bar
hold off
box off
title({'Simulations'})
ylabel({'Frequency'})
%xlabel({'Collapsed observations per localization'})
legend('localizations','simulated localizations')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(ah1,'LineWidth',3,'FontSize',18)

figure(2)
ah2=axes;
figurematrix=[cell2mat(N_true),...
             cell2mat(N_collapsed),...
             cell2mat(N_observed)];

figuregroup=[(ones(size(cell2mat(N_true)))),...
             (ones(size(cell2mat(N_collapsed))))+1,...
             (ones(size(cell2mat(N_observed))))+2]; 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Simulations'})
ylabel({'Counts'})
xlabel({''})
set(gcf,'color',[1 1 1])
set(ah2,'XTickLabel',{'Simulated localization ','Localizations','Observations'},...
     'XTick',[1 2 3])
set(ah2,'LineWidth',3,'FontSize',18)

figure(3)
ah3=axes;
hist(cell2mat(N_near2truedist),40);% cyan dots vs. blue dots
box off
title({'Simulations'})
ylabel({'Frequency'})
xlabel({'nearest to true localization'})
set(gcf,'color',[1 1 1])
set(ah3,'LineWidth',3,'FontSize',18)

saveas(figure(1),'Hist_Simulation_Distribution_5.fig')
saveas(figure(1),'Hist_Simulation_Distribution_5.tif')

saveas(figure(2),'Sim_STATS_Counting_5.fig')
saveas(figure(2),'Sim_STATS_Counting_5.tif')

saveas(figure(3),'Hist_Dist_5.fig')
saveas(figure(3),'Hist_Dist_5.tif')

%%  Yst Ctrl save ROIs
clear 
close all
clc

pixel2nm = 16000/150;%convert pixel size into nm
x_size=1000;
y_size=1000;
DataDirs={...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl1/FrameConectResults',...% Jia's
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl2/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl2/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl2/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl2/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl2/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl3/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl/YstCtrl4/FrameConectResults',...
};

% save ROI for each cell
for ii=1:length(DataDirs)
    DataDir=DataDirs{ii};
    cd(DataDir)
    load('FrameConect.mat')
    
    X = double(SRD.Results.X) .* pixel2nm;   % nm
    Y = double(SRD.Results.Y) .* pixel2nm;   % nm
    
    X_STD = double(SRD.Results.X_STD) .* pixel2nm;   % nm
    Y_STD = double(SRD.Results.Y_STD) .* pixel2nm;   % nm

    [index_ROI,ROI,n_ROIs] = get_ROI(X, Y, x_size, y_size);% hit Q on keyboard to quit
    
    XROIs{ii}=X(index_ROI{n_ROIs});
    YROIs{ii}=Y(index_ROI{n_ROIs});
    
    XROIs_STD{ii}=X_STD(index_ROI{n_ROIs});
    YROIs_STD{ii}=Y_STD(index_ROI{n_ROIs});
    
    close all
end
    
YstCtrlROIs.Xs=XROIs;
YstCtrlROIs.Ys=YROIs;
YstCtrlROIs.X_STDs=XROIs_STD;
YstCtrlROIs.Y_STDs=XROIs_STD;

cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl')
save('YstCtrlROIs','YstCtrlROIs')

%% Yst Casp save ROIs
clear 
close all
clc

pixel2nm = 16000/150;%convert pixel size into nm
x_size=1000;
y_size=1000;
DataDirs={...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp1/FrameConectResults',...% Jia's
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp2/FrameConectResults',... % has 2 ysts
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp2/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp3/FrameConectResults',...% has 2 ysts
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp3/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp4/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp5/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp6/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp/YstCasp7/FrameConectResults',...
};
% save ROI for each cell
for ii=1:length(DataDirs)
    DataDir=DataDirs{ii};
    cd(DataDir)
    load('FrameConect.mat')
    
    X = double(SRD.Results.X) .* pixel2nm;   % nm
    Y = double(SRD.Results.Y) .* pixel2nm;   % nm
    
    X_STD = double(SRD.Results.X_STD) .* pixel2nm;   % nm
    Y_STD = double(SRD.Results.Y_STD) .* pixel2nm;   % nm

    [index_ROI,ROI,n_ROIs] = get_ROI(X, Y, x_size, y_size);% hit Q on keyboard to quit
    
    XROIs{ii}=X(index_ROI{n_ROIs});
    YROIs{ii}=Y(index_ROI{n_ROIs});
    
    XROIs_STD{ii}=X_STD(index_ROI{n_ROIs});
    YROIs_STD{ii}=Y_STD(index_ROI{n_ROIs});
    
    close all
end
    
YstCaspROIs.Xs=XROIs;
YstCaspROIs.Ys=YROIs;
YstCaspROIs.X_STDs=XROIs_STD;
YstCaspROIs.Y_STDs=XROIs_STD;

cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp')
save('YstCaspROIs','YstCaspROIs')

%% Hyp Ctrl save ROIs
clear 
close all
clc

pixel2nm = 16000/150;%convert pixel size into nm
x_size=1000;
y_size=1000;
DataDirs={...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl1/FrameConectResults',...% Jia's
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl2/FrameConectResults',... % 
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl3/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl4/FrameConectResults',...% 
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl5/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl6/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl7/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl8/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl/HypCtrl8/FrameConectResults',... % has 2
};
% save ROI for each cell
for ii=1:length(DataDirs)
    DataDir=DataDirs{ii};
    cd(DataDir)
    load('FrameConect.mat')
    
    X = double(SRD.Results.X) .* pixel2nm;   % nm
    Y = double(SRD.Results.Y) .* pixel2nm;   % nm
    
    X_STD = double(SRD.Results.X_STD) .* pixel2nm;   % nm
    Y_STD = double(SRD.Results.Y_STD) .* pixel2nm;   % nm

    [index_ROI,ROI,n_ROIs] = get_ROI(X, Y, x_size, y_size);% hit Q on keyboard to quit
    
    XROIs{ii}=X(index_ROI{n_ROIs});
    YROIs{ii}=Y(index_ROI{n_ROIs});
    
    XROIs_STD{ii}=X_STD(index_ROI{n_ROIs});
    YROIs_STD{ii}=Y_STD(index_ROI{n_ROIs});
    
    close all
end
    
HypCtrlROIs.Xs=XROIs;
HypCtrlROIs.Ys=YROIs;
HypCtrlROIs.X_STDs=XROIs_STD;
HypCtrlROIs.Y_STDs=XROIs_STD;

cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl')
save('HypCtrlROIs','HypCtrlROIs')

%% Hyp Casp save ROIs
clear 
close all
clc

pixel2nm = 16000/150;%convert pixel size into nm
x_size=1000;
y_size=1000;
DataDirs={...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp1/FrameConectResults',...% Jia's
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp4/FrameConectResults',...% 
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp4/FrameConectResults',...% 4
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp4/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp4/FrameConectResults',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp5/FrameConectResults',...
};

%'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp2/FrameConectResults',... % 
%'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp/HypCasp3/FrameConectResults',...
% save ROI for each cell
for ii=1:length(DataDirs)
    DataDir=DataDirs{ii};
    cd(DataDir)
    load('FrameConect.mat')
    
    X = double(SRD.Results.X) .* pixel2nm;   % nm
    Y = double(SRD.Results.Y) .* pixel2nm;   % nm
    
    X_STD = double(SRD.Results.X_STD) .* pixel2nm;   % nm
    Y_STD = double(SRD.Results.Y_STD) .* pixel2nm;   % nm

    [index_ROI,ROI,n_ROIs] = get_ROI(X, Y, x_size, y_size);% hit Q on keyboard to quit
    
    XROIs{ii}=X(index_ROI{n_ROIs});
    YROIs{ii}=Y(index_ROI{n_ROIs});
    
    XROIs_STD{ii}=X_STD(index_ROI{n_ROIs});
    YROIs_STD{ii}=Y_STD(index_ROI{n_ROIs});
    
    close all
end
    
HypCaspROIs.Xs=XROIs;
HypCaspROIs.Ys=YROIs;
HypCaspROIs.X_STDs=XROIs_STD;
HypCaspROIs.Y_STDs=XROIs_STD;

cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp')
save('HypCaspROIs','HypCaspROIs')

%% Yst Ctrl Spatial Stats
clear all
close all
clc

Sigma_Reg = [0 0]; % no driftcorrection
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl')
load('YstCtrlROIs.mat')

nROIs=size(YstCtrlROIs.Xs,2);

for ii=1:nROIs
    close all
SRc = SRcluster;
SRc.A_ROI = 1.0e+06; % nm^2
SRc.clusterSR([YstCtrlROIs.Xs{ii}, YstCtrlROIs.Ys{ii}], [YstCtrlROIs.X_STDs{ii}, YstCtrlROIs.Y_STDs{ii}], Sigma_Reg);
%SRc.E = SRc.chooseCutoff();
SRc.E =50;
SRc.Algorithm = 'DBSCAN_Daszykowski';   
SRc.minPts = 3;
SRclusterFig = SRc.plotSRclusters();
SRc.PlotFigures=false;
[results, analysisFigs] = SRc.analyzeSRclusters();

YstCtrl.XYcollaped{ii}=SRc.XY;
YstCtrl.Clust3plusNum{ii} = SRc.cutoffPlots();

YstCtrl.NumColapsed{ii}=results.numobjs_collapsed_pass1;% number of object collapsed
YstCtrl.singlesFraction{ii}=results.singles_per_total_clusters;% fraction of singles/total number of clusters

YstCtrl.SingleDensity{ii}=results.numdensity(1);% number of cluster(single,double,multiple)/nm^2
YstCtrl.MultipleDensity{ii}=results.numdensity(3);

YstCtrl.MeanRadius2plus{ii}=results.radius_per_2orMore_cluster;% radii of 2 or 3+

YstCtrl.MeanArea3plus{ii}=results.area_per_multiple_cluster;% area of 3+

YstCtrl.MeanDist{ii}=results.min_c2c_dists;% mean nearer neighbor distance of everything
YstCtrl.ObjClust{ii}=results.numobjs_per_multiple_cluster;% number of objects per cluster
YstCtrl.DensObjClust{ii}=results.numobjs_per_area;% 'mean(# objects / area per multiple cluster) (nm^-2)

YstCtrl.NumObjs{ii}=results.numobjs;
YstCtrl.Clusts{ii}=results.numclust;

end
clear analysisFigs ii  SRc SRclusterFig YstCtrlROIs nROIs Sigma_Reg ans results
close all
save YstCtrlSTATSDBSCAN

%% Yst Casp Spatial Stats
clear all
close all
clc

Sigma_Reg = [0 0]; % no driftcorrection
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp')
load('YstCaspROIs.mat')

nROIs=size(YstCaspROIs.Xs,2);

for ii=1:nROIs
    close all
SRc = SRcluster;
SRc.A_ROI = 1.0e+06; % nm^2
SRc.clusterSR([YstCaspROIs.Xs{ii}, YstCaspROIs.Ys{ii}], [YstCaspROIs.X_STDs{ii}, YstCaspROIs.Y_STDs{ii}], Sigma_Reg);
SRc.E = 50;
SRc.Algorithm = 'DBSCAN_Daszykowski';   
SRc.minPts = 3;
SRclusterFig = SRc.plotSRclusters();
SRc.PlotFigures=false;
%SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();

YstCasp.XYcollaped{ii}=SRc.XY;
YstCasp.Clust3plusNum{ii} = SRc.cutoffPlots();
%YstCasp.Es{ii}=SRc.E; % methods used, parameters

YstCasp.NumColapsed{ii}=results.numobjs_collapsed_pass1;% number of object collapsed
YstCasp.singlesFraction{ii}=results.singles_per_total_clusters;% fraction of singles/total number of clusters

YstCasp.SingleDensity{ii}=results.numdensity(1);% number of cluster(single,double,multiple)/nm^2
YstCasp.MultipleDensity{ii}=results.numdensity(3);

YstCasp.MeanRadius2plus{ii}=results.radius_per_2orMore_cluster;% radii of 2 or 3+

YstCasp.MeanArea3plus{ii}=results.area_per_multiple_cluster;% area of 3+

YstCasp.MeanDist{ii}=results.min_c2c_dists;% mean nearest neighbor distance of everything
YstCasp.ObjClust{ii}=results.numobjs_per_multiple_cluster;% number of objects per cluster
YstCasp.DensObjClust{ii}=results.numobjs_per_area;% 'mean(# objects / area per multiple cluster) (nm^-2)

YstCasp.NumObjs{ii}=results.numobjs;
YstCasp.Clusts{ii}=results.numclust;
end

clear analysisFigs ii results SRc SRclusterFig YstCaspROIs nROIs Sigma_Reg ans
close all
save YstCaspSTATSDBSCAN

%% Hyp Ctrl Spatial Stats
clear all
close all
clc

Sigma_Reg = [0 0]; % no driftcorrection
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl')
load('HypCtrlROIs.mat')

nROIs=size(HypCtrlROIs.Xs,2);

for ii=1:nROIs
    close all
SRc = SRcluster;
SRc.A_ROI = 1.0e+06; % nm^2
SRc.clusterSR([HypCtrlROIs.Xs{ii}, HypCtrlROIs.Ys{ii}], [HypCtrlROIs.X_STDs{ii}, HypCtrlROIs.Y_STDs{ii}], Sigma_Reg);
SRc.E = 50;
SRc.Algorithm = 'DBSCAN_Daszykowski';   
SRc.minPts = 3;
SRclusterFig = SRc.plotSRclusters();
SRc.PlotFigures=false;%SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();

HypCtrl.XYcollaped{ii}=SRc.XY;
HypCtrl.Clust3plusNum{ii} = SRc.cutoffPlots();
%HypCtrl.Es{ii}=SRc.E; % methods used, parameters

HypCtrl.NumColapsed{ii}=results.numobjs_collapsed_pass1;% number of object collapsed
HypCtrl.singlesFraction{ii}=results.singles_per_total_clusters;% fraction of singles/total number of clusters

HypCtrl.SingleDensity{ii}=results.numdensity(1);% number of cluster(single,double,multiple)/nm^2
HypCtrl.MultipleDensity{ii}=results.numdensity(3);

HypCtrl.MeanRadius2plus{ii}=results.radius_per_2orMore_cluster;% radii of 2 or 3+

HypCtrl.MeanArea3plus{ii}=results.area_per_multiple_cluster;% area of 3+

HypCtrl.MeanDist{ii}=results.min_c2c_dists;% mean nearest neighbor distance of everything
HypCtrl.ObjClust{ii}=results.numobjs_per_multiple_cluster;% number of objects per cluster
HypCtrl.DensObjClust{ii}=results.numobjs_per_area;% 'mean(# objects / area per multiple cluster) (nm^-2)

HypCtrl.NumObjs{ii}=results.numobjs;
HypCtrl.Clusts{ii}=results.numclust;
end

clear analysisFigs ii results SRc SRclusterFig HypCtrlROIs nROIs Sigma_Reg ans
close all
save HypCtrlSTATSDBSCAN

%% Hyp Casp Spatial Stats
clear all
close all
clc

Sigma_Reg = [0 0]; % no driftcorrection
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp')
load('HypCaspROIs.mat')

nROIs=size(HypCaspROIs.Xs,2);

for ii=1:nROIs
    close all
SRc = SRcluster;
SRc.A_ROI = 1.0e+06; % nm^2
SRc.clusterSR([HypCaspROIs.Xs{ii}, HypCaspROIs.Ys{ii}], [HypCaspROIs.X_STDs{ii}, HypCaspROIs.Y_STDs{ii}], Sigma_Reg);
SRc.E =50;
SRc.Algorithm = 'DBSCAN_Daszykowski';   
SRc.minPts = 3;
SRclusterFig = SRc.plotSRclusters();
SRc.PlotFigures=false;
%SRclusterFig = SRc.plotSRclusters();
[results, analysisFigs] = SRc.analyzeSRclusters();

HypCasp.XYcollaped{ii}=SRc.XY;
HypCasp.Clust3plusNum{ii} = SRc.cutoffPlots();
%HypCasp.Es{ii}=SRc.E; % methods used, parameters

HypCasp.NumColapsed{ii}=results.numobjs_collapsed_pass1;% number of object collapsed
HypCasp.singlesFraction{ii}=results.singles_per_total_clusters;% fraction of singles/total number of clusters

HypCasp.SingleDensity{ii}=results.numdensity(1);% number of cluster(single,double,multiple)/nm^2
HypCasp.MultipleDensity{ii}=results.numdensity(3);

HypCasp.MeanRadius2plus{ii}=results.radius_per_2orMore_cluster;% radii of 2 or 3+

HypCasp.MeanArea3plus{ii}=results.area_per_multiple_cluster;% area of 3+

HypCasp.MeanDist{ii}=results.min_c2c_dists;% mean nearest neighbor distance of everything
HypCasp.ObjClust{ii}=results.numobjs_per_multiple_cluster;% number of objects per cluster
HypCasp.DensObjClust{ii}=results.numobjs_per_area;% 'mean(# objects / area per multiple cluster) (nm^-2)

HypCasp.NumObjs{ii}=results.numobjs;
HypCasp.Clusts{ii}=results.numclust;
end

clear analysisFigs ii results SRc SRclusterFig HypCaspROIs nROIs Sigma_Reg ans
close all
save HypCaspSTATSDBSCAN

%% Combine all results and plots
close all
clear
clc

% Yst Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl')
load('YstCtrlSTATSDBSCAN.mat')

% Yst Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp')
load('YstCaspSTATSDBSCAN.mat')

% Hyp Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl')
load('HypCtrlSTATSDBSCAN.mat')

% Hyp Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp')
load('HypCaspSTATSDBSCAN.mat')

% number collapsed
figure(1)
aNumColapsd=axes;
figurematrix=[cell2mat(YstCtrl.NumColapsed), cell2mat(YstCasp.NumColapsed), cell2mat(HypCtrl.NumColapsed), cell2mat(HypCasp.NumColapsed)];
figuregroup=[(ones(size(cell2mat(YstCtrl.NumColapsed)))) (ones(size(cell2mat(YstCasp.NumColapsed)))+1)...
     (ones(size(cell2mat(HypCtrl.NumColapsed)))+3) (ones(size(cell2mat(HypCasp.NumColapsed)))+4)]; 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Number of Observations Collapsed'})
ylabel({'Collapsed observations per localization'})
set(aNumColapsd,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aNumColapsd,'LineWidth',3,'FontSize',18)

% Mean of Radii for multiples
figure(2)
aRad2plus=axes;
figurematrix=[cell2mat(YstCtrl.MeanRadius2plus), cell2mat(YstCasp.MeanRadius2plus),...
    cell2mat(HypCtrl.MeanRadius2plus), cell2mat(HypCasp.MeanRadius2plus)];
figuregroup=[(ones(size(cell2mat(YstCtrl.MeanRadius2plus)))) (ones(size(cell2mat(YstCasp.MeanRadius2plus)))+1)...
     (ones(size(cell2mat(HypCtrl.MeanRadius2plus)))+3) (ones(size(cell2mat(HypCasp.MeanRadius2plus)))+4)]; 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Mean Cluster Radii'})
ylabel({'Equivalent Multi-exposure Radius (nm)'})
set(aRad2plus,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aRad2plus,'LineWidth',3,'FontSize',18)

% Mean of nearest neighbor distance . Singlets+Multiple Distance
figure(3)
aClustDist=axes;

figurematrix=[cell2mat(YstCtrl.MeanDist),...
    cell2mat([YstCasp.MeanDist]),...
    cell2mat(HypCtrl.MeanDist),...
    cell2mat(HypCasp.MeanDist)];

figuregroup=[(ones(size(cell2mat(YstCtrl.MeanDist)))),...
    (ones(size(cell2mat(YstCasp.MeanDist)))+1),...
     (ones(size(cell2mat(HypCtrl.MeanDist)))+3),...
     (ones(size(cell2mat(HypCasp.MeanDist)))+4)]; 
 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Mean of Nearest Neighbor Distance (nm)'})
ylabel({'Mean of Nearest Neighbor Distance (nm)'})
set(aClustDist,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aClustDist,'LineWidth',3,'FontSize',18)

%percentage of Singlets/total cluster
figure(4)
aSingleFract=axes;
figurematrix=[cell2mat(YstCtrl.singlesFraction), cell2mat(YstCasp.singlesFraction),...
    cell2mat(HypCtrl.singlesFraction), cell2mat(HypCasp.singlesFraction)];
figuregroup=[(ones(size(cell2mat(YstCtrl.singlesFraction)))) (ones(size(cell2mat(YstCasp.singlesFraction)))+1)...
     (ones(size(cell2mat(HypCtrl.singlesFraction)))+3) (ones(size(cell2mat(HypCasp.singlesFraction)))+4)]; 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Fraction of Singlet Glucan Exposures'})
ylabel({'Fraction of Singlet Glucan Exposures'})
set(aSingleFract,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aSingleFract,'LineWidth',3,'FontSize',18)

%Density of localizations per ROIs
figure(5)
aObjDens=axes;
figurematrix=[cellfun(@sum,YstCtrl.NumObjs),...
    cellfun(@sum,YstCasp.NumObjs),...
    cellfun(@sum,HypCtrl.NumObjs),...
    cellfun(@sum,HypCasp.NumObjs)];

figuregroup=[(ones(size(cellfun(@sum,YstCtrl.NumObjs)))),...
    (ones(size(cellfun(@sum,YstCasp.NumObjs)))+1)...
     (ones(size(cellfun(@sum,HypCtrl.NumObjs)))+3),...
     (ones(size(cellfun(@sum,HypCasp.NumObjs)))+4)];
 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Localization Density (nm^-^2)'})
ylabel({'Localization Density (nm^-^2)'})
set(aObjDens,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aObjDens,'LineWidth',3,'FontSize',18)


%Density of glucan exposure per ROIs
figure(6)
aExpDens=axes;

figurematrix=[mean(cell2mat([YstCtrl.SingleDensity',YstCtrl.MultipleDensity']),2)',...
    mean(cell2mat([YstCasp.SingleDensity',YstCasp.MultipleDensity']),2)',...
    mean(cell2mat([HypCtrl.SingleDensity',HypCtrl.MultipleDensity']),2)',...
    mean(cell2mat([HypCasp.SingleDensity',HypCasp.MultipleDensity']),2)'];

figuregroup=[(ones(size(mean(cell2mat([YstCtrl.SingleDensity',YstCtrl.MultipleDensity']),2)'))),...
    (ones(size(mean(cell2mat([YstCasp.SingleDensity',YstCasp.MultipleDensity']),2)'))+1),...
     (ones(size(mean(cell2mat([HypCtrl.SingleDensity',HypCtrl.MultipleDensity']),2)'))+3),...
     (ones(size( mean(cell2mat([HypCasp.SingleDensity',HypCasp.MultipleDensity']),2)'))+4)]; 
 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Density of Glucan Exposure'})
ylabel({'Glucan Exposure Density (nm^-^2)'})
set(aExpDens,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aExpDens,'LineWidth',3,'FontSize',18)


%number of objects per cluster. 
figure(7)
aObjClust=axes;

figurematrix=[cell2mat(YstCtrl.ObjClust),...
             cell2mat(YstCasp.ObjClust),...
             cell2mat(HypCtrl.ObjClust),...
             cell2mat(HypCasp.ObjClust)];

figuregroup=[(ones(size(cell2mat(YstCtrl.ObjClust)))),...
             (ones(size(cell2mat(YstCasp.ObjClust))))+1,...
             (ones(size(cell2mat(HypCtrl.ObjClust))))+3,...
             (ones(size(cell2mat(HypCasp.ObjClust))))+4]; 

boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Number of Localizations per Multi-Exposure'})
ylabel({'Number of Localizations per Multi-Exposure'})
set(aObjClust,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aObjClust,'LineWidth',3,'FontSize',18)

%density of objects per cluster. Multi-exposure localization density (nm^-2)
figure(8)
aDensObjClust=axes;

figurematrix=[cell2mat(YstCtrl.DensObjClust),...
             cell2mat(YstCasp.DensObjClust),...
             cell2mat(HypCtrl.DensObjClust),...
             cell2mat(HypCasp.DensObjClust)];

figuregroup=[(ones(size(cell2mat(YstCtrl.DensObjClust)))),...
             (ones(size(cell2mat(YstCasp.DensObjClust))))+1,...
             (ones(size(cell2mat(HypCtrl.DensObjClust))))+3,...
             (ones(size(cell2mat(HypCasp.DensObjClust))))+4]; 

boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Multi-exposure localization density (nm^-^2)'})
ylabel({'Multi-exposure Localization Density (nm^-^2)'})
set(aDensObjClust,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aDensObjClust,'LineWidth',3,'FontSize',18)

%Density of singlet glucan exposure per ROIs
figure(9)
aSingDens=axes;

figurematrix=[cell2mat(YstCtrl.SingleDensity')',...
    cell2mat(YstCasp.SingleDensity')',...
    cell2mat(HypCtrl.SingleDensity')',...
    cell2mat(HypCasp.SingleDensity')'];

figuregroup=[(ones(size(cell2mat(YstCtrl.SingleDensity')'))),...
    (ones(size(cell2mat(YstCasp.SingleDensity')'))+1),...
     (ones(size(cell2mat(HypCtrl.SingleDensity')'))+3),...
     (ones(size(cell2mat(HypCasp.SingleDensity')'))+4)]; 
 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Density of Singlet-Glucan Exposure'})
ylabel({'Singlet-Glucan Exposure Density (nm^-^2)'})
set(aSingDens,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aSingDens,'LineWidth',3,'FontSize',18)

%Density of singlet glucan exposure per ROIs
figure(10)
aMulDens=axes;

figurematrix=[cell2mat(YstCtrl.MultipleDensity')',...
    cell2mat(YstCasp.MultipleDensity')',...
    cell2mat(HypCtrl.MultipleDensity')',...
    cell2mat(HypCasp.MultipleDensity')'];

figuregroup=[(ones(size(cell2mat(YstCtrl.MultipleDensity')'))),...
    (ones(size(cell2mat(YstCasp.MultipleDensity')'))+1),...
     (ones(size(cell2mat(HypCtrl.MultipleDensity')'))+3),...
     (ones(size(cell2mat(HypCasp.MultipleDensity')'))+4)]; 
 
boxplot(figurematrix, figuregroup,'colors',[0 0 0],'symbol','k+')
box off
title({'Density of Multi-Glucan Exposure'})
ylabel({'Multi-Glucan Exposure Density (nm^-^2)'})
set(aMulDens,'color', [1 1 1],'XTickLabel',{'Y','Y','H','H'},...
    'XTick',[1 2 3 4])
set(gcf,'color',[1 1 1])
set(aMulDens,'LineWidth',3,'FontSize',18)


%% 
%save multiple figures to the same folder
savedir = '/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Fig5ClusteringRev';
cd(savedir)

saveas(figure(1),'NumColapsed.fig')
saveas(figure(1),'NumColapsed.tif')

saveas(figure(2),'RadiiMultiples.fig')
saveas(figure(2),'RadiiMultiples.tif')

saveas(figure(3),'DistSingleMultiple.fig')
saveas(figure(3),'DistSingleMultiple.tif')

saveas(figure(4),'SingleperTotal.fig')
saveas(figure(4),'SingleperTotal.tif')

saveas(figure(5),'DensLocal.fig')
saveas(figure(5),'DensLocal.tif')

saveas(figure(6),'DensExp.fig')
saveas(figure(6),'DensExp.tif')

saveas(figure(7),'Obj_Cluster.fig')
saveas(figure(7),'Obj_Cluster.tif')

saveas(figure(8),'DensObj_Cluster.fig')
saveas(figure(8),'DensObj_Cluster.tif')

saveas(figure(9),'SingleDens_ROI.fig')
saveas(figure(9),'SingleDens_ROI.tif')

saveas(figure(10),'MultiDens_ROI.fig')
saveas(figure(10),'MultiDens_ROI.tif')

%% Whitney test
close all
clear
clc

% Yst Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl')
load('YstCtrlSTATSDBSCAN.mat')

% Yst Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp')
load('YstCaspSTATSDBSCAN.mat')

% Hyp Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl')
load('HypCtrlSTATSDBSCAN.mat')

% Hyp Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp')
load('HypCaspSTATSDBSCAN.mat')

% Radii of multiple
[MeanRadi.YY]=ranksum(cellfun(@median,YstCtrl.MeanRadius2plus),cellfun(@median,YstCasp.MeanRadius2plus));
[MeanRadi.YH]=ranksum(cellfun(@median,YstCtrl.MeanRadius2plus),cellfun(@median,HypCtrl.MeanRadius2plus));
[MeanRadi.HH]=ranksum(cellfun(@median,HypCtrl.MeanRadius2plus),cellfun(@median,HypCasp.MeanRadius2plus));

%Distance of single+multiple
[MeanDist.YY]=ranksum(cellfun(@median,YstCtrl.MeanDist),...
                              cellfun(@median,YstCasp.MeanDist));
[MeanDist.YH]=ranksum(cellfun(@median,YstCtrl.MeanDist),...
                              cellfun(@median,HypCtrl.MeanDist));
[MeanDist.HH]=ranksum(cellfun(@median,HypCtrl.MeanDist),...
                              cellfun(@median,HypCasp.MeanDist));
                          
%Density of locolizations                          
[DensLoc.YY]=ranksum(cellfun(@median,YstCtrl.NumObjs),cellfun(@median,YstCasp.NumObjs));                          
[DensLoc.YH]=ranksum(cellfun(@median,YstCtrl.NumObjs),cellfun(@median,HypCtrl.NumObjs));                            
[DensLoc.HH]=ranksum(cellfun(@median,HypCtrl.NumObjs),cellfun(@median,HypCasp.NumObjs));   

%Density of exposure
[DensExp.YY]=ranksum(cellfun(@median,[YstCtrl.SingleDensity,YstCtrl.MultipleDensity]),...
                              cellfun(@median,[YstCasp.SingleDensity,YstCasp.MultipleDensity]));
[DensExp.YH]=ranksum(cellfun(@median,[YstCtrl.SingleDensity,YstCtrl.MultipleDensity]),...
                              cellfun(@median,[HypCtrl.SingleDensity,HypCtrl.MultipleDensity]));
[DensExp.HH]=ranksum(cellfun(@median,[HypCtrl.SingleDensity,HypCtrl.MultipleDensity]),...
                              cellfun(@median,[HypCasp.SingleDensity,HypCasp.MultipleDensity]));

%localizations/multi-gluca exposure
[DensLocClust.YY]=ranksum(cellfun(@median,YstCtrl.DensObjClust),cellfun(@median,YstCasp.DensObjClust));                          
[DensLocClust.YH]=ranksum(cellfun(@median,YstCtrl.DensObjClust),cellfun(@median,HypCtrl.DensObjClust));                            
[DensLocClust.HH]=ranksum(cellfun(@median,HypCtrl.DensObjClust),cellfun(@median,HypCasp.DensObjClust));   

%singlet-exposure localization density, singlets/ROI
[DensLocSingle.YY]=ranksum(cellfun(@median,YstCtrl.SingleDensity),cellfun(@median,YstCasp.SingleDensity));                          
[DensLocSingle.YH]=ranksum(cellfun(@median,YstCtrl.SingleDensity),cellfun(@median,HypCtrl.SingleDensity));                            
[DensLocSingle.HH]=ranksum(cellfun(@median,HypCtrl.SingleDensity),cellfun(@median,HypCasp.SingleDensity)); 

%multi-exposure localization density, multiples/ROI
[DensLocmulti.YY]=ranksum(cellfun(@median,YstCtrl.MultipleDensity),cellfun(@median,YstCasp.MultipleDensity));                          
[DensLocmulti.YH]=ranksum(cellfun(@median,YstCtrl.MultipleDensity),cellfun(@median,HypCtrl.MultipleDensity));                            
[DensLocmulti.HH]=ranksum(cellfun(@median,HypCtrl.MultipleDensity),cellfun(@median,HypCasp.MultipleDensity)); 

%Number of localization per multi-exposure
[NumLocClust.YY]=ranksum(cellfun(@median,YstCtrl.ObjClust),cellfun(@median,YstCasp.ObjClust));                          
[NumLocClust.YH]=ranksum(cellfun(@median,YstCtrl.ObjClust),cellfun(@median,HypCtrl.ObjClust));                            
[NumLocClust.HH]=ranksum(cellfun(@median,HypCtrl.ObjClust),cellfun(@median,HypCasp.ObjClust));   
                                                  
%Singlets fraction
[SingletFrac.YY]=ranksum(cellfun(@median,YstCtrl.singlesFraction),cellfun(@median,YstCasp.singlesFraction));                          
[SingletFrac.YH]=ranksum(cellfun(@median,YstCtrl.singlesFraction),cellfun(@median,HypCtrl.singlesFraction));                            
[SingletFrac.HH]=ranksum(cellfun(@median,HypCtrl.singlesFraction),cellfun(@median,HypCasp.singlesFraction));                             
%% Ripley's, new analysis, Yst Ctrl
close all
clear
clc

% Yst Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl')
load('YstCtrlSTATSDBSCAN.mat')

nROIs=size(YstCtrl.Clusts,2);
for ii=1:nROIs
    
    cutoff=1000;
    rate=200;
    area=1000000;
    % experimental data
    ROIpts=[YstCtrl.XYcollaped{1,ii}(:,1),YstCtrl.XYcollaped{1,ii}(:,2)];
    YstCtrlRip.K{ii}=SR_Ripley_Jia(cutoff,rate, area,ROIpts)';
    % simulation
    points=rand(size(ROIpts))*sqrt(area);
    YstCtrlRip.KSim{ii}=SR_Ripley_Jia(cutoff,rate, area,points)';
    
    clear  points cutoff rate area ROIpts
end
clear ii ROIpts points cutoff rate area nROIs YstCtrlROIs
close all
save YstCtrlRip

%% Ripley's, Yst Casp
close all
clear
clc

% Yst Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp')
load('YstCaspSTATSDBSCAN.mat')

nROIs=size(YstCasp.Clusts,2);
for ii=1:nROIs
    
    cutoff=1000;
    rate=200;
    area=1000000;
    % experimental data
    ROIpts=[YstCasp.XYcollaped{1,ii}(:,1),YstCasp.XYcollaped{1,ii}(:,2)];
    YstCaspRip.K{ii}=SR_Ripley_Jia(cutoff,rate, area,ROIpts)';
    % simulation
    points=rand(size(ROIpts))*sqrt(area);
    YstCaspRip.KSim{ii}=SR_Ripley_Jia(cutoff,rate, area,points)';
    
    clear  points cutoff rate area ROIpts
end
clear ii ROIpts points cutoff rate area nROIs YstCaspROIs
close all
save YstCaspRip

%% Ripley's, Hyp Ctrl
close all
clear
clc

% Yst Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl')
load('HypCtrlSTATSDBSCAN.mat')

nROIs=size(HypCtrl.Clusts,2);
for ii=1:nROIs
    
    cutoff=1000;
    rate=200;
    area=1000000;
    % experimental data
    ROIpts=[HypCtrl.XYcollaped{1,ii}(:,1),HypCtrl.XYcollaped{1,ii}(:,2)];
    HypCtrlRip.K{ii}=SR_Ripley_Jia(cutoff,rate, area,ROIpts)';
    % simulation
    points=rand(size(ROIpts))*sqrt(area);
    HypCtrlRip.KSim{ii}=SR_Ripley_Jia(cutoff,rate, area,points)';
    
    clear  points cutoff rate area ROIpts
end
clear ii ROIpts points cutoff rate area nROIs HypCtrlROIs
close all
save HypCtrlRip

%% Ripley's, Hyp Casp
close all
clear
clc

% Yst Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp')
load('HypCaspSTATSDBSCAN.mat')

nROIs=size(HypCasp.Clusts,2);
for ii=1:nROIs
    
    cutoff=1000;
    rate=200;
    area=1000000;
    % experimental data
    ROIpts=[HypCasp.XYcollaped{1,ii}(:,1),HypCasp.XYcollaped{1,ii}(:,2)];
    HypCaspRip.K{ii}=SR_Ripley_Jia(cutoff,rate, area,ROIpts)';
    % simulation
    
    points=rand(size(ROIpts))*sqrt(area);
    HypCaspRip.KSim{ii}=SR_Ripley_Jia(cutoff,rate, area,points)';
    
    clear  points cutoff rate area ROIpts
end
clear ii ROIpts points cutoff rate area nROIs HypCaspROIs
close all
save HypCaspRip

%% combined STATS of Ripley's
close all
clear
clc

% Yst Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCtrl')
load('YstCtrlRip.mat')

% Yst Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/YstCasp')
load('YstCaspRip.mat')

% Hyp Ctrl
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCtrl')
load('HypCtrlRip.mat')

% Hyp Casp
cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/HypCasp')
load('HypCaspRip.mat')

%Yst Ctrl
figure(1)
KYCtrl=axes;
plot(0:5:1000,mean(cell2mat(YstCtrlRip.K),2),'r','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(YstCtrlRip.K),2),std(cell2mat(YstCtrlRip.K),0,2),...
         'r','linewidth',1)
hold on
plot(0:5:1000,mean(cell2mat(YstCtrlRip.KSim),2),'b','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(YstCtrlRip.KSim),2),std(cell2mat(YstCtrlRip.KSim),0,2),...
         'b','linewidth',1)
hold off
box off
title({'YCtrl'})
xlim([0 1100])
ylim([0 inf])
ylabel({'K(r)'})
xlabel({'r(nm)'})
legend({'Y -','', 'Random',''},'location','northwest')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(KYCtrl,'LineWidth',3,'FontSize',18)

%Yst Casp
figure(2)
KYCasp=axes;
plot(0:5:1000,mean(cell2mat(YstCaspRip.K),2),'r','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(YstCaspRip.K),2),std(cell2mat(YstCaspRip.K),0,2),...
         'r','linewidth',1)
hold on
plot(0:5:1000,mean(cell2mat(YstCaspRip.KSim),2),'b','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(YstCaspRip.KSim),2),std(cell2mat(YstCaspRip.KSim),0,2),...
         'b','linewidth',1)
hold off
box off
title({'YCasp'})
xlim([0 1100])
ylabel({'K(r)'})
xlabel({'r(nm)'})
legend({'Y +','', 'Random',''},'location','northwest')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(KYCasp,'LineWidth',3,'FontSize',18)

%Hyp Ctrl
figure(3)
KHCtrl=axes;
plot(0:5:1000,mean(cell2mat(HypCtrlRip.K),2),'r','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(HypCtrlRip.K),2),std(cell2mat(HypCtrlRip.K),0,2),...
         'r','linewidth',1)
hold on
plot(0:5:1000,mean(cell2mat(HypCtrlRip.KSim),2),'b','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(HypCtrlRip.KSim),2),std(cell2mat(HypCtrlRip.KSim),0,2),...
         'b','linewidth',1)
hold off
box off
title({'HCtrl'})
xlim([0 1100])
ylabel({'K(r)'})
xlabel({'r(nm)'})
legend({'H -','', 'Random',''},'location','northwest')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(KHCtrl,'LineWidth',3,'FontSize',18)

%Hyp Casp
figure(4)
KHCasp=axes;
plot(0:5:1000,mean(cell2mat(HypCaspRip.K),2),'r','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(HypCaspRip.K),2),std(cell2mat(HypCaspRip.K),0,2),...
         'r','linewidth',1)
hold on
plot(0:5:1000,mean(cell2mat(HypCaspRip.KSim),2),'b','linewidth',3)
hold on
errorbar(0:5:1000,mean(cell2mat(HypCaspRip.KSim),2),std(cell2mat(HypCaspRip.KSim),0,2),...
         'b','linewidth',1)
hold off
box off
title({'HCasp'})
xlim([0 1100])
ylabel({'K(r)'})
xlabel({'r(nm)'})
legend({'H +','', 'Random',''},'location','northwest')
legend BOXOFF
set(gcf,'color',[1 1 1])
set(KHCasp,'LineWidth',3,'FontSize',18)


%% 
%save multiple figures to the same folder
savedir = '/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/FigRipley';
cd(savedir)

saveas(figure(1),'YstCtrl.fig')
saveas(figure(1),'YstCtrl.tif')

saveas(figure(2),'YstCasp.fig')
saveas(figure(2),'YstCasp.tif')

saveas(figure(3),'HypCtrl.fig')
saveas(figure(3),'HypCtrl.tif')

saveas(figure(4),'HypCasp.fig')
saveas(figure(4),'HypCasp.tif')

saveas(figure(5),'YCtrlP.fig')
saveas(figure(5),'YCtrlP.tif')

saveas(figure(6),'YCaspP.fig')
saveas(figure(6),'YCaspP.tif')

saveas(figure(7),'HCtrlP.fig')
saveas(figure(7),'HCtrlP.tif')

saveas(figure(8),'HCaspP.fig')
saveas(figure(8),'HCaspP.tif')
%% Single probe on glass
clear 
close all
clc

pixel2nm = 16000/150;%convert pixel size into nm
x_size=3000;
y_size=3000;
DataDirs={...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass/G1/Results4',...% Jia's
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass/G1/Results4',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass/G1/Results4',...% has 3 ROIs
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass/G2/Results4',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass/G2/Results4',...
'/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass/G2/Results4',...% has 3 ROIs
};
% save ROI for each cell
for ii=1:length(DataDirs)
    DataDir=DataDirs{ii};
    cd(DataDir)
    load('Fit4.mat')
    
    X = double(SRD.Results.X) .* pixel2nm;   % nm
    Y = double(SRD.Results.Y) .* pixel2nm;   % nm
    
    X_STD = double(SRD.Results.X_STD) .* pixel2nm;   % nm
    Y_STD = double(SRD.Results.Y_STD) .* pixel2nm;   % nm

    [index_ROI,ROI,n_ROIs] = get_ROI(X, Y, x_size, y_size);% hit Q on keyboard to quit
    
    GlassROIs.Xs{ii}=X(index_ROI{n_ROIs});
    GlassROIs.Ys{ii}=Y(index_ROI{n_ROIs});
    
    GlassROIs.X_STDs{ii}=X_STD(index_ROI{n_ROIs});
    GlassROIs.Y_STDs{ii}=Y_STD(index_ROI{n_ROIs});
    
    close all
end

Sigma_Reg = [0 0]; % no driftcorrection
nROIs=size(GlassROIs.Xs,2);

for ii=1:nROIs
    close all
SRc = SRcluster;
SRc.A_ROI = 1.0e+06; % nm^2
SRc.clusterSR([GlassROIs.Xs{ii}, GlassROIs.Ys{ii}], [GlassROIs.X_STDs{ii}, GlassROIs.Y_STDs{ii}], Sigma_Reg);
SRc.E = 50;
SRc.Algorithm = 'DBSCAN_Daszykowski';   
SRc.minPts = 3;
SRclusterFig = SRc.plotSRclusters();
SRc.PlotFigures=false;
[results, analysisFigs] = SRc.analyzeSRclusters();
% 
% Glass.XYcollaped{ii}=SRc.XY;
% 
% Glass.Clust3plusNum{ii} = SRc.cutoffPlots();

Glass.NumColapsed{ii}=results.numobjs_collapsed;% number of object collapsed

% Glass.SingleDensity{ii}=results.numdensity(1);% number of cluster(single,double,multiple)/nm^2
% 
% Glass.MultipleDensity{ii}=results.numdensity(3);
% 
% Glass.MeanRadius2plus{ii}=results.radius_per_2orMore_cluster;% radii of 2 or 3+
% 
% Glass.MeanSingleDist{ii}=results.min_c2c_dist(1);% mean distance of point-point(single)
% 
% Glass.MeanMultipleDist{ii}=results.min_c2c_dist(3);%mean distance of centroid of 3+ clusters
% 
% Glass.NumObjs{ii}=results.numobjs;
% Glass.Clusts{ii}=results.numclust;
end


cd('/Users/jialin/Research/ExpData/Manuscript/SR3rd-Revision/Glass')
save('GlassROIs','GlassROIs')
save('GlassSTATSDBSCAN','Glass')

figure
aGlas=axes;
boxplot(cell2mat(Glass.NumColapsed),'colors',[0 0 0],'symbol','k+')
box off
title({'Glass'})
ylabel({'Number of Observation Collapsed'})
xlabel({''})
set(gcf,'color',[1 1 1])
set(aGlas,'LineWidth',3,'FontSize',18)

saveas(figure(1),'GlasColapsed.fig')
saveas(figure(1),'GlasColapsed.tif')

%% 
% do t test for each sample, probablities vs. r
% [~,ProbYCtrl]=ttest2(cell2mat(YstCtrlRip.K),cell2mat(YstCtrlRip.KSim),'dim',2);
% [~,ProbYCasp]=ttest2(cell2mat(YstCaspRip.K),cell2mat(YstCaspRip.KSim),'dim',2);
% [~,ProbHCtrl]=ttest2(cell2mat(HypCtrlRip.K),cell2mat(HypCtrlRip.KSim),'dim',2);
% [~,ProbHCasp]=ttest2(cell2mat(HypCaspRip.K),cell2mat(HypCaspRip.KSim),'dim',2);
% 
% figure(5)
% aPYCtrl=axes;
% plot(0:5:1000,ProbYCtrl,'b','linewidth',3)
% box off
% title({'YCtrl'})
% xlim([0 1100])
% ylabel({'Probabilities'})
% xlabel({'r(nm)'})
% set(gcf,'color',[1 1 1])
% set(aPYCtrl,'LineWidth',3,'FontSize',18)
% 
% figure(6)
% aPYCasp=axes;
% plot(0:5:1000,ProbYCasp,'b','linewidth',3)
% box off
% title({'YCasp'})
% xlim([0 1100])
% ylabel({'Probabilities'})
% xlabel({'r(nm)'})
% set(gcf,'color',[1 1 1])
% set(aPYCasp,'LineWidth',3,'FontSize',18)
% 
% figure(7)
% aPHCtrl=axes;
% plot(0:5:1000,ProbHCtrl,'b','linewidth',3)
% box off
% title({'HCtrl'})
% xlim([0 1100])
% ylabel({'Probabilities'})
% xlabel({'r(nm)'})
% set(gcf,'color',[1 1 1])
% set(aPHCtrl,'LineWidth',3,'FontSize',18)
% 
% figure(8)
% aPHCasp=axes;
% plot(0:5:1000,ProbHCasp,'b','linewidth',3)
% box off
% title({'HCasp'})
% xlim([0 1100])
% ylabel({'Probabilities'})
% xlabel({'r(nm)'})
% set(gcf,'color',[1 1 1])
% set(aPHCasp,'LineWidth',3,'FontSize',18)



















