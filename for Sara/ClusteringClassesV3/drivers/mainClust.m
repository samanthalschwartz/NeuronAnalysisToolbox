clear all
close all

ROIs = false;
%ROIs = true;

pixel2nm = 16000/150;
x_size = 1000;
y_size = 1000;

infile = '#0001-2015-10-8-14-13-5/#0001-2015-10-8-14-13-5_Results.mat';
%infile = '#0001-2015-10-8-15-6-37/#0001-2015-10-8-15-6-37_Results.mat';
%infile = '#0001-2016-2-2-10-6-25/#0001-2016-2-2-10-6-25_Results.mat';
%infile = '~/local/projects/STMC/Matt/Fibular_Structures/#0001-2015-10-8-14-13-5/#0001-2015-10-8-14-13-5_Results.mat';

c = Clustering();
RT = ROITools();

if ROIs
   [n_ROIs, RoI, Sigma_Reg] = RT.getROI(infile, x_size, y_size);

   X = RoI{1}.X;
   Y = RoI{1}.Y;
   X_STD = RoI{1}.X_STD;
   Y_STD = RoI{1}.Y_STD;
else
   load(infile);

   %SRtest.DriftCorrect = 1;
   %SRtest.FrameConnect = 1;
   %SRtest.Threshold;
   Sigma_Reg = std(SRD.DriftCorrect_XYShift) .* pixel2nm;  %nm
   if isnan(Sigma_Reg)
      Sigma_Reg = [10, 10];
   end

   X = double(SRD.Results.X) .* pixel2nm;   % nm
   Y = double(SRD.Results.Y) .* pixel2nm;   % nm
   X_STD = double(SRD.Results.X_STD) .* pixel2nm;   % nm
   Y_STD = double(SRD.Results.Y_STD) .* pixel2nm;   % nm
end

basefile = regexprep(infile, '_Results.mat', '');

XY = [X, Y];
E = 30;
minPts = 3;
options = 'O';

c.Plotting = true;
c.Alpha = 2;
%c.Valgorithm = [1, 2, 3];
c.Valgorithm = 2;
shrinkFactor = 0.5;
%algorithms = {'DBSCAN_Daszykowski', 'Getis', 'Hierarchical', 'Voronoi'};
%algorithms = {'DBSCAN_Daszykowski', 'Hierarchical', 'Voronoi'};
algorithms = {'Voronoi'};
for i = 1 : length(algorithms)
   algorithm = algorithms{i};
   [nC, C, centers, ptsI] = c.cluster(algorithm, XY, E, minPts);
   fprintf('%s number of clusters = %d\n', algorithm, nC);
   results = c.clusterStats(XY, C, centers, shrinkFactor)
%  save(sprintf('%s_%s.mat', basefile, algorithm), 'results', 'c');
   clusterFig = c.plotClusters(XY, C, centers, ptsI, algorithm, ...
                               shrinkFactor, options);
   showm(clusterFig);
%  saveas(clusterFig, sprintf('%s_%s.fig', basefile, algorithm));
end

SRc = SRcluster();
%SRc.PvalueStatistics = true;
SRc.ShrinkFactor = shrinkFactor;
SRc.clusterSR([X, Y], [X_STD, Y_STD], Sigma_Reg);
%SRc.Algorithm = 'DBSCAN_Daszykowski';   SRc.minPts = minPts;   SRc.E = E;
%SRc.PlotFigures = false;
SRclusterFig = SRc.plotSRclusters();
%saveas(SRclusterFig, sprintf('%s_%s.fig', basefile, 'SRcluster'));
SRcollapseFig = SRc.plotSRcollapse();
%saveas(SRcollapseFig, sprintf('%s_%s.fig', basefile, 'SRcollapse'));
[results, analysisFigs] = SRc.analyzeSRclusters();
%save(sprintf('%s_%s.mat', basefile, 'SRcluster'), 'results', 'SRc');
for i = 1 : length(analysisFigs)
%  saveas(analysisFigs{i}, sprintf('%s_%s%1d.fig', basefile, 'SR', i));
end
