savedir = 'Z:\Sam\IgE SR';
RT = ROITools();
RT.ROI_sizes = [4000, 4000];


clear dat;
dat.files_0ng = {'Z:\Genevieve\Rotation\14-08-01_AF647-IgE\Results\IgE-AF647_well1-Resting_1#0020-2014-8-1-15-1-46.mat',...
    'Z:\Genevieve\Rotation\14-08-01_AF647-IgE\Results\IgE-AF647_well1-Resting_2#0020-2014-8-1-15-13-48.mat'};
dat.files_1ng = {'Z:\Genevieve\Rotation\14-08-01_AF647-IgE\Results\AF647-IgE_well2-0_0001ug-mL_1#0020-2014-8-1-14-0-37.mat',...
    'Z:\Genevieve\Rotation\14-08-01_AF647-IgE\Results\IgE-AF647_well2-0_0001ug-mL_3#0020-2014-8-1-14-41-46.mat'};
dat.files_10ng = {'Z:\Genevieve\Rotation\14-08-01_AF647-IgE\Results\AF647-IgE_well3-0_001ug-mL_1#0020-2014-8-1-13-0-40.mat',...
    'Z:\Genevieve\Rotation\14-08-01_AF647-IgE\Results\AF647-IgE_well3-0_001ug-mL_2#0020-2014-8-1-13-13-36.mat'};
dat.files_100ng = {'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well4-0_1ug-mL_1#0020-2014-7-30-18-24-41.mat',...
    'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well4-0_1ug-mL_2#0020-2014-7-30-18-47-55.mat',...
    'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well4-0_1ug-mL_3#0020-2014-7-30-19-5-42.mat'};
dat.files_1000ng = {'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well5-1ug-mL_1#0020-2014-7-30-17-27-50.mat',...
    'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well5-1ug-mL_2#0020-2014-7-30-17-41-35.mat',...
    'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well5-1ug-mL_3#0020-2014-7-30-18-1-36.mat'};
dat.files_10000ng = {'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well6-10ug-mL_1#0020-2014-7-30-19-27-3.mat',...
    'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well6-10ug-mL_3#0020-2014-7-30-20-0-42.mat',...
    'Z:\Genevieve\Rotation\14-07-30_AF647-IgE\Results\AF647-IgE_well6-10ug-mL_4#0020-2014-7-30-20-19-56.mat'};
cond_names = fieldnames(dat);
clear n_ROIs;
clear RoI;
clear Sigma_Reg;
for ii = 1:numel(cond_names)
    cond = getfield(dat,cond_names{ii});
    for kk = 1:numel(cond)
        infile = cond{kk};
        [a, b, c] = RT.getROI({infile});
        if  kk == 1
            n_ROIs{ii} = a;
            RoI{ii} = b;
            Sigma_Reg{ii} = c;
%             roiname{ii} = cond{kk};
        else  
            n_ROIs{ii} = [n_ROIs{ii} a];
            RoI{ii} = [RoI{ii} b];
            Sigma_Reg{ii} = [Sigma_Reg{ii} c];
%             roiname{ii} = [roiname{ii} cond{kk}];
        end
    end
end
save(fullfile(savedir,'ROIinfo'),'n_ROIs','RoI','Sigma_Reg','dat','RT');
%% now cluster things...
E = 30;
minPts = 3;
options = 'O';
cl = Clustering();
cl.Plotting = true;
shrinkFactor = 0.5;
%algorithms = {'DBSCAN_Daszykowski', 'Getis', 'Hierarchical', 'Voronoi'};
algorithm = 'Hierarchical';
% clear nC;
% clear C;
% clear centers;
% clear ptsI;
clear results;
for ii = 1:numel(cond_names) %loop through conditions
    for jj = 1 : n_ROIs{ii} %loop through all the ROIs associated with this condition
        X = RoI{ii}{jj}.X{1};
        Y = RoI{ii}{jj}.Y{1};
        XY = [X, Y];
        [nC, C, centers, ptsI] = cl.cluster(algorithm, XY, E, minPts);
        results{ii}{jj} = cl.clusterStats(XY, C, centers, shrinkFactor);
        clusterFig = cl.plotClusters(XY, C, centers, ptsI, algorithm, ...
            shrinkFactor, options);
        figure(clusterFig); 
        saveas(gcf,fullfile(savedir,['ClusterPlot_' cond_names{ii} 'ROI_' num2str(jj)]),'png');
        saveas(gcf,fullfile(savedir,['ClusterPlot_' cond_names{ii} 'ROI_' num2str(jj)]),'fig');
        close(gcf);
    end
end
save(fullfile(savedir,'Results'),'results')
%%
clear areas;
for rr = 1:numel(results)
temp = cellmap(@(x) x.areas, results{rr});
areas{rr} = cell2mat(temp);
end
CLarea = figure;
condnames = cellmap(@(x) [strrep(x(7:end),'_','-') '/mL'],cond_names);
plotSpread(areas,'xNames',condnames,'binWidth',.05,'showMM',4)
title('Cluster Area'); ylabel('Cluster Area (nm^2)');

clear equivr;
for rr = 1:numel(results)
temp = cellmap(@(x) x.equiv_radii, results{rr});
equivr{rr} = cell2mat(temp);
end
CLarea = figure;
condnames = cellmap(@(x) [strrep(x(7:end),'_','-') '/mL'],cond_names);
plotSpread(equivr,'xNames',condnames,'binWidth',.05,'showMM',4)
title('Equiv Radii'); ylabel('Equiv Radii of Cluster (nm)');

clear n_pts_per_area;
for rr = 1:numel(results)
temp = cellmap(@(x) x.n_pts_per_area, results{rr});
n_pts_per_area{rr} = cell2mat(temp);
end
n_pts_per_area = figure;
condnames = cellmap(@(x) [strrep(x(7:end),'_','-') '/mL'],cond_names);
plotSpread(equivr,'xNames',condnames,'binWidth',.05,'showMM',4)
title('n-pts-per-area')

clear compactness;
for rr = 1:numel(results)
temp = cellmap(@(x) x.compactness, results{rr});
compactness{rr} = cell2mat(temp);
end
compactness = figure;
condnames = cellmap(@(x) [strrep(x(7:end),'_','-') '/mL'],cond_names);
plotSpread(equivr,'xNames',condnames,'binWidth',.05,'showMM',4)
title('compactness')

clear min_e2e_dists;
for rr = 1:numel(results)
temp = cellmap(@(x) x.min_e2e_dists, results{rr});
min_e2e_dists{rr} = cell2mat(temp);
end
fmin_e2e_dists = figure;
condnames = cellmap(@(x) [strrep(x(7:end),'_','-') '/mL'],cond_names);
plotSpread(min_e2e_dists,'xNames',condnames,'binWidth',.05,'showMM',4)
title('compactness')

clear frac_clust;
for rr = 1:numel(results)
temp = cellmap(@(x) x.n_clustered/x.n_points, results{rr});
frac_clust{rr} = cell2mat(temp);
end
f = figure;
condnames = cellmap(@(x) [strrep(x(7:end),'_','-') '/mL'],cond_names);
plotSpread(frac_clust,'xNames',condnames,'binWidth',.05,'showMM',4)
title('Fraction of localizations in a cluster')
