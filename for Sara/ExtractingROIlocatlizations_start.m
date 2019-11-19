% inputs:
% - SR .mat file with localizations
% - image j ROIs
% - scale factor (5 for gaussian or 10 for pixel)
% - save directory
load('D:\WDKennedyLabHDDBackup\Projects\Project Cry2Olig-Gephyrin\SR\190815 Gaba_Geph_RIM\results\GIB GabaRIM\SMLM_GIB_GABA_RIM+light_GABA647_2_SR_GIB_GABA_RIM+light_RIM568_2_SR.mat');
roiname = 'D:\WDKennedyLabHDDBackup\Projects\Project Cry2Olig-Gephyrin\SR\190815 Gaba_Geph_RIM\results\images\GIB+light_GABA647_RIM568_2\RoiSet.zip';
savedir = 'D:\WDKennedyLabHDDBackup\Projects\Project Cry2Olig-Gephyrin\SR\190815 Gaba_Geph_RIM\results\images\GIB+light_GABA647_RIM568_2_Stuff';
% load in value
X_ch1 = data.SML_data.ch1.position_x;
Y_ch1 = data.SML_data.ch1.position_y;

X_ch2 = data.SML_data.ch2.position_x;
Y_ch2 = data.SML_data.ch2.position_y;


% load in ROIS
[sROI] = ReadImageJROI(roiname);
figure; nrow = ceil(numel(sROI)/2); ncol = 2;
scalefactor = 5;
for ii = 1:numel(sROI)
    % get bounds from sROI
    if numel(sROI) == 1
        ROIrect = sROI.vnRectBounds;
    else
        ROIrect = sROI{ii}.vnRectBounds;
    end
    % find x-y min/max
    ymin = ROIrect(1)*scalefactor;
    ymax = ROIrect(3)*scalefactor;
    xmin = ROIrect(2)*scalefactor;
    xmax = ROIrect(4)*scalefactor;
    % get all points within bounds for ch1
    x_ch1_id = X_ch1(:,1) >= xmin & X_ch1(:,1) < xmax;
    y_ch1_id = Y_ch1(:,1) >= ymin & Y_ch1(:,1) < ymax;
    goodids_ch1 = x_ch1_id & y_ch1_id;
    Xs_ch1 = X_ch1(goodids_ch1);
    Ys_ch1 = Y_ch1(goodids_ch1);
    % get all points within bounds for ch2
    x_ch2_id = X_ch2(:,1) >= xmin & X_ch2(:,1) < xmax;
    y_ch2_id = Y_ch2(:,1) >= ymin & Y_ch2(:,1) < ymax;
    goodids_ch2 = x_ch2_id & y_ch2_id;
    Xs_ch2 = X_ch2(goodids_ch2);
    Ys_ch2 = Y_ch2(goodids_ch2);
    % plot points as scatter plot
    h = subplot(nrow,ncol,ii); hold on;
%     h = figure; hold on;
    scatter(h,Xs_ch2,-Ys_ch2, '.m');
    scatter(h,Xs_ch1,-Ys_ch1, '.c');
    title(['ROI # ' num2str(ii)]);
    input.E = 20;
    input.minPts = 10;
    input.algorithm = 'DBSCAN_Tran';
    input.savestr = num2str(ii);
    [c1,c2,f1,f2] = plotUNMClustering(Xs_ch1,Ys_ch1,Xs_ch2,Ys_ch2,input,savedir);
    close(f1); close(f2);
end