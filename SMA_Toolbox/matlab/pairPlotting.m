tsaDir='C:\OHSU\Analysis2015\2015-02-17-pairs\TrackSegmentAnalysis';

% tsaFile='0121C01 BQ P_MMStack.ome-ROI1.tsa';
% pair=[1 2; 5 1];

% tsaFile='0207C02 68m P_MMStack.ome-ROI1.tsa';
% pair=[6 1; 7 1];

% tsaFile='0217C01 50m uC_MMStack.ome-ROI1.tsa';
% pair=[1 1; 3 1];

tsaFile='0217C01 55m P_MMStack.ome-ROI1.tsa';
pair=[2 1; 3 1];


[~,tsaBase,~]=fileparts(tsaFile)
tsaPath=fullfile(tsaDir,tsaFile);
tsa=TrackSegmentAnalysis(tsaPath);
% trackID1, segID1; trackID2, segID2
t1=tsa.trackLocalizationsTable(pair(1,1),pair(1,2));
t2=tsa.trackLocalizationsTable(pair(2,1),pair(2,2));
frameT=tsa.frameT;

minFrame = max(t1.FrameIdx(1), t2.FrameIdx(2));
maxFrame = min(t1.FrameIdx(end), t2.FrameIdx(end));
Ts = (minFrame:maxFrame).*frameT;
d = sqrt((interp1(t1.t,t1.x,Ts)- interp1(t2.t,t2.x,Ts)).^2 + (interp1(t1.t,t1.y,Ts)- interp1(t2.t,t2.y,Ts)).^2);

f=figure();
plot(Ts,d);
xlabel('Time (s)');
ylabel('Distance (um)');
title(tsaFile,'interpreter', 'none');
saveFileName = sprintf('%s-pair-distance.fig',tsaBase);
savePath=fullfile(tsaDir,'Pairs');
if ~exist(savePath,'dir')
    mkdir(savePath);
end
savePath=fullfile(tsaDir,'Pairs',saveFileName);
savefig(f,savePath);

f=tsa.view3DSegmentSequence(pair(:,1),pair(:,2));
saveFileName = sprintf('%s-pair-plot.fig',tsaBase);
savePath=fullfile(tsaDir,'Pairs',saveFileName);
savefig(f,savePath);

