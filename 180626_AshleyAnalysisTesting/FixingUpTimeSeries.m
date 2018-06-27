%% fixing up files
% load the data and group together
listdir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\180626 Results and Info';
% filelistname = 'all_NL1_global.mat';
filelistname = 'Ashley062018_all_GluA1_Global_allMs.mat';
load(fullfile(listdir,filelistname));
% savedir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\180626 Results and Info\';
for dd = 1:3
distval = dd;
allintensities = nan(numel(allMs),1);
% range = 1:10;
for m = 1:numel(allMs)
    vals = allMs{m}.areanormintensity(distval,:);
    if size(vals,2) > size(allintensities,2)
        addonsize =  size(vals,2) - size(allintensities,2);
       allintensities = cat(2,allintensities,nan(numel(allMs),addonsize));
    end
    allintensities(m,1:size(vals,2)) = vals;
end

% fix up allintensities
% preframes = 1:6; 


newarray = nan(size(allintensities));
thesevals = allintensities(1:9,:);
thesevals(:,1:6) = [];
newarray(1:9,1:size(thesevals,2) )  = thesevals;

clear thesevals
thesevals = allintensities(10:end,:);
blankcols = nan(size(thesevals,1),1);
caton = [thesevals(:,4),blankcols,thesevals(:,5),blankcols,thesevals(:,6)];
thesevals = [caton, thesevals(:,7:end)];
newarray(10:end,1:size(thesevals,2)) = thesevals;

newarray(12,:) = [];
newarray(14,:) = [];

xvalspre = -4:0;
xvalspost = [1:(size(newarray,2)-size(xvalspre,2))]*2;
xvals = [xvalspre,xvalspost];
% figure; plot(allintensities'); hold on;
figure; plot(xvals,newarray','Color',[0 0.7 0.7]); hold on;
plot(xvals,nanmean(newarray,1),'Color','k','LineWidth',2)
xlim([-5 100]);
% ylim([0 100]);
end

%%
% load the data and group together
listdir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\180626 Results and Info';
filelistname = 'all_NL1_global_allMs.mat';
% filelistname = 'Ashley062018_all_GluA1_Global.mat';
load(fullfile(listdir,filelistname));
for dd = 1:3
distval = dd;
allintensities = nan(numel(allMs),1);
% range = 1:10;
for m = 1:numel(allMs)
    vals = allMs{m}.areanormintensity(distval,:);
    if size(vals,2) > size(allintensities,2)
        addonsize =  size(vals,2) - size(allintensities,2);
       allintensities = cat(2,allintensities,nan(numel(allMs),addonsize));
    end
    allintensities(m,1:size(vals,2)) = vals;
end

% fix up allintensities
% preframes = 1:6; 


newarray = nan(size(allintensities));

thesevals = allintensities(1:3,:);
thesevals(:,1) = [];
newarray(1:3,1:size(thesevals,2) )  = thesevals;

clear thesevals
inc = 4:7;
thesevals = allintensities(inc,:);
thesevals(:,1:2) = [];
newarray(inc,1:size(thesevals,2) )  = thesevals;

clear thesevals
inc = 8:15;
thesevals = allintensities(inc,:);
thesevals(:,1:6) = [];
newarray(inc,1:size(thesevals,2) )  = thesevals;

clear thesevals
inc = 16:18;
thesevals = allintensities(inc,:);
thesevals(:,1:5) = [];
newarray(inc,1:size(thesevals,2) )  = thesevals;

clear thesevals
incr = 19:22;
thesevals = allintensities(incr,:);
thesevals(:,1:10) = [];
newarray(incr,1:size(thesevals,2) )  = thesevals;

clear thesevals
incr = 23:24;
thesevals = allintensities(incr,:);
newarray(incr,1:size(thesevals,2) )  = thesevals;

xvalspre = -4:0;
xvalspost = [1:(size(newarray,2)-size(xvalspre,2))]*2;
xvals = [xvalspre,xvalspost];
% figure; plot(allintensities'); hold on;
figure; plot(xvals,newarray','Color',[0.7 0 0.7]); hold on;
plot(xvals,nanmean(newarray,1),'Color','k','LineWidth',2)
xlim([-5 100]);
end
