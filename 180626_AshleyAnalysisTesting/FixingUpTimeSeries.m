disttitles = {'Mean Intensity <50 \mum from Soma','Mean Intensity 50-100 \mum from Soma',...
    'Mean Intensity 100-200 \mum from Soma'};
savestr = {'50','50-100','100-200'};
listdir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\180626 Results and Info';

%% fixing up files
% load the data and group together
savetypestr = 'GluA1';
% filelistname = 'all_NL1_global.mat';
filelistname = 'Ashley062018_all_GluA1_Global_allMs.mat';
load(fullfile(listdir,filelistname));
% savedir = 'E:\Sam\Data\MJK_zapERtrap_for_sam\180626 Results and Info\';
GluA1means = [];
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
GluA1means(dd).mean = nanmean(newarray,1);
GluA1means(dd).ste = nanstd(newarray,1)./sqrt(size(newarray,1));
% GluA1means(dd).ste = nanstd(newarray,1);

GluA1means(dd).x = xvals;
% figure; plot(allintensities'); hold on;
f = figure; plot(xvals,newarray','Color',[0 0.7 0.7]); hold on;
plot(xvals,GluA1means(dd).mean,'Color','k','LineWidth',2)
xlim([-5 100]);
set(gca,'FontSize',16)
title(disttitles{dd},'FontSize',16)
xlabel('Time (min)','FontSize',16);
ylabel('Cargo Intensity Density (AU)','FontSize',16);
saveas(f,fullfile(listdir,[savetypestr '-' savestr{dd}]),'fig');
saveas(f,fullfile(listdir,[savetypestr '-' savestr{dd}]),'png');
close(f)
% ylim([0 100]);
end

%%
% load the data and group together
savetypestr = 'NL1';

filelistname = 'all_NL1_global_allMs.mat';
NL1means = [];

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
NL1means(dd).mean = nanmean(newarray,1);
NL1means(dd).ste = nanstd(newarray,1)./sqrt(size(newarray,1));
% NL1means(dd).ste = nanstd(newarray,1);

NL1means(dd).x = xvals;

% figure; plot(allintensities'); hold on;
f = figure; plot(xvals,newarray','Color',[0.7 0 0.7]); hold on;
plot(xvals,NL1means(dd).mean,'Color','k','LineWidth',2)
title(disttitles{dd},'FontSize',16)
set(gca,'FontSize',16);
xlabel('Time (min)','FontSize',16);
ylabel('Cargo Intensity Density (AU)','FontSize',16);
xlim([-5 100]);
saveas(f,fullfile(listdir,[savetypestr '-' savestr{dd}]),'fig');
saveas(f,fullfile(listdir,[savetypestr '-' savestr{dd}]),'png');
close(f)
end
%% tfr data
% load the data and group together
savetypestr = 'TfR';

filelistname = 'Ashley062018_all_TfR_Global_allMs.mat';
TfRmeans = [];


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

thesevals = allintensities;
thesevals(:,1) = [];
newarray(:,1:size(thesevals,2) )  = thesevals;

xvalspre = -4:0;
xvalspost = [1:(size(newarray,2)-size(xvalspre,2))]*2;
xvals = [xvalspre,xvalspost];
TfRmeans(dd).mean = nanmean(newarray,1);
TfRmeans(dd).ste = nanstd(newarray,1)./sqrt(size(newarray,1));
% TfRmeans(dd).ste = nanstd(newarray,1);

TfRmeans(dd).x = xvals;

% figure; plot(allintensities'); hold on;
f = figure; plot(xvals,newarray','Color',[0.7 0.7 0]); hold on;
plot(xvals,TfRmeans(dd).mean,'Color','k','LineWidth',2)
title(disttitles{dd},'FontSize',16);
set(gca,'FontSize',16);
xlabel('Time (min)','FontSize',16);
ylabel('Cargo Intensity Density (AU)','FontSize',16);
xlim([-5 100]);
saveas(f,fullfile(listdir,[savetypestr '-' savestr{dd}]),'fig');
saveas(f,fullfile(listdir,[savetypestr '-' savestr{dd}]),'png');
close(f)

end

%% 
% plot means with stderr
savetypestr = 'MeanVals';
legstr = {'GluA1','NL1','TfR'};
cols = lines();
for mm = 1:3
    f = figure; hold on;
    % GluA1
    e1 = errorbar(GluA1means(mm).x,GluA1means(mm).mean,GluA1means(mm).ste,'Color',cols(1,:));
    e1.CapSize = 0;
    plot(GluA1means(mm).x,GluA1means(mm).mean,'LineWidth',2,'Color',cols(1,:),'DisplayName',legstr{1});
    % NL1
    e2 = errorbar(NL1means(mm).x,NL1means(mm).mean,NL1means(mm).ste,'Color',cols(2,:));
    e2.CapSize = 0;
    plot(NL1means(mm).x,NL1means(mm). mean,'LineWidth',2,'Color',cols(2,:),'DisplayName',legstr{2});
    
    e3 = errorbar(TfRmeans(mm).x,TfRmeans(mm).mean,TfRmeans(mm).ste,'Color',cols(3,:));
     e3.CapSize = 0;
    plot(TfRmeans(mm).x,TfRmeans(mm).mean,'LineWidth',2,'Color',cols(3,:),'DisplayName',legstr{3});
    set(gca,'FontSize',16)
    set(gca,'xlim',[-5 100])
    title(disttitles{dd},'FontSize',16)
    xlabel('Time (min)','FontSize',16);
    ylabel('Cargo Intensity Density (AU)','FontSize',16);
    legend('show'); 
    legend('Location','nw')
    saveas(f,fullfile(listdir,[savetypestr '-' savestr{mm}]),'fig');
    saveas(f,fullfile(listdir,[savetypestr '-' savestr{mm}]),'png');
    
    set(gca,'xlim',[-5 62])
    saveas(f,fullfile(listdir,[savetypestr '-' savestr{mm} '-TfRscale']),'fig');
    saveas(f,fullfile(listdir,[savetypestr '-' savestr{mm} '-TfRscale']),'png');
    close(f)
end

