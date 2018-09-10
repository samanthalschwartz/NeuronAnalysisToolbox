% first load in all the data (this is slowish)
datadir = 'E:\Matt Becker Data (For Review)\SEPfiles\-TeNT';
% files = uipickfiles('Prompt','Pick Files','FilterSpec',datadir);
files_noT = dir(fullfile(datadir,'*.mat'));
alltraces1_noT = cell(numel(files_noT),1);
alltraces2_noT = cell(numel(files_noT),1);
wb = waitbar(0,'Calculating -TenT');
for ff = 1:numel(files_noT)
    clear obj;
    load(fullfile(datadir,files_noT(ff).name));
    obj.calculate_OldROIsintensities();
    disp('Saving file....');
    save(fullfile(datadir,files_noT(ff).name),'obj');
    disp('Done');
    trace1 = obj.oldROI1_sep_sums;
    trace2 = obj.oldROI2_sep_sums;
    alltraces1_noT{ff} = trace1./mean(trace1(:,1:3),2);
    alltraces2_noT{ff} = trace2./mean(trace2(:,1:3),2);
    waitbar(ff/numel(files_noT),wb);
end
close(wb);

datadir = 'E:\Matt Becker Data (For Review)\SEPfiles\+TeNT';
% files = uipickfiles('Prompt','Pick Files','FilterSpec',datadir);
files_T = dir(fullfile(datadir,'*.mat'));
alltraces1_T = cell(numel(files_T),1);
alltraces2_T = cell(numel(files_T),1);
wb = waitbar(0,'Calculating +TenT');
for ff = 1:numel(files_T)
    clear obj;
    load(fullfile(datadir,files_T(ff).name));
    obj.calculate_OldROIsintensities();
    disp('Saving file....');
    save(fullfile(datadir,files_T(ff).name),'obj');
    disp('Done');
    trace1 = obj.oldROI1_sep_sums;
    trace2 = obj.oldROI2_sep_sums;
    if ~isempty(trace1)
        alltraces1_T{ff} = trace1./mean(trace1(:,1:3),2);
    end
    if ~isempty(trace2)
        alltraces2_T{ff} = trace2./mean(trace2(:,1:3),2);
    end
    waitbar(ff/numel(files_T),wb);
end
close(wb);
%%
test1 = cell2mat(alltraces1_T);
hm = SEP.calculate_HeatMap(test1);
test1 = cell2mat(alltraces1_noT);
SEP.calculate_HeatMap(test1)
%%

test1 = cell2mat(alltraces2_T);
f = figure;
trace = test1;%./mean(test2(:,1:3),2);
allsum = sum(trace,2);
[~, ordx] = sort(allsum, 'descend');
ord_trace = trace(ordx,:);
goodvalsTm = movmean(ord_trace([15:end],:),4,2);
hm1 = heatmap(goodvalsTm);
hm1.GridVisible = 'off';
hm1.Colormap = hot(50);
hm1.ColorLimits = [0 2.9];
set(gcf,'Position',[2249         120         995         849]);

test2 = cell2mat(alltraces2_noT);
f = figure;
trace = test2;%./mean(test2(:,1:3),2);
allsum = sum(trace,2);
[~, ordx] = sort(allsum, 'descend');
ord_trace = trace(ordx,:);
goodvalsnTm = movmean(ord_trace,4,2);
hm2 = heatmap(goodvalsnTm);
hm2.GridVisible = 'off';
hm2.Colormap = hot(50);
hm2.ColorLimits = [0 2.9];
set(gcf,'Position',[2249         120         995         849]);

setval = 3;
hm1.ColorLimits = [0 setval];
hm2.ColorLimits = [0 setval];
%% cdf plots
figuredirname = 'E:\Matt Becker Data (For Review)\Figures\oldROIs';
allvalsnoT = cell2mat(alltraces2_noT);
allvalsT = cell2mat(alltraces2_T);
% allvalsT = allvalsT([15:end],:);




minval = 1.3;
id_noT = sum(allvalsnoT>minval,2)>2;
id_T = sum(allvalsT>minval,2)>2;


maxval = 2;

timewindow = 1:20;
newT = sort(allvalsT(id_T,timewindow),2,'descend');
newT(newT(:,1)>maxval,:) = [];
newnoT = sort(allvalsnoT(id_noT,timewindow),2,'descend');
newnoT(newnoT(:,1)>maxval,:) = [];
figure; hold on;
[F,X] = ecdf(newT(:,2));
plot(X,F,'LineWidth',3,'LineStyle','-.','Color',[0.2 0.2 0.2])
[G,Y] = ecdf(newnoT(:,2));
plot(Y,G,'LineWidth',3,'LineStyle','-','Color','k')
set(gca,'FontSize',16)
xlabel('Max Fold Increase (within first 25 frames)','FontSize',16)
ylabel('Cumulative Probability','FontSize',16);
xlim([0.5 maxval])
ylim([0 1])
[H,P] = kstest2(newT(:,2),newnoT(:,2))
% saveas(gcf,fullfile(figuredirname,'maxCDF'),'emf');
% saveas(gcf,fullfile(figuredirname,'maxCDF'),'fig');
% close(gcf);

%---

newT = mean(allvalsT(id_T,end-2:end),2);
newT(newT>maxval) = [];
newnoT = mean(allvalsnoT(id_noT,end-2:end),2);
newnoT(newnoT>maxval) = [];
figure; hold on;
[F,X] = ecdf(newT);
plot(X,F,'LineWidth',3,'LineStyle','-.','Color',[0.2 0.2 0.2])
[G,Y] = ecdf(newnoT);
plot(Y,G,'LineWidth',3,'LineStyle','-','Color','k')
set(gca,'FontSize',16)
xlim([0.5 maxval]); ylim([0 1]);
xlabel('F/F0 Final','FontSize',16)
ylabel('Cumulative Probability','FontSize',16);
[H,P] = kstest2(newT,newnoT)

% saveas(gcf,fullfile(figuredirname,'endCDF'),'emf');
% saveas(gcf,fullfile(figuredirname,'endCDF'),'fig');
% close(gcf);
%---

newT = sort(allvalsT(id_T,timewindow),2,'descend');
% newT(newT(:,1)>4,:) = [];
% newTend = allvalsT(id_T,end),2,'descend');
% newT(newT(:,1)>4,:) = [];
newnoT = sort(allvalsnoT(id_noT,timewindow),2,'descend');
% newnoT(newnoT(:,1)>4,:) = [];
figure; hold on;
% [F,X] = ecdf(newT(:,2)-allvalsT(id_T,end));
Tval =  allvalsT(id_T,end-1)-newT(:,2);
Tval(Tval<-2) = [];
[F,X] = ecdf(Tval);
plot(X,F,'LineWidth',3,'LineStyle','-.','Color',[0.2 0.2 0.2])
noTval = allvalsnoT(id_noT,end-1)-newnoT(:,2);
noTval(noTval<-2) = [];
[G,Y] = ecdf(noTval);
plot(Y,G,'LineWidth',3,'LineStyle','-','Color','k')
set(gca,'FontSize',16)
xlabel('Post Peak Decrease','FontSize',16)
ylabel('Cumulative Probability','FontSize',16);
% xlim([-1 1])
ylim([0 1])
[H,P] = kstest2(noTval,Tval)

% saveas(gcf,fullfile(figuredirname,'postpeakCDF'),'emf');
% saveas(gcf,fullfile(figuredirname,'postpeakCDF'),'fig');
% close(gcf);
% 

