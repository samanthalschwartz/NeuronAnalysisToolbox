% calculating some things!!
filename = uipickfiles('Prompt','Pick all the AshleyFiles to Calculate Distances (can pick more than 1)','FilterSpec','C:\Users\schwsama\Documents\Data\zapERtrap');
if ~iscell(filename)
    filesize = 1;
    filename = {filename};
else
    filesize = numel(filename);
end
alldistances = cell(1,filesize);
wb = waitbar(0,'loading distances...');
for ii = 1:filesize
    clear aa;
load(fullfile(filename{ii}));
alldistances{ii} = aa.distancematrix;
waitbar(ii/filesize,wb);
end
close(wb);
% gluA1distances = alldistances;
% gluA1files = filename
NL1distances = alldistances;
NL1files = filename

%% 
currfilename = NL1files;
currdistances = NL1distances;
% currfilename = gluA1files;
% currdistances = gluA1distances;

savedir = 'C:\Users\schwsama\Documents\Data\zapERtrap\180814 PlotsIII';
if ~exist(savedir)
    mkdir(savedir)
end
max_dist = 150;
min_dist = 5;
min_dist2 = 30;
inc = 30;
% distarry = min_dist:inc:max_dist;
distarry = [min_dist,min_dist2:inc:max_dist];
timebins = 0:4:180;


% max_dist = 150;
% min_dist = 20;
% inc = 20;
% distarry = min_dist:inc:max_dist;
% timebins = 0:180;

for adid = 1:numel(currdistances)
f = figure; hold on;
totevents = sum(currdistances{adid}(:,6)>min_dist);

for distid = 1:numel(distarry)-1
min_i = distarry(distid);
max_i = distarry(distid+1);
idx = currdistances{adid}(:,6)>min_i & currdistances{adid}(:,6)<max_i;
forcompare = currdistances{adid}(idx,5); %time column
[N,edges] = histcounts(forcompare,timebins);
newy = N;
newy(newy==0) = NaN;
x  = edges(1:end-1);
plot(x,newy./totevents,'o','DisplayName',[num2str(min_i) '-' num2str(max_i) ' microns from soma']); hold on;
end
legend('Location','northwest');
% ymax = max(150,max(get(gca,'Ylim')));
xlim([timebins(1),timebins(end)]);
% ylim([0 ymax]);
[path,NAME,~] = fileparts(currfilename{adid});
nametitle = strrep(NAME,'_','-');
title(nametitle);
xlabel('Time after Release (min)');
ylabel('# of Events');
saveas(f,fullfile(savedir,[NAME '_numEventsPlot']),'fig');
saveas(f,fullfile(savedir,[NAME '_numEventsPlot']),'png');
close(f);
end

%% 
% savedir = 'C:\Users\schwsama\Documents\Data\zapERtrap\180814 Plots';
% if ~exist(savedir)
%     mkdir(savedir)
% end

% currdistance = NL1distances;
% titlestr = 'NL-1';
currdistance = gluA1distances(1:end);
titlestr = 'GluA-1';

max_dist = 150;
min_dist = 5;
min_dist2 = 40;
inc = 30;
% distarry = min_dist:inc:max_dist;
distarry = [min_dist,min_dist2:inc:max_dist];
timebins = 0:2:122;
distance_bins = cell(numel(distarry)-1,1);
display_Name = {};
for adid = 1:numel(currdistance)
    % f = figure; hold on;
    totevents = sum(currdistance{adid}(:,6)>min_dist);
    for distid = 1:numel(distarry)-1
        min_i = distarry(distid);
        max_i = distarry(distid+1);
        if adid == 1
           display_Name{distid} = [num2str(min_i) '-' num2str(max_i) ' (\mum from Soma)'];   
        end
        
        idx = currdistance{adid}(:,6)>min_i & currdistance{adid}(:,6)<max_i;
        forcompare = currdistance{adid}(idx,5); %time column
        [N,edges] = histcounts(forcompare,timebins);
%         totevents = sum(N);
        newy = N./totevents;
%         newy = N./totevents;
%         if distid == 1
%             normfactor = max(N(1:end-1));
%         end
%         newy = N./normfactor;
        newy(newy==0) = NaN;
        x = timebins(1:end-1);
%         x  = edges(1:end-1);
        % plot(x,newy,'o','DisplayName',[num2str(min_i) '-' num2str(max_i) ' microns from soma']); hold on;
        distance_bins{distid}(adid,:) = newy;
    end
end

figure
cols = lines(numel(distance_bins));
for ii =1:numel(distance_bins)
    frac_values = nanmean(distance_bins{ii},1);
%     frac_values_std = nanstd(distance_bins{ii},1);
    frac_values_std = nanstd(distance_bins{ii},1)./sqrt(size(distance_bins{ii},2));
%     plot(timebins(1:end-1),frac_values,'-','DisplayName',display_Name{ii}); hold on;
    errorbar(timebins(1:end-1),frac_values,frac_values_std,...
        '-','CapSize',0,'Color',cols(ii,:),'DisplayName',' ');hold on;
    plot(timebins(1:end-1),frac_values,'-','DisplayName',display_Name{ii},'LineWidth',2,'Color',cols(ii,:)); hold on;
    set(gca,'FontSize',16);
    set(gca,'Position',[0.1300 0.1105 0.7750 0.8145]);
    set(gca,'OuterPosition',[0 5.6748e-04 1 0.9994]);
    set(gcf,'Position',[2217 195 872 702])
%     legend('Location','northwest');
% % ymax = max(150,max(get(gca,'Ylim')));
% xlim([timebins(1),timebins(end)]);
ylim([0 .02]);
% [path,NAME,~] = fileparts(filename{adid});
% nametitle = strrep(NAME,'_','-');
% title(nametitle);
xlabel('Time after Release (min)');
ylabel('Fraction of Total Events');
title(titlestr);
% saveas(f,fullfile(savedir,[NAME '_numEventsPlot']),'fig');
% saveas(f,fullfile(savedir,[NAME '_numEventsPlot']),'png');
% close(f);
end
legend('Location','northwest');
%%
excelsavename = 'Z:\Sam\MJK_zapERtrap_for_sam\180816 Results\180816 GluA-1 Results';
for ii = 1:numel(distance_bins)
sheetname = strrep(display_Name{ii},'\mum', 'microns');
    xlswrite(excelsavename,distance_bins{ii},sheetname);
end
xlswrite(excelsavename,timebins(1:end-1),'TimeBins')
