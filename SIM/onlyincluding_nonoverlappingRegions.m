% only include regions with little overlap between pre and post  --- not
% really overlapping in z
close all; clear all;
filepath = uipickfiles('Prompt','Pick Files to Plot','FilterSpec','C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\AB_Controls_040819');
prepostdis_list_all = cell(1,numel(filepath));
numabeta_aroundsynapse_all = cell(1,numel(filepath));
plotvals_all = cell(1,numel(filepath));

for ff = 1:numel(filepath)
    clear obj;
    load(filepath{ff});
    % select the good ids
    okids = nan(numel(obj.results.selectedROIs),1);
    wb = waitbar(0);
    for ii = 1:numel(obj.results.selectedROIs)
        okids(ii) = obj.pre_postOverlap(obj.results.selectedROIs{ii},40);
        waitbar(ii/numel(obj.results.selectedROIs),wb)
    end
    close(wb)
    if sum(okids)==0 || isnan(range(okids))
        continue;
    end
    okids = logical(okids);
    prepostdis_list_all{ff} =  obj.results.prepostdis_list(okids');
    numabeta_aroundsynapse_all{ff} = obj.results.numabeta_aroundsynapse(okids);
    plotvals_all{ff} = obj.results.plotvals(:,okids);
    disp(['File Finished: ' num2str(ff) ]);
end
%%
numabeta = [];
all_dists = [];
plotsvals = [];
figure;hold on;

cols = lines(numel(plotvals_all));
for ii = 1:numel(plotvals_all)
    if isempty(plotvals_all{ii})
        continue;
    end
   all_dists = [all_dists; prepostdis_list_all{ii}'];
   numabeta = [numabeta; numabeta_aroundsynapse_all{ii}'];
   plotsvals = cat(2,plotsvals,plotvals_all{ii});
   scatter(plotvals_all{ii}(1,:),plotvals_all{ii}(2,:),'MarkerFaceColor',cols(ii,:),'MarkerEdgeColor',cols(ii,:))
end
figure; scatter(plotsvals(1,:).*obj.XYpxsize,plotsvals(2,:).*obj.XYpxsize,'*'); pbaspect([1 1 1])
ylabel('Distance from Center of POST-PRE')
xlabel('Distance Perpendicular to Synapse')
xlim([0 20]); ylim([-20 20])
hold on; plot([0:40],repmat(0,1,41),'w')

plotsvals(plotsvals(:,1)==0,:) = [];
[gca,N] = scatter2heatmap(plotsvals(1,:).*obj.XYpxsize,plotsvals(2,:).*obj.XYpxsize,15);
% colormap('hot')

xlim([0 1]); ylim([-1 1]); pbaspect([1 1 1])
hold on; plot([0:40],repmat(0,1,41),'w','LineWidth',2)
ylabel('Distance from Center of Pre-Post')
xlabel('Distance Perpendicular to Synapse')

%%
joinchannels('rgb',obj.abeta.mask,obj.ch1.mask,obj.ch2.mask)
for pp = 1:numel(obj.results.selectedROIs(okids))
   vertices = obj.results.selectedROIs{pp};
   p = patch('Vertices',vertices,'EdgeColor',[1 0 0],'Faces',1:size(vertices,1),'FaceAlpha',0);
end


