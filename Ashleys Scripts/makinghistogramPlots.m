%% --- 
% ---------------update these values!!!---------------
close all; clear all
%%
% datapath ='G:\Sam\Data\Ashley\050118 NL1 insertion';
% datafiles = {'cell1_AshleyFile','cell3_AshleyFile','cell4_AshleyFile','cell5_AshleyFile'};
% savename = '050118-NL1-insertion';
% baseline_frames = 1:(7-1);


% datapath = 'G:\Sam\Data\Ashley\DHFR-GFP-GluA1\051718\merges\slip1';
% datafiles = {'1_merge_AshleyFile','2_merge_AshleyFile','3_merge_AshleyFile','4_merge_AshleyFile','5_merge_AshleyFile'};
% savename = 'DHFR-GFP-GluA1-051718-slip1';
% baseline_frames = 1:(11-1);
% 
% % 
% datapath = 'G:\Sam\Data\Ashley\DHFR-GFP-GluA1\051718\merges\slip2';
% datafiles = {'1_merge_AshleyFile','2_merge_AshleyFile','3_merge_AshleyFile',...
% '4_merge_AshleyFile','5_merge_AshleyFile'};
% savename = 'DHFR-GFP-GluA1-051718-slip2';
% baseline_frames = 1:(11-1);

% 
% datapath = 'G:\Sam\Data\Ashley\051318 NL1 insertion';
% datafiles = {'1_merge_AshleyFile','3_merge_AshleyFile',...
% '4_merge_AshleyFile','5_merge_AshleyFile','6_merge_AshleyFile',...
% '7_merge_AshleyFile','8_merge_AshleyFile'};
% savename = '051318-NL1-insertion';
% baseline_frames = 1:(11-1);
% %---

% datapath = 'G:\Sam\Data\Ashley\GluA1release_REs\060118\merges\slip1';
% datafiles = {'1_merge_AshleyFile','2_merge_AshleyFile','3_merge_AshleyFile',...
% '4_merge_AshleyFile','5_merge_AshleyFile','6_merge_AshleyFile'};
% savename = 'GluA1release_REs-060118-slip1';
% baseline_frames = 1:(11-1);

datapath = 'G:\Sam\Data\Ashley\GluA1release_REs\060118\merges\slip2';
datafiles = {'1_merge_AshleyFile','2_merge_AshleyFile'};
savename = 'GluA1release_REs-060118-slip2';
baseline_frames = 1:(11-1);


%%
histogram_savedir = 'G:\Sam\Data\Ashley\180607 Results\histogram distances';
savedir = 'G:\Sam\Data\Ashley\180607 Results';
px2um = 0.2857;
% glua1 info
baseline_minperframe = 1;
postrelease_minperframe = 2;
distances = [0:30:120];

savedhistograms = cell(numel(datafiles),numel(distances));
for ff = 1:numel(datafiles)
curr = load(fullfile(datapath,datafiles{ff},'MinPaths.mat'));
%% to get just post release time course
% %  now you want distVtimemat to represent your results
% % postrelease_mat = distVtimemat(distVtimemat(:,1)>baseline_frames(end),:);

postreleaseids = curr.distMat(:,1)>= (baseline_frames(end)+1);
postrelease_mat = curr.distMat(postreleaseids,:);
% for now get rid of frame = 0 and frame = end
postrelease_mat(postrelease_mat(:,1) == max(postrelease_mat(:,1)),:) = [];
postrelease_mat(isnan(postrelease_mat(:,1)),:) = [];
%convert to microns
postrelease_mat(:,2) = postrelease_mat(:,2) .* px2um;
%convert to mins
postrelease_mat(:,1) = postrelease_mat(:,1) .* postrelease_minperframe;
%% loop at save histogram objects for each time group
for tt = 1:numel(distances)
    if tt < numel(distances)
        binvals = postrelease_mat(postrelease_mat(:,2)>= distances(tt) & postrelease_mat(:,2)< distances(tt+1),:);  
        info.titlestr = (['Events Occuring between '...
            num2str(distances(tt)) ' - ' num2str(distances(tt+1)) ' \mum of Cell Soma']);

    else
        binvals = postrelease_mat(postrelease_mat(:,2)>= distances(tt),:);
        info.titlestr = (['Events Occuring at greater than '...
            num2str(distances(tt)) ' \mum of Cell Soma']);
    end
    histo = histogram(binvals(:,1),'BinEdges',distances); 
    info.xvals = diff(histo.BinEdges)/2 + histo.BinEdges(1:end-1);
    info.yvals = histo.Values./size(postrelease_mat,1); % normalize to total events
    savedhistograms{ff,tt} = info;    
end
end
save(fullfile(histogram_savedir,['HistogramResults-' savename]),'savedhistograms');
%% make the plot
f = figure; hold on;
for tt = 1:numel(distances)
    subplot(numel(distances),1,tt); hold on;
    for ff = 1:numel(datafiles)
        info = savedhistograms{ff,tt};
        plot(info.xvals,info.yvals','DisplayName',strrep(datafiles{ff},'_','-'));
    end
    if tt ==1
    l = legend(); set(l,'Position', [0.1451    0.7968    0.1535    0.1050])
    end
    title(info.titlestr,'FontSize',12)
    xlabel('Time (min)','FontSize',10);
    ylabel('% Total Cargo Events','FontSize',10);    
    clear info;
end
suptitle(['Cargo Delivery For File: ' datapath])
set(gcf,'Position',[1          41         912        1040]);
saveas(f,fullfile(datapath,[savename '_DistanceHistograms']),'fig');
saveas(f,fullfile(datapath,[savename '_DistanceHistograms']),'png');
saveas(f,fullfile(savedir,[savename '_DistanceHistograms']),'png');
%%
% suptitle(['Cargo Delivery For File: ' datafilename])

% %% ------------- now we can plot some things -----------------
% % break cargo appearances up into times
% 
% bins = zeros(numel(distances),1);
% f = figure; cax = gca;
% for tt = 1:numel(distances)
%     if tt < numel(distances)
%         binvals = distVtimemat(distVtimemat(:,2)>= distances(tt) & distVtimemat(:,2)< distances(tt+1),:);
%         titlestr = (['Events Occuring between '...
%             num2str(distances(tt)) ' - ' num2str(distances(tt+1)) ' \mum of Cell Soma']);
%     else
%         binvals = distVtimemat(distVtimemat(:,2)>= distances(tt),:);
%         titlestr = (['Events Occuring at greater than '...
%             num2str(distances(tt)) ' \mum of Cell Soma']);
%     end
%     subplot(numel(distances),1,tt)
%     histo = histogram(binvals(:,1),10);
%     title(titlestr,'FontSize',12)
%     xlabel('Time (min)','FontSize',10);
%     ylabel('# Cargo Events','FontSize',10);
% end
% suptitle(['Cargo Delivery For File: ' datafilename])

%%
% %% to get entire time course ------------- clean up the data ------------------
% % first copy the distMat to clean things up and convert
% distVtimemat = distMat;
% % for now get rid of frame = 0 and frame = end
% distVtimemat(distVtimemat(:,1) == 0,:) = [];
% distVtimemat(distVtimemat(:,1) == max(distVtimemat(:,1)),:) = [];
% distVtimemat(isnan(distVtimemat(:,1)),:) = [];
% %convert to microns
% distVtimemat(:,2) = distVtimemat(:,2) .* px2um;
% %-- convert baseline and postrelease to minutes
% baselineids = ismember(distVtimemat(:,1),baseline_frames);
% postreleaseids = distVtimemat(:,1)>= (baseline_frames(end)+1);
% distVtimemat(baselineids,1) = distVtimemat(baselineids,1) .* baseline_minperframe;
% distVtimemat(postreleaseids,1) = distVtimemat(postreleaseids,1) .* postrelease_minperframe;
% max(distVtimemat(:,2))