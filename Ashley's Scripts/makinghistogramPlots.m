%% --- 
% ---------------update these values!!!---------------
load('Z:\Ashley\For Sam\050118 NL1 insertion\cell3_AshleyFile\MinPaths.mat');
% datafilename = 'Z:\Ashley\For Sam\051318 NL1 insertion\7_merge_AshleyFile';
% load(fullfile(datafilename,'MinPaths.mat'));
px2um = 0.2857;
% glua1 info
baseline_frames = (1:6)-1;
baseline_minperframe = 2;
postrelease_minperframe = 2;
% nrl1 info
% baseline_frames = (1:(6-1));
% baseline_minperframe = 2;
% ------------- clean up the data ------------------
% first copy the distMat to clean things up and convert
distVtimemat = distMat;
% for now get rid of frame = 0 and frame = end
distVtimemat(distVtimemat(:,1) == 0,:) = [];
distVtimemat(distVtimemat(:,1) == max(distVtimemat(:,1)),:) = [];
distVtimemat(isnan(distVtimemat(:,1)),:) = [];
%convert to microns
distVtimemat(:,2) = distVtimemat(:,2) .* px2um;
%-- convert baseline and postrelease to minutes
baselineids = ismember(distVtimemat(:,1),baseline_frames);
postreleaseids = distVtimemat(:,1)>= (baseline_frames(end)+1);
distVtimemat(baselineids,1) = distVtimemat(baselineids,1) .* baseline_minperframe;
distVtimemat(postreleaseids,1) = distVtimemat(postreleaseids,1) .* postrelease_minperframe;
max(distVtimemat(:,2))
%% 
%  now you want distVtimemat to represent your results

%% ------------- now we can plot some things -----------------
% break cargo appearances up into times
distances = [0:20:100];
bins = zeros(numel(distances),1);
f = figure; cax = gca;
for tt = 1:numel(distances)
    if tt < numel(distances)
        binvals = distVtimemat(distVtimemat(:,2)>= distances(tt) & distVtimemat(:,2)< distances(tt+1),:);
        titlestr = (['Events Occuring between '...
            num2str(distances(tt)) ' - ' num2str(distances(tt+1)) ' \mum of Cell Soma']);
  
    else
        binvals = distVtimemat(distVtimemat(:,2)>= distances(tt),:);
        titlestr = (['Events Occuring at greater than '...
            num2str(distances(tt)) ' \mum of Cell Soma']);
    end
    subplot(numel(distances),1,tt)
    histo = histogram(binvals(:,1),10);
    title(titlestr,'FontSize',12)
    xlabel('Time (min)','FontSize',10);
    ylabel('# Cargo Events','FontSize',10);
end
% suptitle(['Cargo Delivery For File: ' datafilename])

