%% loop through and save line plots
close all; clear all;
datafolder = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\results\moving masks';
files = dir(fullfile(datafolder,'*.mat'));
minval = 1.2;

for ff = 1:numel(files)
    clear in;
    filename = files(ff).name(1:end-4);
    in = load(fullfile(datafolder,filename));
    % -- can start from here after loading saved file
    nr_ids = sum(in.results.normraw>minval,1);
    nm_ids = sum(in.results.normmean>minval,1);
    % figure; plot(in.results.normraw(:,nr_ids>0)); hold on; plot(mean(in.results.normraw(:,nr_ids>0),2),'--','LineWidth',2,'Color','k');
    f1 = figure; plot(in.results.normmean(:,nm_ids>0));hold on; plot(mean(in.results.normmean(:,nm_ids>0),2),'--','LineWidth',2,'Color','k');
    ylim([0.5 3]);
    xlabel('Frame Number');
    ylabel('Fold Increase in GluA1 Fluorescence at Homer1 Puncta');
    title([filename ': Only Includes Puncta with  > ' num2str(minval) ' Fold Increase']);
    saveas(f1, fullfile(datafolder,[filename '_AllTraces']),'fig');
    saveas(gcf, fullfile(datafolder,[filename '_AllTraces']),'png');
    close all;
    figure;  plot(mean(in.results.normmean(:,nm_ids>0),2),'--','LineWidth',2,'Color','k');
    ylim([0.9 1.35]);
    xlabel('Frame Number');
    ylabel('Fold Increase in GluA1 Fluorescence at Homer1 Puncta');
    title([filename ': Only Includes Puncta with  > ' num2str(minval) ' Fold Increase']);
    saveas(gcf, fullfile(datafolder,[filename '_MeanTrace']),'fig');
    saveas(gcf, fullfile(datafolder,[filename '_MeanTrace']),'png');
    close all;
end
%%
%% now find 'good' rois and plot ontop of 488 image the perimeter of the mask
gdtest = sum(in.results.normmean>minval,1);
gdids = goodid(gdtest>0);
%% 
% pickednorms = norms(gdtest>0,:);
% figure; plot(pickednorms');
% figure; plot(mean(pickednorms)');
[f,x] = ecdf(max(in.results.normmean,[],1));
figure; plot(x,f);
title('Cell 6')
%%
datafolder = 'G:\Sam\Data\180516_GluA1mCh-HomerGFP\results';
savename = 'CD Plots';
files = dir(fullfile(datafolder,'*.mat'));
imout = figure(); hold on;
for ff = 1:numel(files)
    clear in;
    in = load(fullfile(datafolder,files(ff).name));
    [f,x] = ecdf(max(in.results.normmean,[],1));
    plot(x,f,'DisplayName',files(ff).name,'LineWidth',2);
end
xlim([1 2]);
xlabel('Max Fold Increase in GluA1 Fluroscence at Homer1 Puncta');
ylabel(['Cumulative Probability']);
legend();
saveas(gcf,fullfile(datafolder,savename),'png')
saveas(gcf,fullfile(datafolder,savename),'fig')
