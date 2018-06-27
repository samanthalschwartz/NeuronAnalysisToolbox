listdir = 'G:\Sam\Data\MJK_zapERtrap_for_sam';
% filelistname = 'all_NL1_global.mat';
filelistname = 'Ashley062018_all_TfR_Global.mat';
zdrivepath = 'Z:\Sam\MJK_zapERtrap_for_sam';
savedir = fullfile(listdir,'180626 Results and Info');

thefiles = load(fullfile(listdir,filelistname));
allMs = cell(1,numel(thefiles.files));
wb = waitbar(0,'Calculating intensity density');
for ff = 1:numel(thefiles.files)  
    waitbar(ff/numel(thefiles.files),wb);
    clear aa;
    filepath = strrep(thefiles.files{ff},zdrivepath,listdir);
    load(filepath);
    M = aa.plotDensityperTime();
    allMs{ff} = M;
end
close(wb)

save(fullfile(savedir,[filelistname(1:end-4) '_allMs']),'allMs');
%% group the M.areanormintensity and plot

distval = 3;
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
% figure; plot(allintensities'); hold on;
figure; plot(allintensities','Color',[0 0.7 0.7]); hold on;
plot(nanmean(allintensities,1),'Color','k','LineWidth',2)
% legend(thefiles.files(range))