%% for just quickly plotting the results of one file
figure; hold on;
index = 7;
for mm = 1:size(allMs{index}.intensity,1)
    windowsz = 10;
%     vals = movmean(diff(allMs{1}.intensity(mm,:)),windowsz);
    vals = allMs{index}.intensity(mm,:);
    plot(vals);
end
legend({'1','2','3'})

%% this is for calculating and grouping the allMs files - change the load and save filepaths
h = waitbar(0,'Making some mean intensity calculations...');
clear files;
load('Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localDendrites.mat');
allMs = {};
for ff = 1:numel(files)
    clear aa;
    load(files{ff});
   M = aa.plotDensityperTime; 
   waitbar(ff/numel(files))
   allMs{ff} =M;
end
close(h);

savename = 'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localDendrites_allMs.mat';
save(savename,'allMs')

%% this is for grouping all the data from allMs cell array
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat';
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat';
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\all_GluA1_global_allMs.mat';
filename = 'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat'
% 
[FILEPATH,name,EXT]  = fileparts(filename);
imtitlename = strrep(name,'_','-');
load(filename)% distance 0-100 uM from soma
% cols = [autumn(numel(allMs)/3); summer(numel(allMs)/3); spring(numel(allMs)/3)];

cols = hsv(numel(allMs));
f = figure; fa = gca;

id = 1;
subplot(3,1,id); hold on;
sumvals = cellmap(@(x) sum(x.intensity(x.intensity>0)),allMs);
output = cellmap(@(x,y) x.intensity(id,:)/y,allMs,sumvals);
for oo = 1:numel(output)
    output{oo}(output{oo}<0) = 0;
   plot(output{oo},'Color',cols(oo,:));    
%    ylim([0 20e6]);
%    ylim([0 2e1]);
end
title('Mean Intensity <100 \mum from Soma')

id = 2;
subplot(3,1,2); hold on;
sumvals = cellmap(@(x) sum(x.intensity(x.intensity>0)),allMs);
output = cellmap(@(x,y) x.intensity(id,:)/y,allMs,sumvals);
for oo = 1:numel(output)
     output{oo}(output{oo}<0) = 0;
   plot(output{oo},'Color',cols(oo,:));  
%    ylim([0 20e5]);
end
title('Mean Intensity 100-200 \mum from Soma')

id = 3;
subplot(3,1,3); hold on;
sumvals = cellmap(@(x) sum(x.intensity(x.intensity>0)),allMs);
output = cellmap(@(x,y) x.intensity(id,:)/y,allMs,sumvals);
for oo = 1:numel(output)
     output{oo}(output{oo}<0) = 0;
   plot(output{oo},'Color',cols(oo,:));  
%    ylim([0 20e5]);
end
title('Mean Intensity >300 \mum from Soma')

l = legend(string(1:numel(allMs)),'Location','northwest');

hout=suptitle(['Mean Delivery Density for Cargo: ' imtitlename ]);

set(f,'Position',[59   301   904   688]);
set(l,'Position',[0.0195    0.3030    0.0697    0.6374]);
%% this is for grouping and plotting the mean intensity with the standard error
filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat',...
    'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat',...
    'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat',...
    'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\all_GluA1_global_allMs.mat';
titlenames = {'Mean Intensity <100 \mum from Soma','Mean Intensity 100-200 \mum from Soma',...
    'Mean Intensity >300 \mum from Soma'};


for ff = 1:numel(filename)
    f = figure; fa = gca;
[FILEPATH,name,EXT]  = fileparts(filename{ff});
imtitlename = strrep(name,'_','-');
load(filename{ff})% distance 0-100 uM from soma
% cols = [autumn(numel(allMs)/3); summer(numel(allMs)/3); spring(numel(allMs)/3)];

cols = hsv(numel(allMs));

for id = 1:3
subplot(3,1,id); hold on;
sumvals = cellmap(@(x) sum(x.intensity(x.intensity>0)),allMs);
output = cellmap(@(x,y) x.intensity(id,:)/y,allMs,sumvals);
% output = cellmap(@(x) x.intensity(id,:),allMs);
mxtime = max(cellfun(@(x) size(x,2), output))
resultsarr = nan(numel(output),mxtime);
for oo = 1:numel(output)
    currvals = output{oo};
    resultsarr(oo,1:size(currvals,2)) = currvals;
end
means = nanmean(resultsarr,1);
countnans = ~isnan(resultsarr);
cnts = sum(countnans,1);
errs = nanstd(resultsarr,1)./sqrt(cnts);
errorbar(means,errs)   
title(titlenames{id})
end
hout=suptitle(['Mean Delivery Density for Cargo: ' imtitlename ]);

end
%%
% l = legend(string(1:numel(allMs)),'Location','northwest');


set(f,'Position',[59   301   904   688]);
set(l,'Position',[0.0195    0.3030    0.0697    0.6374]);

%% this is for plotting intensities - normalizing by the last couple max frames
%% this is for grouping and plotting the mean intensity with the standard error
filename = {'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_global_allMs.mat',...
    'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_TfR_Global_allMs.mat',...
    'Z:\Sam\MJK_zapERtrap_for_sam\Ashley062018_all_GluA1_Global_allMs.mat',...
    'Z:\Sam\MJK_zapERtrap_for_sam\all_NL1_localSoma_allMs'};
% filename = 'Z:\Sam\MJK_zapERtrap_for_sam\all_GluA1_global_allMs.mat';
titlenames = {'Mean Intensity <100 \mum from Soma','Mean Intensity 100-200 \mum from Soma',...
    'Mean Intensity >300 \mum from Soma'};


for ff = 1:numel(filename)
    f = figure; fa = gca;
[FILEPATH,name,EXT]  = fileparts(filename{ff});
imtitlename = strrep(name,'_','-');
load(filename{ff})% distance 0-100 uM from soma
% cols = [autumn(numel(allMs)/3); summer(numel(allMs)/3); spring(numel(allMs)/3)];

cols = hsv(numel(allMs));

for id = 1:3
subplot(3,1,id); hold on;
sumvals = cellmap(@(x) sum(x.intensity(x.intensity>0)),allMs);
output = cellmap(@(x,y) x.intensity(id,:)/y,allMs,sumvals);
% output = cellmap(@(x) x.intensity(id,:),allMs);
mxtime = max(cellfun(@(x) size(x,2), output))
resultsarr = nan(numel(output),mxtime);
for oo = 1:numel(output)
    currvals = output{oo};
    resultsarr(oo,1:size(currvals,2)) = currvals;
end
means = nanmean(resultsarr,1);
countnans = ~isnan(resultsarr);
cnts = sum(countnans,1);
errs = nanstd(resultsarr,1)./sqrt(cnts);
errorbar(means,errs)   
title(titlenames{id})
end
hout=suptitle(['Mean Delivery Density for Cargo: ' imtitlename ]);

end
