filename = uipickfiles('Prompt','Pick all the AshleyFiles to Calculate Distances (can pick more than 1)','FilterSpec','C:\Users\schwsama\Documents\Data\zapERtrap');
if ~iscell(filename)
    filesize = 1;
    filename = {filename};
else
    filesize = numel(filename);
end
wb = waitbar(0,'making heatmaps...');
for ii = 1:filesize
clear aa;
load(fullfile(filename{ii}));
h = aa.plot_cargo_minFrame(180);
saveas(h,fullfile([filename{ii}(1:end-4) '_timeHeatMap']),'fig');
saveas(h,fullfile([filename{ii}(1:end-4) '_timeHeatMap']),'png');
close(h);
waitbar(ii/filesize,wb);
end
close(wb);



