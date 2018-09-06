% calling SEP analysis from SEP class
% filepath = 'E:\Matt Becker Data (For Review)\SEPGlua1_mch';
filepath = 'E:\Matt Becker Data (For Review)\SEPGlua1_mch2aTeNT';
savedir = fullfile(filepath,'SEPfiles');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
files = uipickfiles('Prompt','Pick Files','FilterSpec',filepath);


for ff= 1:numel(files)
filename = files{ff};
display(['Analyzing File: ' filename]);
[path,name,ext] = fileparts(filename);
% load and calculate things
sp = SEP(); sp.filepath = path; sp.filename = [name, ext];
disp('Loading Images...');
sp.loadimages;
disp('Masking Images...');
sp.make_mask_cellfill;
sp.make_mask_sep;
sp.make_mask_sep_fixed;
disp('Calculating individual mask intensities...');
sp.calculate_sepintensities;
disp('Calculating overall mask intensities...');
sp.calculateSpineShaftIntensities;
sp.loadIJROIs1;
sp.loadIJROIs2;
% save and plot things
disp('Saving SEP File...');
sp.saveSEP(savepath);
h = sp.viewSEPMaskFixedwOldRois();
% saveas(h,fullfile(savepath,[sp.filename(1:end-4) '_allROIs']),'fig');
saveas(h,fullfile(savedir,[sp.filename(1:end-4) '_allROIs']),'png');
close(h);
f = sp.calculate_sepIHeatMap;
set(f,'Position', [2146,49,1193,935]);
% saveas(f,fullfile(savepath,[sp.filename(1:end-4) '_heatmap']),'fig');
saveas(f,fullfile(savedir,[sp.filename(1:end-4) '_heatmap']),'png');
close(f);
end


