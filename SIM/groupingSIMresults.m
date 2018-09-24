topdir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\080317\500 nM Abeta_PSD95_gephyrin';
ids = [3 1 2];
channelorderingstr = {'chABeta','PSD95i','Gephyrin'}; % channel abeta, channel 1, channel2
savedir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files';
if ~exist(savedir)
    mkdir(savedir);
end
% select the directory to choose files from and the corresponding
% parameters for it
% pick the files and then make and save SIM files
files = uipickfiles('FilterSpec',topdir,'Prompt',['Pick Files for ' channelorderingstr{1} ,...
    '-' channelorderingstr{2} '-' channelorderingstr{3}]);
bins = 0:31;
results1 = zeros(size(bins,2),numel(files));
results2 = zeros(size(bins,2),numel(files));
results_cell1 = cell(numel(files));
results_cell2 = cell(numel(files));


for ff = 1:numel(files)
    clear lb mask1_distance mask2_distance image cellmask
    disp(['Analyzing File: ' files]);
    s  = SIM();
    s.channelordering = ids;
    s.channelorderingstr = channelorderingstr;
    s.loadNDfile(files{ff});
    s.make_masks();
    s.make_distancemasks();
    lb = label(s.abeta.mask,1);
    mask1_distance = s.ch1.distance_mask;
    mask2_distance = s.ch2.distance_mask;
    image = s.abeta.image;
    cellmask = s.cellmask;
    results_cell1{ff} = s.calculateRadialDensity(mask1_distance,image,cellmask,bins);
    results_cell2{ff} = s.calculateRadialDensity(mask2_distance,image,cellmask,bins);
    results1(:,ff) = results_cell1{ff}.radial_density;
    results2(:,ff) = results_cell2{ff}.radial_density;
end
%%
raddist1 = mean(results1,2);
raddist2 = mean(results2,2);

raddist1_st = std(results1')./(1.96.*sqrt(size(results1,2)));
raddist2_st = std(results2')./(1.96.*sqrt(size(results2,2)));

figure;
errorbar(bins,raddist1,raddist1_st);
hold on;
errorbar(bins,raddist2,raddist2_st);
%%
topdir = 'G:\SIM\180704 Coverslip3';
ids = [1 2 3];
channelorderingstr = {'Gephyrin-GFP','anti-GFP','anti-Gephyrin'}; % channel abeta, channel 1, channel2
savedir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files';
if ~exist(savedir)
    mkdir(savedir);
end
% select the directory to choose files from and the corresponding
% parameters for it
% pick the files and then make and save SIM files
files = uipickfiles('FilterSpec',topdir,'Prompt',['Pick Files for ' channelorderingstr{1} ,...
    '-' channelorderingstr{2} '-' channelorderingstr{3}]);
bins = 0:31;
results1 = zeros(size(bins,2),numel(files));
results2 = zeros(size(bins,2),numel(files));
results_cell1 = cell(numel(files));
results_cell2 = cell(numel(files));


for ff = 1:numel(files)
    clear lb mask1_distance mask2_distance image cellmask
    disp(['Analyzing File: ' files]);
    s  = SIM();
    s.channelordering = ids;
    s.channelorderingstr = channelorderingstr;
    s.loadNDfile(files{ff});
    s.make_masks();
    s.make_distancemasks();
    lb = label(s.abeta.mask,1);
    mask1_distance = s.ch1.distance_mask;
    mask2_distance = s.ch2.distance_mask;
    image = s.abeta.image;
    cellmask = s.cellmask;
    results_cell1{ff} = s.calculateRadialDensity(mask1_distance,image,cellmask,bins);
    results_cell2{ff} = s.calculateRadialDensity(mask2_distance,image,cellmask,bins);
    results1(:,ff) = results_cell1{ff}.radial_density;
    results2(:,ff) = results_cell2{ff}.radial_density;
end
raddist1 = mean(results1,2);
raddist2 = mean(results2,2);

raddist1_st = std(results1')./(1.96.*sqrt(size(results1,2)));
raddist2_st = std(results2')./(1.96.*sqrt(size(results2,2)));

figure;
errorbar(bins,raddist1,raddist1_st);
hold on;
errorbar(bins,raddist2,raddist2_st);
legend({channelorderingstr{2},channelorderingstr{3}});