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
results3 = zeros(size(bins,2),numel(files));
results4 = zeros(size(bins,2),numel(files));
% 
results_cell1 = cell(numel(files));
results_cell2 = cell(numel(files));
results_cell3 = cell(numel(files));
results_cell4 = cell(numel(files));

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
    results_cell3{ff} = s.calculateRadialDensity(mask1_distance,s.ch1.image,cellmask,bins);
    results_cell4{ff} = s.calculateRadialDensity(mask2_distance,s.ch2.image,cellmask,bins);
    results1(:,ff) = results_cell1{ff}.radial_density;
    results2(:,ff) = results_cell2{ff}.radial_density;
    results3(:,ff) = results_cell3{ff}.radial_density;
    results4(:,ff) = results_cell4{ff}.radial_density;
end
save(fullfile('C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\results_radialdensity\psd95ibGephResults'),...
    'results_cell1','results_cell2','results_cell4','results_cell3','results1','results2','results3','results4');


% bar graph
dat4bar1 = results1(1:5,:);
dat4bar2 = results2(1:5,:);

mean1 = mean(dat4bar1(:));
mean2 = mean(dat4bar2(:));
st1 = std(dat4bar1(:)./(1.96.*sqrt(numel(dat4bar1))));
st2 = std(dat4bar2(:)./(1.96.*sqrt(numel(dat4bar2))));
leg={'A\beta from PSD95ib','A\beta from Gephyrin'};
c = categorical([leg(1);leg(2)]);
figure; 
bar(c,[mean1,mean2]);
hold on;
errorbar(c,[mean1,mean2],[st1,st2],'.b')
ylim([1 1.5])
 ylabel('Normalized A\beta Intensity')
 set(gca,'FontSize',14)

%%
raddist1 = mean(results1,2);
raddist2 = mean(results2,2);
raddist3 = mean(results3,2);
raddist4 = mean(results4,2);


raddist1_st = std(results1')./(1.96.*sqrt(numel(results1)));
raddist2_st = std(results2')./(1.96.*sqrt(numel(results2)));
raddist3_st = std(results3')./(1.96.*sqrt(numel(results3)));
raddist4_st = std(results4')./(1.96.*sqrt(numel(results4)));

bins = bins.*s.XYpxsize;
colors = lines(10);
leg={'ABeta from PSD95ib','ABeta from Gephyrin','Gephyrin from PSD95ib','PSD95ib from Gephyrin'};

figure; hold on;
errorbar(bins,raddist1,raddist1_st,'Color',colors(1,:),'DisplayName','95% CI'); 
plot(bins,raddist1,'Color',colors(1,:),'LineWidth',2,'DisplayName',leg{1});

hold on;
errorbar(bins,raddist2,raddist2_st,'Color',colors(2,:),'DisplayName','95% CI');
plot(bins,raddist2,'Color',colors(2,:),'LineWidth',2,'DisplayName',leg{2});

hold on;
errorbar(bins,raddist3,raddist3_st,'Color',colors(3,:),'DisplayName','95% CI');
plot(bins,raddist3,'Color',colors(3,:),'LineWidth',2,'DisplayName',leg{3});
hold on;
errorbar(bins,raddist4,raddist4_st,'Color',colors(4,:),'DisplayName','95% CI');
plot(bins,raddist4,'Color',colors(4,:),'LineWidth',2,'DisplayName',leg{4});


legend;

xlabel('Distance (\mum)');
% ylabel('Normalized Intensity Density');
set(gca,'FontSize',14)
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