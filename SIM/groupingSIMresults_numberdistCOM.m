
%%
topdir = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071917\Abeta_PSD95_bassoon';
ids = [3 1 2];
channelorderingstr = {'chABeta','PSD95ib','Gephyrin'}; % channel abeta, channel 1, channel2
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
    
    mask1_distance = s.ch1.distance_mask;
    mask2_distance = s.ch2.distance_mask;
    image = s.abeta.image;
    lbimage = s.abeta.labeled_mask;
    cellmask = s.cellmask;
    
    results_cell1{ff} = s.calculateNumberDensityCOM(mask1_distance,lbimage,image,cellmask,bins);
    results_cell2{ff} = s.calculateNumberDensityCOM(mask2_distance,lbimage,image,cellmask,bins);
    results1(:,ff) = (results_cell1{ff}.radialnumber./results_cell1{ff}.nummask1)./(results_cell1{ff}.totalnumber/results_cell1{ff}.nummask1);
    results2(:,ff) = (results_cell2{ff}.radialnumber./results_cell2{ff}.nummask1)./(results_cell2{ff}.totalnumber/results_cell2{ff}.nummask1);
end
raddist1 = mean(results1,2);
raddist2 = mean(results2,2);

raddist1_st = std(results1')./(1.96.*sqrt(size(results1,2)));
raddist2_st = std(results2')./(1.96.*sqrt(size(results2,2)));

figure;
leg={'ABeta from PSD95ib','ABeta from Gephyrin'};
colors = lines(10);
hold on;
errorbar(bins,raddist1,raddist1_st,'Color',colors(1,:),'DisplayName','95% CI');
plot(bins,raddist1,'Color',colors(1,:),'LineWidth',2,'DisplayName',leg{1});
hold on;
errorbar(bins,raddist2,raddist2_st,'Color',colors(2,:),'DisplayName','95% CI');
plot(bins,raddist2,'Color',colors(2,:),'LineWidth',2,'DisplayName',leg{2});


ylabel('Normalized #of AB/Gephyrin');
xlabel('Distance in pixels');
set(gca,'FontSize',14)
legend;

save(fullfile('C:\Users\KennedyLab\Documents\Hannah\SIM data\SIM_files\results_numberdensityCOM\psd95ibGephResults'),...
    'results_cell1','results_cell2','results1','results2')
%%
results_num1 = zeros(size(bins,2),size(results_cell1,2));
results_num2 = zeros(size(bins,2),size(results_cell1,2));
bins = 0:31;

for ff = 1:size(results_cell1,2)
results_num1(:,ff) = (results_cell1{ff}.radialnumber./results_cell1{ff}.volume)./(results_cell1{ff}.totalnumber./results_cell1{ff}.totalvolume);
results_num2(:,ff) = (results_cell2{ff}.radialnumber./results_cell2{ff}.volume)./(results_cell1{ff}.totalnumber./results_cell1{ff}.totalvolume);
end
raddist1 = mean(results_num1,2);
raddist2 = mean(results_num2,2);

raddist1_st = std(results_num1')./(1.96.*sqrt(size(results_num1,2)));
raddist2_st = std(results_num2')./(1.96.*sqrt(size(results_num2,2)));

figure;
leg={'ABeta from PSD95ib','ABeta from Gephyrin'};
colors = lines(10);
hold on;
errorbar(bins.*s.XYpxsize,raddist1,raddist1_st,'Color',colors(1,:),'DisplayName','95% CI');
plot(bins.*s.XYpxsize,raddist1,'Color',colors(1,:),'LineWidth',2,'DisplayName',leg{1});
hold on;
errorbar(bins.*s.XYpxsize,raddist2,raddist2_st,'Color',colors(2,:),'DisplayName','95% CI');
plot(bins.*s.XYpxsize,raddist2,'Color',colors(2,:),'LineWidth',2,'DisplayName',leg{2});


ylabel('Fold Increase in Number of A\beta Puncta');
xlabel('Distance in \mum');
set(gca,'FontSize',14)
legend;

