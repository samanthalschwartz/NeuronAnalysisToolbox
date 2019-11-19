clear all
close all

pc = PairCorr();
pc.Results = 'results';

[x, y] = textread('../data/9021_5.txt',  '%*u %u %u %*u', 'headerlines', 1);
XY_5 =  [x, y];
[x, y] = textread('../data/9021_10.txt', '%*u %u %u %*u', 'headerlines', 1);
XY_10 = [x, y];
pc.ROI = [0, 0, 7400, 6000];

results_pcc  = pc.pair_correlation(2.7559, '9021', XY_5, XY_10)
results_pac1 = pc.pair_correlation(2.7559, '9021_5',  XY_5)
results_pac2 = pc.pair_correlation(2.7559, '9021_10', XY_10)

results_Vpcc = ...
   pc.pair_correlation_Veatch(XY_5, XY_10, 2.7559, '9021', 'cross');
results_Vpac1 = ...
   pc.pair_correlation_Veatch(XY_5,  [], 2.7559, '9021_5', 'auto');
results_Vpac2 = ...
   pc.pair_correlation_Veatch(XY_10, [], 2.7559, '9021_10','auto');

load('../../analysis/methods/data/a_3variables.mat');
XY = regs(1).points;
ROI = regs(1).ROI;
pc.ROI = ROI;

results_pac  = pc.pair_correlation(10, 'a_3variables', XY)

results_Vpac = pc.pair_correlation_Veatch(XY, [], 10, 'a_3variables', 'auto');

pc.ROI = [0, 0, 2625, 2625];
[x1, y1] = ...
   textread('../data/B20_A431_2min_6.txt',  '%*n %n %n', 'headerlines', 1);
[x2, y2] = ...
   textread('../data/B20_A431_2min_12.txt', '%*n %n %n', 'headerlines', 1);
XY1 = [x1, y1];
XY2 = [x2, y2];
clear x1 y1 x2 y2

results_pcc  = pc.pair_correlation(10, 'B20_A431_2min', XY1, XY2)
results_pac1 = pc.pair_correlation(10, 'B20_A431_2min_6',  XY1)
results_pac2 = pc.pair_correlation(10, 'B20_A431_2min_12', XY2)

results_Vpcc = ...
   pc.pair_correlation_Veatch(XY1, XY2, 10, 'B20_A431_2min', 'cross');
results_Vpac1 = ...
   pc.pair_correlation_Veatch(XY1, [], 10, 'B20_A431_2min_6',  'auto');
results_Vpac2 = ...
   pc.pair_correlation_Veatch(XY2, [], 10, 'B20_A431_2min_12', 'auto');
