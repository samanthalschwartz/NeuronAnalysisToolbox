clear all
close all

XY = load('../data/threeD.csv');

XY1 = XY(1:2:end, :);
XY2 = XY(2:2:end, :);

pc = PairCorr();
pc.Results = 'results';
pc.Rmax_axis = 100;

results_pcc  = pc.pair_correlation(10, 'threeD', XY1, XY2)
results_pac1 = pc.pair_correlation(10, 'threeD_1', XY1)
results_pac2 = pc.pair_correlation(10, 'threeD_2', XY2)

results_Vpcc  = pc.pair_correlation_Veatch(XY1, XY2, 10, 'threeD', 'cross');
results_Vpac1 = pc.pair_correlation_Veatch(XY1, [], 10, 'threeD_1', 'auto');
results_Vpac2 = pc.pair_correlation_Veatch(XY2, [], 10, 'threeD_2', 'auto');
