rt = ROITools();
rt.getROI([10000, 20000; 20000, 40000; 30000, 60000], '3x2 matrix');
rt.getROI([10000, 20000, 10, 10; ...
           20000, 40000, 10, 10; ...
           30000, 60000, 10, 10], '3x4 matrix');
load('data/Fit4');
rt.getROI(SRtest, 'SR_demo');
rt.getROI(SRtest.Results, 'SRtest.Results');
clear SRtest
rt.getROI('data/Fit4', 'SRtest file');
rt.getROI({'/home/wester/local/projects/STMC/Farzin/autophagsomes/SR_demo_Results_Label_01.mat', '/home/wester/local/projects/STMC/Farzin/autophagsomes/SR_demo_Results_Label_02.mat'}, 'SMD file');
load('/home/wester/local/projects/STMC/cluster/HeLa-Tub1-2017-8-28-13-56-9_Results.mat');
rt.getROI(SMASR, 'SMA_SR');
rt.getROI(SMASR.SMD', 'SMASR.SMD');
clear SMASR
rt.getROI({'/home/wester/local/projects/STMC/cluster/HeLa-Tub1-2017-8-28-13-56-9_Results.mat'}, 'SMA_SR file');
