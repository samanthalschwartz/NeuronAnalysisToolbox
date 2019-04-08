ids = [3 1 2];
channelorderingstr = {'chABeta','PSD95ib','Gephyrin'}; % channel abeta, channel 1, channel2
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\080317\500 nM Abeta_PSD95_gephyrin\PSD95ib488_gephyrin561_Abeta647_006_Reconstructed.nd2';
% ids = [2 1 3];
% channelorderingstr = {'chABeta','PSD95','Bassoon'}; % channel abeta, channel 1, channel2
% filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\071117\PSD95_Abeta_Bassoon\PSD95488_Abeta561_Bassoon647_003_Reconstructed.nd2';

sor  = SIM();
sor.channelordering = ids;
sor.channelorderingstr = channelorderingstr;
sor.loadNDfile(filepath);
sor.make_masks;
sor.make_distancemasks;

test1 = SIM.calculateNumberDensityCOM(sor.ch1.distance_mask,sor.abeta.labeled_mask,sor.abeta.image,sor.cellmask,0:20);
test2 = SIM.calculateNumberDensityCOM(sor.ch2.distance_mask,sor.abeta.labeled_mask,sor.abeta.image,sor.cellmask,0:20);
figure;
plot(test1.d,(test1.radialnumber./test1.nummask1)./(test1.totalnumber/test1.nummask1),'r');
hold on;
plot(test2.d,(test2.radialnumber./test2.nummask1)./(test2.totalnumber/test2.nummask1),'g');
xlim([0 20])


figure;
plot(test1.d,(test1.radialnumber),'r');
hold on;
plot(test2.d,(test2.radialnumber),'g');
xlim([2 20])



figure;
plot(test1.d,(test1.cumulativeradialnumber./test1.nummask1)./(test1.totalnumber/test1.nummask1),'r');
hold on;
plot(test2.d,(test2.cumulativeradialnumber./test2.nummask1)./(test2.totalnumber/test2.nummask1),'g');
xlim([2 20])
