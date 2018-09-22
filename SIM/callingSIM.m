
%         channelordering(1) = abeta channel in image
%         channelordering(2) = channel in image corresponding to 'ch1';
%         channelordering(3) = channel in image corresponding to 'ch2';

ids = [3 1 2];
channelorderingstr = {'chABeta','PSD95i','Gephyrin'}; % channel abeta, channel 1, channel2
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\080317\500 nM Abeta_PSD95_gephyrin\PSD95ib488_gephyrin561_Abeta647_001_Reconstructed.nd2';

ids = [2 1 3];
filepath = 'G:\SIM\180704 Coverslip3\Coverslip3_Reconstructed-cell1pos3.nd2';
channelorderingstr = {'antiGFP','Gephyrin','antiGephyrin'}; % channel abeta, channel 1, channel2
s  = SIM();
s.channelordering = ids;
s.channelorderingstr = channelorderingstr;
s.loadNDfile(filepath);
% s.make_masks();
% s.make_distancemasks();
% s.calculateDensity_abetaINch1();
% s.calculateDensity_abetaINch2();
s.make_masks;
s.make_cellmask;
s.make_maskABeta;
s.ch1.mask = s.make_maskSynapseMarkerLow(s.ch1.image);
s.ch2.mask = s.make_maskSynapseMarkerLow(s.ch2.image);
s.make_distancemasks();
s.calculateDensities_abetaINch1();
s.calculateDensities_abetaINch2();
s.measure_allthings();


% s.dothething;
%%
% make some plots
output = s.calculateDensities_abetaINch(s.ch2,s.ch2.mask,s.ch2.image,s.cellmask,0:30)
figure; plot(output.radialdensity_norm(:,1),output.radialdensity_norm(:,2),'LineWidth',3);
xlabel('Distance from Gephyrin Mask'); ylabel('Relative Intensity Density (AU/voxel)');
set(gca,'FontSize',16)
hold on;
plot(s.ch1.abetaCh.radialdensity_norm(:,1),s.ch1.abetaCh.radialdensity_norm(:,2),'LineWidth',3);

figure; plot(output.cumulative_radialnumberdensity_norm_raw(:,1),output.cumulative_radialnumberdensity_norm_raw(:,2),'LineWidth',3);
xlabel('Distance from Gephyrin Mask'); ylabel('Relative Number of ABeta Aggregates (#/voxel)');
set(gca,'FontSize',16)
hold on;
plot(s.ch1.abetaCh.radialdensity_norm(:,1),s.ch1.abetaCh.radialdensity_norm(:,2),'LineWidth',3);

%%
% shift image and calculate density
newimage = zeros(size(s.ch1.image,1),10,size(s.ch1.image,3));
shiftch1 = cat(2,s.ch1.image(:,5:end,:),newimage);
output_shift = s.calculateDensities_abetaINch(s.ch1,s.ch1.mask,shiftch1,s.cellmask,0:30);
figure; plot(output_shift.radialdensity_norm(:,1),output_shift.radialdensity_norm(:,2),'LineWidth',3);
