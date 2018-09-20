
%         channelordering(1) = abeta channel in image
%         channelordering(2) = channel in image corresponding to 'ch1';
%         channelordering(3) = channel in image corresponding to 'ch2';

ids = [2 1 3];
channelorderingstr = {'chABeta','Gephyrin','Bassoon'}; % channel abeta, channel 1, channel2
filepath = 'C:\Users\KennedyLab\Documents\Hannah\SIM data\111617\GluA1488_PSD95561_ABeta647_100717\GluA1488_PSD95561_Abeta647_Reconstructed.nd2';

s  = SIM();
s.channelordering = ids;
s.loadNDfile(filepath);
% s.make_masks();
% s.make_distancemasks();
% s.calculateDensity_abetaINch1();
% s.calculateDensity_abetaINch2();
s.dothething;

%% -- calculate background of ~ch1 average and set all pixels less than it to be equal to mean background
s  = SIM();
s.channelordering = ids;
s.loadNDfile(filepath);
s.make_cellmask;
s.make_maskch1;


s.make_masks;
premask = s.ch1.mask;
preimage = s.ch1.image;

m = zeros(1,size(s.cellmask,3));
for ii = 1:size(s.cellmask,3)
   currframe = dip_image(s.ch1.image(:,:,ii));
   currcellmask = squeeze(~s.cellmask(:,:,ii-1));
   currcellmask = bdilation(currcellmask,2);
   m(ii) = single(sum(currframe,currcellmask,[1 2])./sum(currcellmask)); 
   currframe(currframe<m(ii)) = m(ii);
   s.ch1.image(:,:,ii) = single(currframe);
end