 NA = NeuronAnalysis();
%  NA.cutoff.xrange = 15;
%  NA.cutoff.yrange = 15;
%  NA.path_channel_cellfill = 'G:\Sam\050118 NL1 insertion\cell4-C2.tif';
  NA.path_channel_cellfill = 'G:\Sam\051318 GluA1 insertion\cell2_merge-C2.tif';
%  NA.path_channel_DHFR = 'F:\Sam\050118 NL1 insertion\cell2_merge-C2.tif';
 NA.path_channel_DHFR = 'G:\Sam\051318 GluA1 insertion\cell2_merge-C3.tif';
 NA.loadimgfrompaths();
  NA.cutoff.xrange = 195:size(NA.channel_cellfill,1)-1;
 NA.cutoff.yrange = 0:size(NA.channel_cellfill,1)-195;
 NA.cropchannels();
 NA.make_masks();
 
 h1 = NA.maskoverlay_cellfill;
 h2 = NA.maskoverlay_DHFR;
