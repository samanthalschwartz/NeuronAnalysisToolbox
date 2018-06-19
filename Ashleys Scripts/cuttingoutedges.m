% load in the Ashley File - variable name aa
close all; clear all;
filename = 'Z:\Sam\MJK_zapERtrap_for_sam\AMB_previous\DHFR-GFP-GluA1\051718\merges\slip1\5_merge_AshleyFile.mat';
load(filename);
% find the crap
dipshow(aa.surfaceCargo.image,'log');
badframes = 71:74;
badmask = (aa.surfaceCargo.image == 0);
test = GeneralAnalysis.bwmorph_timeseries(badmask,'thicken',6);
aa.surfaceCargo.mask = aa.surfaceCargo.mask.*single(~test);
save(filename,'aa','-append')


