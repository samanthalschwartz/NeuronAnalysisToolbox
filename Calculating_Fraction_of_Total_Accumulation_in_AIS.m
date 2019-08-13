close all; clear all
% load in Ashley File
SurfaceCargoInAIS = aa.cleanedcargomask.*aa.cellFill.AIS_mask.*aa.surfaceCargo.image;
SurfaceCargoTotal = aa.cleanedcargomask.*aa.surfaceCargo.image;
sum(SurfaceCargoInAIS(:))/sum(SurfaceCargoTotal(:))

