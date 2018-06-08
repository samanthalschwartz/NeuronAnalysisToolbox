%% Gephyrin
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\intrabody-gabaA2\intrabody-gabaA2_2_FRAP488_20180601_122100 PM';
str = '*w0001.tif';
% str = '*w0000.tif';
fp.load(datadir,str);
fp.ROIs = {[110, 203],[186, 188],[382, 328],[261, 301]};
fp.inc = 8;
% fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 10;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();

%% GABA
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\intrabody-gabaA2\intrabody-gabaA2_2_FRAP561_20180601_122036 PM';
str = '*w0000.tif';
fp.load(datadir,str);
fp.ROIs = {[110, 203],[186, 188],[382, 328],[261, 301]};
fp.inc = 8;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 10;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%% gamma 2 gephryn
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\intrabody-gabaG2\intrabody-gabaG2_2_FRAP488-longer_20180602_80535 PM';
str = '*w0001.tif';
fp.load(datadir,str);
fp.ROIs = {[116, 454],[237, 173],[233, 454],[318, 377]};
fp.inc = 15;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%%
%% gamma 2 gaba same cell - still blecahed by 488
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\intrabody-gabaG2\intrabody-gabaG2_2_FRAP488-longer_20180602_80535 PM';
str = '*w0000.tif';
fp.load(datadir,str);
fp.ROIs = {[116, 454],[237, 173],[233, 454],[318, 377]};
fp.inc = 8;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%%
fp = FRAP(); %super zipping around Gaba but not intrabody

datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\intrabody-gabaG2\intrabody-gabaG2_2_FRAP561-longer_20180602_70655 PM';
str = '*w0000.tif';
% str = '*w0001.tif';
fp.load(datadir,str);
fp.ROIs = {[83, 252],[175, 97],[241, 382],[171, 278]};
fp.inc = 8;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%% gaba g2 shorter
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\intrabody-gabaG2\intrabody-gabaG2_2_FRAP561-shorter_20180602_70818 PM';
str = '*w0000.tif';
fp.load(datadir,str);
fp.ROIs = {[83, 252],[175, 97],[241, 382],[171, 278]};
fp.inc = 15;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%% cry2 gaba
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\cy2gephryn-gabaA2\intrabody-gabaG2_2_bluelight_FRAP561_20180602_92322 PM';
str = '*t00*.tif';
fp.load(datadir,str);
fp.ROIs = {[305, 205],[49, 38],[330, 317],[407, 294]};
fp.inc = 8;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%%
%% cry2 gaba
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\cy2gephryn-gabaA2\intrabody-gabaG2_3_bluelight_FRAP561_20180602_101332 PM';
str = '*w0000.tif';
fp.load(datadir,str);
% fp.ROIs = {[305, 205],[49, 38],[330, 317],[407, 294]};
fp.inc = 15;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
%% cry2 gaba g2
fp = FRAP();
datadir = 'G:\Sam\Data\180601 GABAR + Gephryn\cy2gephryn-gabaG2\cry2gephrin-gabaG2_1_bluelight_561FRAP_20180602_111114 PM';
str = '*w0000.tif';
fp.load(datadir,str);
fp.ROIs = {[305, 205],[49, 38],[330, 317],[407, 294]};
fp.inc = 15;
fp.showROIs;
% fp.adjustROIs;
fp.bleachframes = 5;
fp.framerate = 30;

fp.getBackgroundROI;
fp.getControlROI;
% fp.setROIs()
% [a,b] = fp.makeRecoveryCurve();
fp.setROIs;
fp.plotRecoveryCurves();
