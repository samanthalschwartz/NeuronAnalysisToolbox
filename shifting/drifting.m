filedir = 'I:\Projects\Project Cry2Olig-Gephyrin\180724 GephIntraOlig-GabamCh\coverslip1-olig-mChGaba\process\combined\'
filenameZproj = 'gephOlig-488light_TimeSeries_Zproj_s';
filename = 'gephOlig-488light_TimeSeries_s';
ff=8
uiopen(fullfile(filedir,[filenameZproj num2str(ff) '.tif']),1);
[img_out,sv_arr] = GeneralAnalysis.timedriftCorrect(image);
LibTiff(img_out,fullfile(filedir,[filenameZproj num2str(ff) '_driftcorrected.tif']))
uiopen(fullfile(filedir,[filename num2str(ff) '.tif']),1);