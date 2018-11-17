filepath = 'G:\FromMicroscopeComputer\Sam\181114 NL1-GFP\new drug 1500nM';
% filestr = 'NL1-GFP_newzap_+AB_postrelease_timecourse_w3640_s2_t*';
filestr = 'NL1-GFP_newzap_+ABtimelapse_w2640_s5_t*'
savename = fullfile('G:\FromMicroscopeComputer\Sam\181114 NL1-GFP',[filestr(1:end-1) '_timeseries']);
option = 'maxproj';
im_array = GeneralAnalysis.loadtiffseries(filepath,filestr,option);
% im_array_z = GeneralAnalysis.loadtiffseries(filepath,filestr);
LibTiff(im_array,savename)


