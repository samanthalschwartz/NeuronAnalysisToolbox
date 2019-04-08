filepath = 'C:\Users\schwsama\Documents\Data\181115 HaloPdisplay';
% filestr = 'NL1-GFP_newzap_+AB_postrelease_timecourse_w3640_s2_t*';
for ii = 2:9
filestr = ['TfR-Halo_mChCellFill_6nMhalo660_longExp_timecourse_w2640_s' num2str(ii) '_t*'];
savename = fullfile('C:\Users\schwsama\Documents\Data\181115 HaloPdisplay',[filestr(1:end-1) '_timeseries']);
option = 'maxproj';
im_array = GeneralAnalysis.loadtiffseries(filepath,filestr,option);
% im_array_z = GeneralAnalysis.loadtiffseries(filepath,filestr);
LibTiff(im_array,savename);
disp(['Finished ' num2str(ii)]);
end

filepath = 'Z:\Sam\Data!\181121 Halo-pDisplay Kinetics\pDisplay_Halo_10nM\pDisplay_Halo_10nMAF660_timecourse_fr4dye_w2640';
files = dir(fullfile(filepath,'pDisplay*'));
files = natsortfiles(files);
option = 'maxproj';
savename = fullfile(filepath,[files(1).name(1:end-4) '_timeseries']);
for ff = 1:numel(files)
    filestr = files(ff).name;
    currim = GeneralAnalysis.loadtiffseries(filepath,filestr,option);
    if ff == 1
        im_array = zeros(size(currim,1),size(currim,2),numel(files));
    end
    im_array(:,:,ff) = currim;
    disp(['Finished ' num2str(ff)]);
end
LibTiff(im_array,savename);


