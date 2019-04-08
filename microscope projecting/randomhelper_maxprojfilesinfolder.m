savefoldername = 'G:\FromMicroscopeComputer\190118 CamKII virus test\300nL_I206K_postlight';
currname = '300nL_I206K_postlight_w1561_s3_t';
filesavefolder = fullfile(savefoldername,currname);
% filesavefolder = 'G:\FromMicroscopeComputer\190118 CamKII virus test\300nL_I206K_+488_50perpower\300nL_I206K_+488_50perpower_w2561_s2_t';

im_array = KennedyLabMicroscopeData.loadtiffseries(filesavefolder,[currname '*'],'maxproj');
        %         [image2save,~] = GeneralAnalysis.timedriftCorrect(im_array);
        image2save = im_array;
        GeneralAnalysis.LibTiff(image2save,fullfile(savefoldername,erase(currname,'*')));