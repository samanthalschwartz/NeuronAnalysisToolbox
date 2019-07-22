% topdir = 'F:\GephIntrabody Project\180724 GephIntraOlig-GabamCh\coverslip2-olig-mChGaba';
clear all
topdir = 'F:\FromMicroscopeComputer\180808 GlyR-GephOlig';
KennedyLabMicroscopeData.moveMetaMorphThumbs(topdir,'thumb')
remainingfiles = dir(fullfile(topdir,'*.tif'));
while ~isempty(remainingfiles)
    namestr = remainingfiles(1).name;
    % namestr = 'TfR-Halo_mChCellFill_6nMhalo660_longExp_timecourse_w1561_s4_t6';
    chstr = '_w  
    % different folder
    timeid = cellfun(@(x) ~isempty(strfind(x,timestr)),matches);
    timepos = out(timeid);
    if isempty(timepos) % if there aren't multiple times, just move file into savefoldername
        im_array = KennedyLabMicroscopeData.loadtiffseries(topdir,namestr,'maxproj');
        maxfold = fullfile(topdir,'maxfold');
        if ~exist(maxfold)
            mkdir(maxfold);
        end
        GeneralAnalysis.LibTiff(im_array,fullfile(fullfile(maxfold,[namestr(1:end-4) '_maxproj'])));
        srcfile = fullfile(topdir,namestr);
        destfile = fullfile(savefoldername,namestr);
        movefile(srcfile,destfile)
    else
        if isempty(chid) % there are multiple times and not multiple channels
            currname = namestr(1:timepos+1);
        else % there are multiple times and multiple channels
            if timepos<chpos % Andor File
                currname =  [namestr(1:timepos+1),'*',namestr(chpos:end-4)];
            else % Metamorph File
                currname = namestr(1:timepos+1);
            end
        end
        filesavefolder = fullfile(savefoldername,erase(currname,'*'));
        if ~exist(filesavefolder)
            mkdir(filesavefolder)
        end
        files = dir2cell(topdir,[currname '*']);
        for ff = 1:numel(files)
            [filepath,name,ext] = fileparts(files{ff});
%             try
                srcfile = fullfile(filepath,[name,ext]);
                destfile = fullfile(filesavefolder,[name,ext]);
                if length(destfile)>255
%                     if ispc
% %                       srcfile = strrep(srcfile,'+','-');
%                         destfile = namestr(min(regexp(namestr,'_[wt]')):end);
%                     end
                end
                movefile(srcfile,destfile)
%             catch
%                 display(['no file called ' srcfile]);
%             end
        end
        im_array = KennedyLabMicroscopeData.loadtiffseries(filesavefolder,[currname '*'],'maxproj');
        %         [image2save,~] = GeneralAnalysis.timedriftCorrect(im_array);
        image2save = im_array;
        GeneralAnalysis.LibTiff(image2save,fullfile(savefoldername,erase(currname,'*')));
    end
    remainingfiles = dir(fullfile(topdir,'*.tif'));
end
