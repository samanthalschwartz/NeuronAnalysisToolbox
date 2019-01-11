% filedir = 'C:\Users\sammy\Desktop\Brooke SIM\SIM_Files\071117';
% filedir = 'C:\Users\sammy\Desktop\Brooke SIM\SIM_Files\072617';
% filedir = 'C:\Users\sammy\Desktop\Brooke SIM\SIM_Files\071917';
filedirs = {'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071117',...
    'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\071917',...
     'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\072617',...
      'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\080317',...
       'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\111617',...
       'C:\Users\KennedyLab\Dropbox\Shared with Hannah\SIM data\SIM_Files\112117',...
   };
%%
for dd = 1:numel(filedirs)
    filedir = filedirs{dd};
    files = dir(fullfile(filedir,'*.mat'));
    for ff = 1:numel(files)
        clear obj
        saveinfo = fullfile(filedir,files(ff).name);
        load(saveinfo);
        %----- delete these after calling once
        obj.ch1.rawimage = obj.ch1.image;
        obj.ch2.rawimage = obj.ch2.image;
        obj.abeta.rawimage = obj.abeta.image;
        %-----
        obj.savepath = saveinfo;
        obj.setimage;
        obj.make_laplacemasks;
        obj.save;
    end
end

%%
wb = waitbar(0);
for dd = 1:numel(filedirs)
    filedir = filedirs{dd};
    files = dir(fullfile(filedir,'*.mat'));
    for ff = 1:numel(files)
        clear obj
        saveinfo = fullfile(filedir,files(ff).name);
        load(saveinfo);
        obj.make_cellmask;
        obj.make_maskchAB;
        obj.savepath = saveinfo;
        obj.save
    end
    waitbar(dd/numel(filedirs),wb)
end