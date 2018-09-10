function uiopen(filename,direct)
%uiopen  overload uiopen to open files
%
% MarK J. Olah (mjo@cs.unm.edu)
% 12/15/2014

[~, ~, ext] = fileparts(filename);
orig_dir=cd;

switch ext
    case '.hsi'
        obj=HSData(filename);
        name=getObjName('hsd');
    case '.hsdata'
        obj=HSData(filename);
        name=getObjName('hsd');
    case '.spt'
        obj=SPT(filename);
        name=getObjName('spt');
    case '.tsa'
        obj=TrackSegmentAnalysis(filename);
        name=getObjName('tsa');
    case '.rpt'
        obj=RPT(filename);
        name=getObjName('rpt');
    case '.hsrpt'
        obj=HSRPT(filename);
        name=getObjName('hsrpt');
    case '.spdata'
        obj=SPData(filename);
        name=getObjName('spd');
%     case '.mat'
%         obj=SPData(filename);
%         name=getObjName('spd');
    case {'.tif','.TIF'}
        obj = loadtiff(filename);
        name = 'image';
        if ndims(obj) == 4
            cmap = {red,green,blue};
            % then it is a color image
            for ii = 1:size(obj,4)
                if ~isempty(find(obj(:,:,:,ii),1))
                    dipshow(obj(:,:,:,ii),'lin');
                    dipmapping(gcf,'Colormap',cmap{ii});
                end
            end
        end
    case '.reganalysis'
        obj = RegistrationAnalysis(filename);
        name = getObjName('reg');
    case '.regpairs'
        obj = RegisteredPairAnalysis(filename);
        name = getObjName('regPairs');
    otherwise
        try
            [next_uiopen_filename, current_uiopen_dir]=get_next_uiopen();
            next_uiopen_dir=fileparts(next_uiopen_filename);
            cd(next_uiopen_dir);
            uiopen(filename,direct);
        catch err
            disp(getReport(err));
            fprintf('Uiopen: %s failuire!\n',current_uiopen_dir);
        end
        cd(orig_dir);
        return
end
assignin('base',name, obj);
if ~strcmp(ext,'.tif') &&  ~strcmp(ext,'.TIF') 
    evalin('base',sprintf('disp(%s)',name));
end
% if strcmp(ext,'.tif')
%     GeneralAnalysis.displaytiff(obj);
% else
% evalin('base',sprintf('disp(%s)',name));
% end
fprintf('Auto Loaded File:"%s"\nFile contained %s object.\nSaved as workspace variable: "%s"\n',...
        filename, class(obj), name);
if ismethod(obj,'gui') %If there is an obj.gui method, then run it
    evalin('base',sprintf('%s.gui;',name));
end

    function [uiopen_filename,current_uiopen_dir]=get_next_uiopen()
        % find next shadowed uiopen
        alluiopen = which('uiopen','-all');
        alluiopen_dirs=cellmap(@fileparts,alluiopen);
        currentFile = dbstack('-completenames');
        currentFile(~strcmp('uiopen',{currentFile.name})) = [];
        current_uiopen_dir=fileparts(currentFile(1).file);
        if ispc()
            spch=';'; %in windows paths are seperated by ;
        else
            spch=':'; %in linux by :
        end
        ps=strsplit(path,spch);
        curr_idx=find(strcmp(current_uiopen_dir,ps),1,'first');
        for i=curr_idx+1:length(ps)
            next_idx=find( strcmp(ps(i),alluiopen_dirs),1,'first');
            if next_idx
                uiopen_filename=alluiopen{next_idx};
                return 
            end
        end
        uiopen_filename='';
    end

    function objName=getObjName(baseName)
        if ~evalin('base', sprintf('exist(''%s'',''var'')',baseName))
            objName=baseName;
        else
            idx=1;
            objName=sprintf('%s%i',baseName,idx);
            while evalin('base', sprintf('exist(''%s'',''var'')',objName))
                idx=idx+1;
                objName=sprintf('%s%i',baseName,idx);
            end
        end
    end

end

