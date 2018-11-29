classdef KennedyLabMicroscopeData < handle
    properties
        test = [];
    end
    
    methods (Static)
        function moveMetaMorphThumbs(startdir,idstring)
            if nargin<2
                idstring = 'thumb';
            end
            if nargin<1
                startdir = uigetdir(pwd,'Select a directory to move all ''thumb'' files to a new folder');
            end
            thumbdir = fullfile(startdir,idstring);
            if ~exist(thumbdir,'dir')
                mkdir(thumbdir)
            end
            allfiles = dir(fullfile(startdir));
            bools = cell2mat(arrayfun(@(x) (~isempty(strfind(x.name,idstring)) & ~x.isdir),allfiles,'UniformOutput',false));
            movefiles = allfiles(bools);
            for mm = 1:size(movefiles,1)
                movefile(fullfile(startdir,movefiles(mm).name),thumbdir);
            end
        end
        
        function boolarray = findfiles(allfiles,inputstr)
           boolarray = cell2mat(arrayfun(@(x) (~isempty(strfind(x.name,inputstr)) & ~x.isdir),allfiles,'UniformOutput',false));   
        end
        
    function im_array = loadtiffseries(filepath,filestr,option)
        % this function loads in a tiff series to make a movie
        % inputs:
        %       option is an action to take on the files as they come in.
        %           option = 'z-project';
        %           option = '';
        if nargin == 0 || isempty(filepath)
            filepath = pwd;
        end
        if nargin > 0 && isempty(filestr) 
            [files, filepath] = uigetfile(filepath,'*.*','Multiselect','on');
        elseif nargin >= 2
            files = dir2cell(filepath,filestr);
        end
        if nargin>2
            switch option
                case 'maxproj'
                    maxproj = @(x)(max(x,[],3));
                    img_operation = maxproj;
                case 'sumproj'
                    maxproj = @(x)(sum(x,[],3));
                    img_operation = maxproj;
            end
        end
        % make sure file names are in the correct order
        % remove files with string 'thumb'
        outids  = cell2mat( cellfun(@(x) contains(x,'thumb'),files,'UniformOutput',false) );
        files(outids) = []; 
        im_array = KennedyLabMicroscopeData.loadtiffseries_fromfiles(files,img_operation);
%         files = natsortfiles(files);
%         % first determine image size
%         path = fullfile(filepath,files{1});
%         oimg = loadtiff(path);
%         if nargin>2
%             oimg = img_operation(oimg);
%         end
%         im_array = zeros([size(oimg),numel(files)]);
%         if numel(files)>1
%             img_nd = ndims(im_array);
%             otherdims = repmat({':'},1,img_nd-1);
%             im_array(otherdims{:}, 1) = oimg;
%         else
%             im_array = oimg;
%         end
%         wb = waitbar(0,'Loading Files...');
%         for ff = 2:numel(files)
%             path = fullfile(filepath,files{ff});
%             oimg = loadtiff(path);
%             if nargin>2
%                 oimg = img_operation(oimg);
%             end
%             im_array(otherdims{:}, ff) = oimg;
%             waitbar(ff/numel(files),wb);
%         end
%         close(wb);
    end    
        
    function im_array = loadtiffseries_fromfiles(files,img_operation)
        if nargin>2
            if isstring(img_operation)
                switch option
                    case 'maxproj'
                        maxproj = @(x)(max(x,[],3));
                        img_operation = maxproj;
                    case 'sumproj'
                        maxproj = @(x)(sum(x,[],3));
                        img_operation = maxproj;
                end
            end
        end
        files = natsortfiles(files);
        % first determine image size
        oimg = loadtiff(files{1});
        if nargin>1
            oimg = img_operation(oimg);
        end
        im_array = zeros([size(oimg),numel(files)]);
        if numel(files)>1
            img_nd = ndims(im_array);
            otherdims = repmat({':'},1,img_nd-1);
            im_array(otherdims{:}, 1) = oimg;
        else
            im_array = oimg;
        end
        wb = waitbar(0,'Loading Files...');
        for ff = 2:numel(files)
            oimg = loadtiff(files{ff});
            if nargin>1
                oimg = img_operation(oimg);
            end
            im_array(otherdims{:}, ff) = oimg;
            waitbar(ff/numel(files),wb);
        end
        close(wb); 
    end
        
        
    end
    
    
end
