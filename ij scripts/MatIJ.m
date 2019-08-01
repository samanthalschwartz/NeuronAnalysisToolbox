classdef MatIJ < handle
    % wrapper class for some additional useful functions when using ImageJ and MIJ functions in Matlab
    properties
        test = [];
    end
    
    methods (Static)
        function imp = showImage(I,varargin)
            % display Matlab image as ImageJ with window name the same as
            % variable name
            windowname = matlab.lang.makeValidName(inputname(1)); % this is to get the variable name
            if nargin > 1
                imp = othercopytoImagePlus(I,windowname,varargin{1});
            else
                strs = MatIJ.guessDimStr(I);
                if ~isempty(strs)
                    imp = othercopytoImagePlus(I,windowname,strs);
                else % just let default channel settings take over since too much going on
                    imp = othercopytoImagePlus(I,windowname);
                end
            end
            imp.show;
        end
        function imp = showMask(I,varargin)
            I = single(I);
            windowname = matlab.lang.makeValidName(inputname(1)); % this is to get the variable name
            imp = othercopytoImagePlus(I,windowname,varargin{1});
            imp.show;
        end
        function imnames = getOpenImageWindowNames()
            imnames = cell(ij.WindowManager.getImageTitles());
        end
        function closeAllIJWindows()
            ij.WindowManager.closeAllWindows();
        end
        function img = imgfromimplus(implus)
            % currently this only works for xyztc - if there are more
            % dimensions, need to fix this.....
            % returns Matlab Matrix image in forMat --xyztc--- from imageJ image plus (composite) ij params
            dims = implus.getDimensions'; xy_dims = dims(1:2);
            nch = implus.getNChannels;
            nframes = implus.getNFrames;
            nslice = implus.getNSlices;
            ids = dims==nch|dims==xy_dims(1)|dims==xy_dims(2)|dims==nframes|dims==nslice;
            rem_dims = dims(~ids);
            final_dims = [xy_dims, nslice, nframes, nch, rem_dims];
            img = zeros(final_dims);
            % z first
            % t second
            % c third
            stack = implus.getImageStack();
            for zz = 1:nslice
                for tt = 1:nframes
                    for cc = 1:nch
                        currframe = img(:,:,zz,tt,cc);
                        stckid = implus.getStackIndex(cc,zz,tt);
                        test = stack.getProcessor(stckid); % this is the 2d image
                        vals = test.getPixels();
                        currframe(:) = vals(:);
                        currframe = permute(currframe,[2 1]);
                        img(:,:,zz,tt,cc) = currframe;
                        
                    end
                end
            end
        end
        
        function img = pullimage(input)
            % returns image in forMat xyztc by asking for properties from image
            % in image j.
            if nargin>0
                if ischar(input)
                    disp('Pulling IJ image into Matlab....');
                    implus = ij.WindowManager.getImage(input);
                    ij.WindowManager.setCurrentWindow(ij.WindowManager.getWindow(input));
                    img = MatIJ.imgfromimplus(implus);
                    disp('Done');
                end
            end
        end
        
        function img = getImage()
            imnames = ij.getImageTitles();
            if nargout == 1
                [indx,tf] = listdlg('ListString',imnames,'SelectionMode','single');
                if tf
                    img = pullimage(imnames{indx});
                else
                    img = [];
                end
            else
                [indx,tf] = listdlg('ListString',imnames);
                if tf
                    for ii = 1:numel(indx)
                        img = pullimage(imnames{indx(ii)});
                        varname = MatIJ.cleanfilename(imnames{indx(ii)});
                        assignin('base',varname, imag);
                    end
                    img = numel(indx);
                else
                    img = [];
                end
            end
        end
        
        
        function outname = cleanfilename(inname)
            oo = strsplit(inname,{'.','-'});
            prename = oo{1};
            if ~isempty(str2num(prename(1)))
                outname = ['im_' prename];
            else
                outname = prename;
            end
        end
        function saveImage(inputimagename,savename,dimstr)
            % save an image open in the image j window to tiff or save a Matlab image
            % as a tif.
            % user is prompted to identify open image if none is specified
            % user is prompted for save file name is none is specified
            % if saving a tiff from a Matlab variable - first use showImage
            % function to display in ImageJ and then save it. optional
            % dimension str for call can be provided ie: 'xyct' or 'xyctz'
            if isempty(nargin)
                imnames = MatIJ.getOpenImageWindowNames();
                [~,tf] = listdlg('ListString',imnames);
                if ~tf
                    return;
                end
            elseif nargin<2
                [file,path,~] = uiputfile('*.tif');
                savename = fullfile(path,file);
            end
            if ischar(inputimagename)
                try
                    ij.WindowManager.setCurrentWindow(inputimagename)
                    ij.IJ.saveAs("Tiff",savename);
                catch
                    disp(['No image open named ' inputimagename]);
                end
            elseif exist(inputimagename,'var') && ndims(inputimagenames)>1
                if nargin>2 && ischar(dimstr)
                    imp = MatIJ.showImage(inputimagename,dimstr);
                else
                    imp = MatIJ.showImage(inputimagename);
                end
                IJ.saveAsTiff(imp,savename);
            end
            
        end
        function dims = setDims(imgPlus)
            
            
        end
        
        function strs = guessDimStr(I)
            strs = 'XY';
            n_dims = ndims(I);
            dims = size(I);
            
            switch n_dims
                case 2
                    strs = 'XY';
                case 3
                    if dims(3) > 4 % assume if dimension is more than 4 than prob not the channels
                        strs = 'XYZ';
                    else
                        strs = 'XYC';
                    end
                case 4
                    [B,ids] = sort(dims(3:end),'ascend'); % get the smallest valued dimension to assess if it's a color
                    if B(1) > 4 % assume there is no color channel
                        chs = {'Z','T'};
                        for ii = ids
                            strs = [strs,chs{ii}];
                        end
                    else
                        chs = {'C','T'}; % assume color channel is smaller than time channel dimension
                        for ii = ids
                            strs = [strs,chs{ii}];
                        end
                    end
                case 5
                    [~,ids] = sort(dims(3:end),'ascend');
                    chs = {'C','Z','T'}; %order of increasing dimension size (channel is smallest, then z then time)
                    for ii = ids
                        strs = [strs,chs{ii}];
                    end
                otherwise
                    strs = [];
            end
        end
    end
end