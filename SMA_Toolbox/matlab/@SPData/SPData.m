% SPData.m
% Mark J. Olah (mjo@cs.unm.edu)
% 10/13/14
%

classdef SPData < BaseData
    %
    % For now we handle only 2D EMCCD camera data
    properties (Constant=true, Hidden=true)
        %Inherited from Pickle
        saveFileExt='.spdata'; % The extension used for saved files
        SaveableDataFormats={'*.spdata', 'SPData (.spdata)'};
        LoadableDataFormats={'*.mat;*.tif;*.ics;*.spdata;*.spt','All Loadable Sources';...
                             '*.spdata','SPData file'; '*.spt', 'SPT object files (.spt)';...
                             '*.mat','Mat captures (.mat)';'*.tif', 'Tiff Capture format (.tif)';'*.ics', 'ICS Capture format (.ics)'};
   
        
        RawDataFormats={'.mat', '.tif', '.ics'};
        SaveableRawDataFormats={'*.mat','MAT CCD Capture dipimage data (.mat)'};
        LoadableRawDataFormats={'*.mat;*.tif;*.ics','CCD Capture data (.mat, .tif, .ics)'};
    end

    properties 
        size; % [Y X T] in pixels
        pixelSize; % [sizeX, sizeY] in microns
        psf; % [sigmaX, sigmaY] in pixels
        frameT; % frame time step in s

        CCDBackground; % The EMCCD camera background constant
        CCDGain; % The EMCCD camera gain constant in [e-/ADU]
       
        Paths=struct(... %All filenames are relative paths from the workingDir
             'saveFile',[],... %filename(w/extension) of the .spdata file where this will save to
             'rawData',[],... % base name for raw image file (w/ extension).
             'CCDBackground',[],...  %filename(w/extension) where the bgMean information was read
             'CCDGainCal',[],... %filename(w/extension) where the gain calibration (beads) information was read
             'SPT','SPT',... % Path to files for SPT objects using this data
             'SPTParams',[],...
             'RPT','RPT',...% Path to files for RPT objects using this data
             'TrackSegmentAnalysis','TrackSegmentAnalysis'... % Path to files for TrackSegmetnAnalysis objects using this data
             );
        AcquisitionParams=struct(...
            'acqDateTime', [],...       %Vector datetime of sequence acquisition
            'gainCalDateTime', [] ...  %Vector datetime of gain calibration
            );

        sumImage;
    end

    properties (Hidden=true)
        version=1; %For future file format version changes
    end

    properties (Dependent=true)
        nFrames; %number of frames
        sizeX;  %number of pixels in the X (up/down) dimension
        sizeY;  %number of pixels in the Y (left/right) dimension
        frameSize; %[sizeY sizeX]
        nFramePixels; % sizeX*sizeY
        nExperimentPixels; % sizeX*sizeY*sizeT
        experimentT; % experiment duration (s) [takes into account globalTBounds
        
        rawDataPath;  % Full path to the raw data file
        rawDataFormat; %The extension of the raw data file
        
        calibrated; %Boolean if calibrated?
    end

    properties (Dependent=true, Hidden=true)
        sequence;%The raw data as a dip_image object.  This conceptually mirrors the 'sequence' name in the .mat file.  
    end

    properties (Access=protected, Transient=true)
        gainCalGuiFig;
        
        %Lazily loaded properties
        CCDBackground_frames_; % [Y X T] 
        CCDGainCal_frames_; % [Y X T] 
        %Lazily loaded property flags
        CCDBackground_frames_loaded=false;
        CCDGainCal_frames_loaded=false;
    end

    methods
        function obj=SPData( varargin )
            % This constructor calls the obj.load method with the arguments
            % given.  Empty input arguments give an empty object, otherwise the load
            % method determines the valid input arguments
            if nargin>0
                obj.load(varargin{:});
            end
        end

        function load(obj, data, varargin)
            %
            % the load method overwrites the stored properties of this object,
            % loading in properties from the 'data' parameter.  The 'data'
            % paremeter can be one of several different types allowing loading from
            % different data sources.  The end goal is to create an object with
            % as much information as possible gathered from the .data source.
            % At a minimum the data source will define the source of rawdata
            % that this SPData object refers to.  The filename of the .spdata
            % file to save to will be in the same directory that the raw data
            % was found.  To save the object to the .spdata file, use the
            % .save() method.  To save the object under a new name or location,
            % use the .saveas() method, which will casue a copy of the data if
            % the directory is different.   In the end, the .spdata file *must*
            % be in the same directory as the raw data source it points to.
            %
            % IN:
            %   data - Some sort of data input [See below for options]
            %   spdatafilename - The filename where the spdata is to be saved
            %               [used only when initializing from non-saved data 
            %                 (matlab array or dip_image object)]
            %   defaultParameters - Optional last argument.  Should be a struct with
            %     values for default parameters.  Use obj.getParamStruct() to get a param
            %     struct suitable for input.
            % 
            % Depending on what 'data' parameter is, the spdatafilename is either 
            % unused, optional, or required.
            %
            % There are 4 cases I-IV of how a SPData object can be loaded.
            % 
            % Case I: Loading from saved .spdata
            %   spdatafilename parameter is not used, we just load from saved
            %   filename.  
            % Applies to data input formats:
            %   (1) .spdata filename
            %
            % Case II: Loading from SPT
            %   The spdatafilename parameter is not used.  The .spdata filename
            %   of this object will match directory and base name  of the .mat file linked
            %   from the SPT object.  If there is not enough information in the
            %   .spt to locate the .mat file the user will be prompted to
            %   loacate the data or an error will be raised.
            % Applies to data input formats:
            %   (2) SPT object
            %   (3) .spt filename
            %
            % Case III: Loading from saved image file
            %   The spdatafilename parameter is not used. The .spdata filename of this object
            %   will match the input directory and base name of the raw data file.
            % Applies to data input formats:
            %   (4) .mat filename (sequence capture structure)
            %   (5) .tif filename
            %   (6) .ics filename or filename sequence
            %
            % Case IV: Loading from raw image object directly
            %   The spdatafilename is required.  This will cause the data to be
            %   saved in .mat format matching the spdatafilename.
            % Applies to data input formats:
            %   (7) dip_image object
            %   (8) matlab array
            %
            
            spdatafile=[];
            if length(varargin)==2
                if ischar(varargin{1})
                    spdatafile=varargin{1};
                else
                    warning('SPData:load','Exepected .spdata file');
                end
                if isstruct(varargin{2})
                    obj.preservedProperties=varargin{2};
                else
                    warning('SPData:load','Exepected struct.');
                end
            elseif length(varargin)==1
                if ischar(varargin{1}) %Arg is an spdatafile
                    spdatafile=varargin{1};
                elseif isstruct(varargin{1}) %Arg is default parameters
                    obj.preservedProperties=varargin{1};
                end
            end
                    
            if nargin==1
                [loadfile,loadpath]=uigetfile(obj.LoadableDataFormats,'Select Data Source',obj.saveFilePath);
                if ~loadfile; return; end
                data=fullfile(loadpath,loadfile);
            end
            if ischar(data)
                %data is a filename
                filepath=data;
                if ~exist(filepath,'file')
                    error('SPData:FileNotFound', 'unable to open file: %s', filepath);
                end
                if ~isempty(spdatafile) 
                    %we do not use the spdatafilename argument if the data
                    %already allows us to infer a location and filename for the
                    %raw data
                    warning('SPData:load','additional .spdata filename argument ignored')
                end
                [workingDir, ~, ext]=fileparts(filepath); %#ok<*PROP>
                switch ext
                    case SPData.saveFileExt %case I       
                        s=load(filepath,'-mat');
                        if obj.initialized
                            obj.savePreservedProperties();
                            pp=obj.preservedProperties;
                            obj.resetObject();
                            obj.preservedProperties=pp;
                        end
                        obj.copyobj(s.obj);
                        obj.workingDir=collapsepath(workingDir);
                        obj.initialized=true;
                        obj.dirty=false;
                    case '.spt' %case II
                        obj.loadSPT(SPT(filepath));
                    case '.mat' %case III
                        obj.loadMatCapture(filepath);
                    case '.tif' %case III
                        obj.loadTifCapture(filepath);
                    case '.ics' %case III
                        error('SPData:load','Load ICS not implemented');
                    otherwise
                        error('SPData:UnknownFileType', 'Unknown file type: %s', ext);
                end
            elseif isa(data,'SPT')
                %case II
                if ~isempty(spdatafile) 
                    warning('SPData:load','additional .spdata filename argument ignored')
                end
                obj.loadSPT(data);
            elseif isnumeric(data)
                %case IV
                obj.loadDipImage(dip_image(data), spdatafile);
            elseif isa(data,'dip_image')
                %case IV
                obj.loadDipImage(data, spdatafile);
            else
                error('SPData:UnexpectedInputType','Cannot load type: %s', class(data));
            end
        end
       
        function saveDataAs(obj, newPath, overwrite)
            % This is how data is renamed
            obj.assertInitialized();
            if nargin == 1
                %Allow selection from a dialog
                % Make a sensible default name in defaultName
                if isempty(obj.rawDataPath)                    
                    defaultName=fullfile(obj.workingDir, [obj.savFileBaseName '.mat']);
                else
                    [rawDataDir, rawDataBaseName, ~] = fileparts(obj.rawDataPath);
                    defaultName=fullfile(rawDataDir, [rawDataBaseName '.mat']);
                end
                [savefile,savepath] = uiputfile(obj.SaveableRawDataFormats,'Save Data As...',defaultName);
                if ~savefile; return; end
                newPath = fullfile(savepath,savefile);
                if exist(newPath,'file')
                    question = sprintf('Data file "%s" exists.  Saving will overwite this file.  Continue with overwrite?', newPath);
                    overwrite = strcmp('Yes',questdlg(question,'Confirm Data Overwrite','No'));
                end
            elseif nargin == 2
                overwrite = false;
            end
            
            oldPath = obj.rawDataPath();
            if exist(newPath,'file') && (strcmp(newPath,oldPath) && ~overwrite)
                return
            end
            if ~obj.raw_frames_loaded
                obj.load_raw_frames
            end
            obj.updateWaitbar(0,'Saving Raw Data...');
            obj.saveRawData(newPath); %write-out
            obj.updateWaitbar(1);
        end


        function convertToMat(obj, newName)
            obj.assertInitialized();
            if strcmp(obj.rawDataFormat, '.mat')
                error('SPData:convertToMat','Data is already in ".mat" format.');
            end
            if nargin == 1
                %Make the default path be a .mat version of the original rawDataFile name
                [path, name, ~]=fileparts(obj.rawDataPath);
                newPath=fullfile(path, [name '.mat']);
                obj.saveDataAs(newPath);
            else
                [newBasePath, newBaseName, ~]=fileparts(newName);
                if ~isempty(newBasePath)
                    newPath=fullfile(newBasePath, [newBaseName '.mat']);
                elseif ~isempty(obj.rawDataPath)
                    [rawDataDir, ~, ~] = fileparts(obj.rawDataPath);
                    newPath=fullfile(rawDataDir, [newBaseName '.mat']);
                else
                    newPath=fullfile(obj.workingDir, [newBaseName '.mat']);
                end
                obj.saveDataAs(newPath);
            end
        end

        function [frames,roi]=getRawFrames(obj,roi_in)
            % Get raw camera data as unaltered uint16 matlab array.  This uses the cached copy
            % so we only ever load once.
            %
            % If we ask for whole dataset, we ignore global time bounds
            % here.  This will allow you to get all the frames including
            % those excluded by the time bounds.  This enable the getFrames
            % method to work correctly, as it must have access to all of
            % the frames to correctly slice into them.
            % In:
            %   roi_in - [optional] A 1x6 ROI or a index into the local ROI cell array
            % Out:
            %   frames - uint16 array size: [Y X T] 
            %   roi - 1x6 array: [xmin, xmax, ymin, ymax, tmin, tmax]
            obj.assertInitialized()
            if nargin==1; roi_in=[]; end
            if ~obj.raw_frames_loaded 
                obj.load_raw_frames() %Try to load.  Will throw exception if there is a problem
            end
            if isempty(roi_in) 
                %Return the whole dataset if no ROI given.  This is different behavior than
                %obj.getFrames() which obeys globalTBounds.
                roi=[1 obj.sizeX 1 obj.sizeY 1 obj.size(3)];
                frames=obj.raw_frames_;
            else
                roi=obj.getROI(roi_in);
                frames=obj.raw_frames_(roi(3):roi(4), roi(1):roi(2), roi(5):roi(6)); %Y X T
            end
        end

        function [frames,roi]=getFrames(obj, roi_in)
            % Get gain corrected frames as matlab single array.  This uses the cached copy of
            % the gain-corrected frames, so we only have to gain correct once.
            % In:
            %   roi_in - A 1x6 ROI or a index into the local ROI cell array
            % Out:
            %   frames - single array size: [Y X T] 
            %   roi - 1x6 array: [xmin, xmax, ymin, ymax, tmin, tmax]
            obj.assertInitialized()
            if nargin==1; roi_in=[]; end
            if ~obj.calibrated
                error('SPData:getFrames','SPData is not calibrated.  Set calibration first.');
            elseif ~obj.frames_loaded 
                obj.load_frames() %Try to load.  Will throw exception if there is a problem
            end
            roi=obj.getROI(roi_in);
            frames=obj.frames_(roi(3):roi(4), roi(1):roi(2), roi(5):roi(6)); %Y X T
        end
        

        function imref = getImRef(obj, roi_in)
            if nargin==1; roi_in=[]; end
            roi=obj.getROI(roi_in);
            sizeX=roi(2)-roi(1)+1;
            sizeY=roi(4)-roi(3)+1;
            pS = obj.pixelSize .* ones(1,2);
            imref = imref2d([sizeY, sizeX], pS(1), pS(2));
        end

        function viewRawFrames(obj, roi_in)
            if nargin==1; roi_in=[]; end
            frames=obj.getRawFrames(roi_in);
            name=['Raw Frames - ' obj.makeROITitleString(roi_in)];
            obj.viewMaximizedDipFig(dip_image(frames),'Name',name);
        end

        function viewFrames(obj, roi_in)
            if nargin==1; roi_in=[]; end
            frames=obj.getFrames(roi_in);
            name=['Calibrated Frames - ' obj.makeROITitleString(roi_in)];
            obj.viewMaximizedDipFig(dip_image(frames),'Name',name);
        end

        function [im,roi]=getSumImage(obj, roi_in)
            if nargin==1; roi_in=[]; end
            roi=obj.getROI(roi_in);
            if roi(5)==1 && roi(6)==obj.nFrames
                im=obj.sumImage(roi(3):roi(4), roi(1):roi(2)); %Y X T
            else
                im=sumImage2D(obj.getFrames(roi));
            end
        end
        
        function imH=plotSumImage(obj, roi_in)
            if nargin==1
                im=obj.sumImage;
            else
                im=obj.getSumImage(roi_in);
            end
            imH=imagesc([.5,size(im,2)-.5],[.5,size(im,1)-.5],im);
            xlabel('X (px)');
            ylabel('Y (px)');
            colormap(gray);
        end

        %% EMCCD Gain Calibration
        function frames=getCCDGainCalFrames(obj)
            % The gain calibration sequence
            if ~obj.CCDGainCal_frames_loaded 
                obj.load_CCDGainCalFrames() %Try to load.  Will throw exception if there is a problem
            end          
            frames=obj.CCDGainCal_frames_;
        end

        function frames=getCCDBackgroundFrames(obj)
            % The Background calibration sequence
            if ~obj.CCDBackground_frames_loaded 
                obj.load_CCDBackgroundFrames() %Try to load.  Will throw exception if there is a problem
            end
            frames=obj.CCDBackground_frames_;
        end

        function setCCDGainCalImage(obj, impath)
            frames=obj.loadRawDipImage(impath);
            if obj.initialized
                %only set paths if we are initialized.  This allows unitialized SPData's
                %to use the gain cal GUI, without having to specify an analysisPath
                [~,imname,ext]=fileparts(impath);
                obj.Paths.CCDGainCal=fullfile('Calibration',[imname ext]);
                GUIBuilder.safeCopyfile(impath, obj.getFilePath('CCDGainCal'),0);
            end
            %Methods aboved will throw an error if there is a problem before we get here
            %Now it is safe to set the frame cache
            obj.CCDGainCal_frames_=frames;
            obj.CCDGainCal_frames_loaded=true;
        end

        function setCCDBackgroundImage(obj, impath)
            frames = obj.loadRawDipImage(impath);
            if obj.initialized
                %only set paths if we are initialized.  This allows unitialized SPData's
                %to use the gain cal GUI, without having to specify an analysisPath
                [~,imname,ext] = fileparts(impath);
                obj.Paths.CCDBackground = fullfile('Calibration',[imname ext]);
                GUIBuilder.safeCopyfile(impath, obj.getFilePath('CCDBackground'),0);
            end
            if all(frames(:,:,1)==0) %For some reason we often capture ablank initial frame in the sequence.
                frames = frames(:,:,2:end);
            end
            %Methods aboved will throw an error if there is a problem before we get here
            %Now it is safe to set the frame cache
            obj.CCDBackground_frames_=frames;
            obj.CCDBackground_frames_loaded=true;
        end

        function clearCCDCalibrationImages(obj)
            obj.Paths.CCDBackground=[];
            obj.Paths.CCDGainCal=[];
            obj.CCDBackground_frames_=[];
            obj.CCDBackground_frames_loaded=false;
            obj.CCDGainCal_frames_=[];
            obj.CCDGainCal_frames_loaded=false;            
        end

        function recalibrateFrames(obj, offset, gain)
            % This recalibrates the frames using the internal obj.CCDBackground and
            % obj.CCDGain values.  This should be called after the values are modified.
            if nargin==1
                offset = obj.CCDBackground;
                gain = obj.CCDGain;
            end
            if ~isempty(offset) && ~isempty(gain) && offset>=0 && gain>0
                 %Get all frames even outside global time bounds
                obj.frames_ = (single(obj.getRawFrames())-offset)*gain;
                obj.frames_loaded = true;
            else
                error('SPData:recalibrateFrames','CCD Parameters not valid');
            end
            if nargin>1
                %Only here are we actually changing any gain settings and actually need to recompute
                %the sumImage and set the dirty flag
                obj.CCDBackground = offset;
                obj.CCDGain = gain;
                obj.sumImage = sumImage2D(obj.frames_);
                obj.dirty = true;
            end
                
        end

        function figHs=calibrate(obj, gain_path, bg_path)
            % Calls cal_readnoise with the saved image files or with those
            % specified as optional parameters.  Sets the CCD gain calibration
            % paramters and gain correct the raw data frames.
            %
            % gain_path - [Optional] path to the gain cal image. Set to empty to
            %             ignore.
            % bg_path - [Optional] path to the background image. Set to empty to
            %             ignore.
            obj.updateWaitbar(0,'Calibration');
            obj.updateWaitbar(0.2,'Loading Background');
            try
                if nargin==3 && ~isempty(bg_path)
                    obj.setCCDBackgroundImage(bg_path)
                end
                obj.updateWaitbar(0.4,'Loading Gain Calibration');
                if nargin>=2 && ~isempty(gain_path)
                    obj.setCCDGainCalImage(gain_path)
                end
                bg_frames=obj.getCCDBackgroundFrames();
                gain_frames=obj.getCCDGainCalFrames();
                if isempty(bg_frames) || isempty(gain_frames)
                    obj.updateWaitbar(1);
                    error('SPData:calibrate','Unable to load calibration images');
                end
                if ismatrix(bg_frames)
                    obj.updateWaitbar(1);
                    error('SPData:calibrate','CCDBackground Image has only a single frame');
                end
                if ismatrix(gain_frames)
                    obj.updateWaitbar(1);
                    error('SPData:calibrate','CCDGainCal Image has only a single frame');
                end
                %Make sizes equal
                if ~all(size(bg_frames)==size(gain_frames))
                    newsize = min(size(bg_frames),size(gain_frames));
                    bg_frames = bg_frames(1:newsize(1),1:newsize(2),1:newsize(3));
                    gain_frames = gain_frames(1:newsize(1),1:newsize(2),1:newsize(3));
                end

                obj.updateWaitbar(0.6,'cal_readnoise');

                warning('off','curvefit:fit:mismatchedOptions');
                fs = findobj('Type','figure','NumberTitle','on'); %record current figs
                out = cal_readnoise(gain_frames, bg_frames);
                figHs = setdiff(findobj('Type','figure','NumberTitle','on'),fs); %Look for new figs
                if ~isempty(figHs) %Reposityion the figures
                    movegui(figHs(1),[30 30]);
                    movegui(figHs(2),[30 600]);
                end
                warning('on','curvefit:fit:mismatchedOptions');

                obj.updateWaitbar(0.8,'Recalibrate Raw Data');
                obj.recalibrateFrames(out(4), out(2));
            catch err
                obj.updateWaitbar(1);
                disp(getReport(err))
                throw(err)
            end
            obj.updateWaitbar(1);
        end

        %% SPT Tracking
        function setSPTParams(obj,sptparams_file)
            paramspath = fullfile(obj.workingDir,'SPT','Params');
            if nargin==1
                sptparams_file = Pickle.selectExistingFileName(paramspath,...
                        '*.mat',{'*.mat','SPT Params Files (.mat)'}, 'Load SPT Params File');
            end
            if ~exist(sptparams_file,'file')
                error('SPData:setSPTParams','Unable to open SPTParams file "%s"',sptparams_file);
            end
            %copy to a local file name
            newpattern = sprintf('SPTparams.%s.mat',obj.saveFileBaseName);
            newfile = Pickle.findUnusedFileName(paramspath,  newpattern);
            GUIBuilder.safeCopyfile(sptparams_file,newfile);
            rpath = relativepath(obj.workingDir, newfile);
            if ~strcmp(obj.Paths.SPTParams, rpath)
                obj.Paths.SPTParams = rpath;
                obj.dirty = true;
            end
        end

        function sptobj=trackSPT(obj, roi_idx, sptFilename)
            % If sptFilename is provided but empty then look for a new file
            % otherwise try to load existing file
            obj.assertInitialized()
            if nargin<2
                roi_idx = 1;
            else
                [~, roi_idx] = obj.getROI(roi_in);
            end
            if isempty(roi_idx) || roi_idx>length(obj.ROI)
                error('SPData:trackSPT','Bad ROI specified: %s',roi_in);
            end
            obj.checkComplete(); % Check that the data entries are complete
            
            if ~strcmp(obj.rawDataFormat,'.mat')
                error('SPData:trackSPT','Raw Data is in "%s" format.  Must convert to .mat for SPT.',obj.rawDataFormat);
            end
            roi = obj.ROI{roi_idx};
            roi_name = obj.ROIname{roi_idx};
            sptpath = obj.getFilePath('SPT');
            paramsfilepath = obj.getFilePath('SPTParams');
            roi_frames = obj.getFrames(roi_idx);
            spt_file_pattern = sprintf('%s_%s.spt',obj.saveFileBaseName, roi_name);
            if nargin == 3 && isempty(sptFilename)
                sptFilename = Pickle.selectUnusedFileName(sptpath, spt_file_pattern,...
                               {'*.spt','SPT Files (.spt)'}, 'Select an New SPT Filename');
                if isempty(sptFilename)
                    sptobj=[];
                    return; 
                end
            elseif nargin == 2 || isempty(sptFilename)
                %Now try to load existing file
                sptFilename = Pickle.selectExistingFileName(sptpath, spt_file_pattern,...
                               {'*.spt','SPT Files (.spt)'}, 'Load an Existing SPT Filename');
                if isempty(sptFilename)
                    sptobj=[];
                    return; 
                end
            end                                
            if exist(sptFilename,'file')
                %open existing file
                sptobj = SPT.loadFile(sptFilename, obj.rawDataPath);
                sptobj.DataFile = obj.rawDataPath;
                sptobj.SaveDir = sptpath;
                sptobj.Data = dip_image(roi_frames); %preload the data
                
                if exist(paramsfilepath,'file')
                    sptobj.ParamSaveFile = paramsfilepath;
                end
            else
                %create new file
                %Initialize New Object from paramters file if it exists
                [~,baseName,~]=fileparts(sptFilename);
                i=regexp(baseName,sprintf('%s-',obj.saveFileBaseName),'end');
                appendix = baseName(i+1:end);
                sptobj = SPT(paramsfilepath);
                sptobj.DataFile = obj.rawDataPath;
                sptobj.ParamsGeneral.SaveAppendix = appendix;
                sptobj.SaveDir = sptpath;
                if ~isempty(obj.psf); sptobj.ParamsGeneral.Psf = obj.psf; end
                if ~isempty(obj.CCDBackground); sptobj.ParamsGeneral.CCDOffset = obj.CCDBackground; end
                if ~isempty(obj.CCDGain); sptobj.ParamsGeneral.Gain = 1/obj.CCDGain; end %We use -e/ADU not ADU/e-
                if ~isempty(obj.frameT); sptobj.ParamsGeneral.TimeStep = obj.frameT; end
                if ~isempty(obj.pixelSize); sptobj.ParamsGeneral.PixelSize = obj.pixelSize; end
                sptobj.ParamsGeneral.Roi = roi(1:4)-1;
                if roi(5)~=1 || roi(6)~=obj.nFrames
                    sptobj.ParamsGeneral.Frames = (roi(5):roi(6));
                end
                if exist(paramsfilepath,'file')
                    sptobj.ParamSaveFile = paramsfilepath;
                end
                sptobj.Data=dip_image(roi_frames); %preload the data
                %These Stats would normally be set in SPT.loadData, but we
                %skip that step.
                sptobj.Stats.Data.size = obj.getROIsize(roi);
                sptobj.Stats.Data.PixelSize = obj.pixelSize;
                sptobj.Stats.Data.TimeStep = obj.frameT;
                sptobj.saveFile();
            end
                       
            if obj.inGui
                currSPTGUIs = findobj('Type','Figure','Tag','SPT.gui');
                sptobj.gui();
                obj.appendOpenFigs(setdiff(findobj('Type','Figure','Tag','SPT.gui'),currSPTGUIs));
            end
        end
        
        
        
        %% RPT Analysis
        function rpt = trackRPT(obj, roi_in, defaultParams, forceNewRPT)
            % Open a new or existing RPT object for this roi
            % [in]
            %  roi_in - This is a roi index or name of existing ROI.
            %  defaultParams - An RPT object or anouther .rpt file to use as default params
            %  forceNewRPT - [default:false] A boolean if we should force the creation of a new RPT even if there is an existing file
            % [out]
            %  rpt - An rpt object.
            if nargin<4
                forceNewRPT=false;
            end
            if nargin<3
                defaultParams=[];
            end
            if ischar(defaultParams) && exist(defaultParams,'file')
                defaultRPT = RPT(defaultParams);
                defaultParams = defaultRPT.getParamStruct();
            elseif isa(defaultParams,'RPT')
                defaultParams = defaultParams.getParamStruct();
            end
            obj.assertInitialized();
            if nargin<2
                roi_idx = 1;
            else
                [~, roi_idx] = obj.getROI(roi_in);
            end
            if isempty(roi_idx)
                error('SPData:trackRPT','Bad ROI specified: %s',roi_in);
            end
            obj.checkComplete(); % Check that the data entries are complete
            
            filePath = obj.getROIFiles(roi_in, 'RPT', true);
            if isempty(filePath) || forceNewRPT % Create new ROI
                if isempty(defaultParams) % Try to look for other .RPT files to copy parameters from
                    allRPT = cellmap(@(roi) obj.getROIFiles(roi, 'RPT'), 1:numel(obj.ROI));
                    allRPT = [allRPT{:}];
                    if ~isempty(allRPT)
                        defaultRPT = RPT(allRPT{1});
                        defaultParams = defaultRPT.getParamStruct();
                    end
                end
                if ~isempty(defaultParams) % Load defult params
                    rpt = RPT();
                    rpt.setPropertyDefaults(defaultParams); % First set default params
                    rpt.load(obj,roi_idx); % Now load spdata and ROI
                else % start fresh
                    rpt = RPT(obj,roi_idx); 
                end
            else
                rpt = RPT(filePath); % load existing RPT
            end
        end
        
        %% GUI Management
        function closeGUI(obj)
            obj.gainCalGuiFig=GUIBuilder.closeCloseableHandles(obj.gainCalGuiFig);
            closeGUI@GUIBuilder(obj)
        end

        
        %% ROI Management
        function roi_size=getROIsize(obj, roi_in)
            % [out]
            %  roi_size: [Y X T]
            roi=obj.getROI(roi_in);
            roi_size=[roi(4)-roi(3)+1, roi(2)-roi(1)+1, roi(6)-roi(5)+1];
        end
        
        function [roi, roi_idx]=getROI(obj, in)
            %
            % Get an roi based on user input.  Allows multiple ways of selecting a ROI.  Monst methods
            % should call this method when the user input is an ROI specification.  This allows the user
            % many ways to specify ROI.
            %
            % in can be:
            %  - empty - Returns entire frame area as ROI [1 sizeX 1 sizeY]
            %  - roi-array format: [xmin xmax ymin ymax]
            %  - roi index
            %  - roi name as string
            obj.assertInitialized();
            roi = in;
            roi_idx = [];
            if nargin==1 || isempty(in) 
                roi=[1 obj.sizeX 1 obj.sizeY obj.globalTBounds];
                return
            elseif ischar(in) % slect based on roi name
                roi_idx = find(strcmpi(in,obj.ROIname),1,'first');
                if isempty(roi_idx)
                    error('SPData:getROI','Unknown ROI name: "%s"',in);
                end
                roi = obj.ROI{roi_idx};
            elseif isscalar(in) && round(in)==in %integer
                if in<=0 || in>length(obj.ROI)
                    error('SPData:getROI','Index out of bounds: %i',in);
                end
                roi_idx = in;
                roi = obj.ROI{roi_idx};
            end
            roi = obj.checkROI(roi);
        end
        
        function setGlobalTBounds(obj,tB)
            obj.assertInitialized();
            if length(tB)==2 && tB(1)>=1 && tB(2)<=obj.nFrames && tB(2)>tB(1)
                if length(obj.globalTBounds)==2 && all(obj.globalTBounds==tB)
                    return %No changes made
                end
                obj.globalTBounds=tB;
                for i=1:length(obj.ROI)
                    obj.ROI{i}(5)=max(tB(1),obj.ROI{i}(5));
                    obj.ROI{i}(6)=max(tB(2),obj.ROI{i}(6));
                end
                obj.updateSumImage();
                obj.dirty=true;
            else
                error('SPData:setGlobablTBounds', 'Bad settings for global time bounds: %s',arr2str(tB));
            end          
        end
        
        %% Dependent property get methods
        function val=get.nFrames(obj)
            %number of time points
            if ~isempty(obj.size);  val=obj.size(3);
            else val=[]; end
        end

        function val=get.sizeX(obj)
            %%number of pixels in the X (up-down) dimension
            if ~isempty(obj.size);  val=obj.size(1);
            else val=[]; end
        end
        
        function val=get.sizeY(obj)
            %number of pixels in the Y (left-right) dimension
            if ~isempty(obj.size);  val=obj.size(2);
            else val=[]; end
        end

        function val=get.frameSize(obj)
            % frame size is [X Y] 
            if ~isempty(obj.size);  val=obj.size(1:2);
            else val=[]; end
        end

        function val=get.nFramePixels(obj)
            % #pixels/frame 
            if ~isempty(obj.size);  val=prod(obj.size(1:2));
            else val=[]; end
        end

        function val=get.nExperimentPixels(obj)
            % #pixels for entire experiment 
            if ~isempty(obj.size);   val=prod(obj.size);
            else val=[]; end
        end

        function val=get.experimentT(obj)
            %time(s) to aquire entire experiment [X Y T] image
            if ~isempty(obj.size) && ~isempty(obj.frameT)
                val=obj.frameT*obj.nFrames;
            else
                val=[];
            end
        end

        function fname=get.rawDataPath(obj)
            if isempty(obj.workingDir) || isempty(obj.Paths.rawData)
                fname=[];
            else
                fname=fullfile(obj.workingDir, obj.Paths.rawData);
            end
        end

        function fmt=get.rawDataFormat(obj)
            p=obj.rawDataPath;
            if isempty(p)
                fmt=[];
            else
                [~,~,fmt]=fileparts(p);
            end
        end

        function val=get.calibrated(obj)
            % True if the data set can be calibated.  Should be overloaded by subclasses
            % w/different camera model.
            val= ~isempty(obj.CCDBackground) && ~isempty(obj.CCDGain);
        end

        function val=get.sequence(obj)
            val=dip_image(obj.getRawFrames());
        end

    end % Public methods 

    methods (Access=protected)
        function roi=checkROI(obj, roi_in, roi_name)
            if isempty(obj.size)
                error('SPData:checkROI', 'size not set.  Cannot add ROI.');
            end
            if nargin>2 && any(cellfun(@(s) strcmp(roi_name,s),obj.ROIname))
                error('SPData:checkROI', 'Already have ROI name "%s"',roi_name);
            end
            if isempty(obj.globalTBounds)
                tBounds=[1 obj.size(3)];
            else
                tBounds=obj.globalTBounds;
            end
            roi=roi_in(:)'; %make a row vector
            if length(roi)==4
                roi=[roi tBounds];
            end
            roi(1)=max(1,roi(1));
            roi(2)=min(obj.sizeX,roi(2));
            roi(3)=max(1,roi(3));
            roi(4)=min(obj.sizeY,roi(4));
            roi(5)=max(tBounds(1),roi(5));
            roi(6)=min(tBounds(2),roi(6));
            if length(roi)~=6 || any(round(roi)~=roi)
                error('SPData:checkROI','Bad ROI format');                 
            elseif roi(2)<roi(1)
                error('SPData:checkROI','Bad ROI X order');
            elseif roi(4)<roi(3)
                error('SPData:checkROI','Bad ROI Y order');
            elseif roi(6)<roi(5)
                error('SPData:checkROI','Bad ROI T order');
            end
        end
               
        function checkComplete(obj)
            %Check all necessary properties are set for tracking and
            %analysis
            obj.assertInitialized();
            if ~obj.calibrated
                error('SPData:checkComplete','Data is not yet calibrated');
            elseif isempty(obj.pixelSize)
                error('SPData:checkComplete','PixelSize is not set.');
            elseif isempty(obj.psf)
                error('SPData:checkComplete','Point Spread function size is not set.');
            elseif isempty(obj.frameT)
                error('SPData:checkComplete','Frame Time is not set.');
            elseif obj.dirty
                error('SPData:checkComplete','Object is diry.  Try saving first.');
            end
        end
        
        function saveRawData(obj, raw_data_path)
            %Internal helper for saving a new raw Data file.  This should not be called directly, it is
            % used by higher level public methods like saveDataAs 
            [~,~,format] = fileparts(raw_data_path);
            switch format
                case '.mat'
                    %avoid globalTimeBounds slicing and get whole dataset
                    sequence=obj.getRawFrames(); %#ok<NASGU>
                    if isfield(obj.AcquisitionParams,'paraCCD')
                        %Check for and copy parameters from the Mat capture format.  These
                        %would only exist if we loaded already from a previous copy.
                        paraCCD=obj.AcquisitionParams.paraCCD; %#ok<NASGU>
                        paraOptics=obj.AcquisitionParams.paraOptics; %#ok<NASGU>
                        paraMisc=obj.AcquisitionParams.paraMisc; %#ok<NASGU>
                        save(raw_data_path, 'sequence','paraCCD','paraOptics','paraMisc','-v7.3');
                    else
                        save(raw_data_path, 'sequence','-v7.3');
                    end
                otherwise
                    error('SPData:saveRawData', 'Error cannot save raw data type "%s"',format);
            end
            obj.Paths.rawData=relativepath(obj.workingDir,raw_data_path);
            obj.dirty=true;
        end

        function updateSumImage(obj)
            obj.assertInitialized();
            if obj.calibrated
                obj.sumImage=sumImage2D(obj.getFrames());
            else
                obj.sumImage=sumImage2D(obj.getRawFrames());
            end
        end

        function load_raw_frames(obj)
            % Helper to internaly load the frames
            % Normally only called once when obj.raw_frames_loaded=false;
            % Calling agian will reload.
            if ~exist(obj.rawDataPath,'file')
                return
            end
            obj.updateWaitbar(0,'Loading raw data');
            try
                obj.raw_frames_loaded=false;
                obj.raw_frames_=obj.loadRawDipImage(obj.rawDataPath);
                obj.updateWaitbar(0.9);
                if isempty(obj.sumImage) % This is in case we have a saved sumImage and then get the rawFrames
                    obj.sumImage=sumImage2D(obj.raw_frames_);
                end
            catch err
                obj.updateWaitbar(1);
                throw(err);
            end
            obj.raw_frames_loaded=true;
            obj.updateWaitbar(1);
        end

        function load_frames(obj)
            % Helper to internaly load the gain-calibrated frames
            % Normally only called once when obj.frames_loaded=false;
            % Calling agian will cause re-calibration.
            if ~obj.calibrated
                return
            end
            obj.frames_loaded=false;
            obj.recalibrateFrames();
            if isempty(obj.frames_) 
                error('SPData:load_frames','Unable to load frames.');
            end
            obj.frames_loaded=true;
        end

        function load_CCDGainCalFrames(obj)
            % Helper to actualy load the frames and store the data in the private property
            filepath=obj.getFilePath('CCDGainCal');
            if ~exist(filepath,'file')
                return
            end
            obj.CCDGainCal_frames_loaded=false;
            obj.CCDGainCal_frames_=obj.loadRawDipImage(filepath);
            obj.CCDGainCal_frames_loaded=true;
        end

        function load_CCDBackgroundFrames(obj)
            % Helper to actualy load the frames and store the data in the private property
            filepath=obj.getFilePath('CCDBackground');
            if ~exist(filepath,'file')
                return
            end
            obj.CCDBackground_frames_loaded=false;
            obj.CCDBackground_frames_=obj.loadRawDipImage(filepath);
            obj.CCDBackground_frames_loaded=true;
        end      

        function loadSPT(obj, sptobj)
            % Load from a .spt file
            if isempty(sptobj.DataFile)
                error('SPData:loadSPT','SPT object has no DataFile property set');
            end
            if ~exist(sptobj.DataFile, 'file')
                warning('SPData:loadSPT','Unable to find raw data file: "%s"', sptobj.DataFile);
                title=sprintf('Select raw data file: %s',sptobj.DataFile);
                sptobj.DataFile=uigetfile('*.mat',title,sptobj.DataPath);
                if ~exist(sptobj.DataFile, 'file')
                    error('SPData:loadSPT','Unable to load raw data file: "%s"', sptobj.DataFile);
                end
            end
            [datapath, datafile, dataext]=fileparts(sptobj.DataFile);
            if ~strcmp(dataext,'.mat')
                error('SPData:loadSPT','SPT object data file is not a .mat');
            end
            %Now we are really overwriting object, so save previous properties.
            if obj.initialized
                obj.savePreservedProperties();
                pp=obj.preservedProperties;
                obj.reset();
                obj.preservedProperties=pp;
            end

            obj.workingDir=collapsepath(datapath);
            obj.Paths.rawData=[datafile '.mat'];
            obj.Paths.saveFile=[datafile SPData.saveFileExt];
            obj.psf=sptobj.ParamsGeneral.Psf;
            obj.pixelSize=sptobj.ParamsGeneral.PixelSize;
            obj.frameT=sptobj.ParamsGeneral.TimeStep;
            obj.CCDBackground=sptobj.ParamsGeneral.CCDOffset;
            obj.CCDGain=sptobj.ParamsGeneral.Gain;
            if ~isempty(sptobj.ParamsGeneral.Roi)
                if ~isempty(sptobj.ParamsGeneral.Frames)
                    fs=sptobj.ParamsGeneral.Frames;
                    obj.addROI([sptobj.ParamsGeneral.Roi+1, min(fs), max(fs)]);
                else
                    obj.addROI([sptobj.ParamsGeneral.Roi+1, 1, obj.nFrames]);
                end
            end
            obj.loadRawDataParameters();
            obj.initialized=true;
            obj.dirty=true;
        end     
    
        function loadMatCapture(obj, capturefile)
            % This allows loading from ".mat" files which are home-brew formats produced by the LidkeLab
            % instrumentation classes.  The "Old"-style .mat files (version 1) were produced by the 
            % LaserControler class.
            % They store the image as a "sequence" variable in dip_image format.  
            % They also include a 'paraCCD', and 'paraMisc' structures with some limited information.  
            % Most of the usefull feilds are however not filled in.  
            %
            % The "New" format (version 2) is saved by SPTSRcollect.  It has a 'Params' value that stores
            % information.  Importantly the sequence variable is stored permuted from the "standard" 
            % orientation and as a unit16 value.
            %
            % Finally there is a 'raw' .mat format, where there is just a 'sequence' variable.  
            % These are created when converting .tif to .mat format.  Since there are no instrumentation
            % objects availible to save there are no auxillary 'paraMisc' or 'Params' variables saved 
            % along with the sequence.   Importantly the image orientation is correct in that it 
            % matches version 1 orientation.
            [filepath, filename, ext]=fileparts(capturefile);
            obj.updateWaitbar(0,sprintf('Loading Raw Data "%s"',[filename ext]));            
            capture=load(capturefile);
            if ~isfield(capture,'sequence')
                obj.updateWaitbar(1); %close waitbar
                error('SPData:loadMatCapture','Capture .mat file has no sequence field.');
            end
            obj.updateWaitbar(0.9,'Loading parameters');
            %Now we need to 'detect' what format the .mat file is and load appropriately.
            %Need to confer with Keith and Sheng to understand what versions actually exist in data
            if isfield(capture,'paraCCD')
                %This is a version 1 .mat produced by LaserControl??
                obj.loadMatCapture_format1(capture);
            elseif isfield(capture,'Params')
                %This is a version 2 .mat produced by SPTRSRcollect
                obj.loadMatCapture_format2(capture);
            elseif isfield(capture,'Version') && capture.Version==3
                %This is a version 3 .mat produced by SPTRSRcollect                
                obj.loadMatCapture_format3(capture);
            else
                %This is a 'raw' .mat with just a 'sequence' variable
                obj.loadSequence(capture.sequence);  
            end
            
            obj.workingDir = collapsepath(filepath);
            obj.Paths.rawData = [filename '.mat'];
            obj.Paths.saveFile = [filename SPData.saveFileExt];

            obj.updateWaitbar(1);
            obj.dirty=true;
        end

        function loadMatCapture_format1(obj, capture)
            % This is called by loadMatCapture to read "version1" formatted files.
            obj.loadSequence(capture.sequence); 
            if isfield(capture,'paraCCD')
                if ~isempty(capture.paraCCD.PixelSize)
                    obj.pixelSize=capture.paraCCD.PixelSize;
                end
                obj.frameT=capture.paraCCD.FrameTime;
                obj.AcquisitionParams.paraCCD=capture.paraCCD;
                obj.AcquisitionParams.paraOptics=capture.paraOptics;
            end
            if isfield(capture,'paraMisc')
                obj.AcquisitionParams.paraMisc=capture.paraMisc;
            end
        end
        
        function loadMatCapture_format2(obj, capture)
            % This is called by loadMatCapture to read "version2" formatted files.
            % A distinguishing feature is that the data is saved in a permuted format so we need to unpermute it
            % to have it work correctly.
            obj.loadSequence( permute(capture.sequence, [2 1 3])); %            
            if isfield(capture,'Params')
                params = capture.Params;
                obj.AcquisitionParams=params; %Save this internally
                if isfield(params,'CameraObj')
                    camera = params.CameraObj;
                    if isfield(camera,'CameraParameters')
                        cameraParams = camera.CameraParameters;
                        if isfield(cameraParams, 'ExposureTime')
                            obj.frameT = cameraParams.ExposureTime;
                        end
                    end
                end
            end
        end
        
        function loadMatCapture_format3(obj, capture)
            % This is called by loadMatCapture to read "version3" formatted files.
            obj.loadSequence(capture.sequence);
            obj.AcquisitionParams=capture.AcquisitionParams;
            obj.frameT = capture.AcquisitionParams.Camera.FrameTime;
        end
            
        function loadTifCapture(obj, capturefile)
            % Prerconditions:
            %  capturefile exists and is .tif file.
            [filepath, filename, ext]=fileparts(capturefile);
            obj.updateWaitbar(0,sprintf('Loading Raw Data "%s"',[filename ext]));
            sequence=readmultitif(capturefile,obj.waitbarH);
            if isempty(sequence)
                obj.updateWaitbar(1); %close waitbar
                error('SPData:loadTifCapture','Unable to readmultitif() on file "%s"',capturefile);
            end
            obj.updateWaitbar(0.99,'Loading parameters');
            obj.loadSequence(sequence); %Do the common sequence loading subroutine
            obj.workingDir=collapsepath(filepath);
            obj.Paths.rawData=[filename '.tif'];
            obj.Paths.saveFile=[filename SPData.saveFileExt];
            obj.updateWaitbar(1);
            obj.dirty=true;
        end
        
        function loadDipImage(obj, im, spdatapath)
            % load from a dipimage.
            % spdatapath.  The full path to the new .spdata file
            [filepath, filename, ext]=fileparts(spdatapath);
            if isempty(filepath) || isempty(filename) || (~isempty(ext) && ~strcmp(ext,SPData.saveFileExt))
                error('SPData:loadDipImage','Invalid .spdata path "%s"',spdatapath);
            end
            obj.loadSequence(im);
            obj.workingDir=collapsepath(filepath);
            obj.Paths.rawData=[filename '.mat'];
            obj.Paths.saveFile=[filename SPData.saveFileExt];
            obj.dataDirty=true; %Mark fact that data is not saved yet.
        end

        function loadSequence(obj, sequence)
            %Load given a sequence of raw data not already associated with a
            %.spdata file.  This is the inner core of the loadDipImage, loadMatCapture and
            %loadTifCapture routines that are all eventaully getting a dip_image sequence
            %in some way.  The common details go here.
            if isempty(sequence)
                error('SPData:loadSequence','Data file contained an empty sequence');
            end
            if obj.initialized
                %Save any previous paramters to use as defaults.
                obj.savePreservedProperties();
                pp=obj.preservedProperties;
                obj.resetObject();
                obj.preservedProperties=pp;
            end
            %Save the raw frames so we don not need to reload
            if isa(sequence,'dip_image')
                %Size [X Y T] is for a dip_image
                obj.size=size(sequence); %#ok<*CPROP>
            else
                %Size [Y X T] is for a matlab array *IMORTANT*
                obj.size=size(sequence);
                obj.size(1:2)=[obj.size(2) obj.size(1)];
            end
            obj.raw_frames_=uint16(sequence);
            obj.raw_frames_loaded=true;
            obj.sumImage=sumImage2D(obj.raw_frames_);
            if ~isfield(obj.preservedProperties,'globalTBounds') || ~isempty(obj.preservedProperties.globalTBounds)
                obj.globalTBounds = [1, obj.size(3)];
            end
            obj.filterValidROI();
            obj.initialized=true;
            obj.dirty=true;
        end

        function frames=loadRawDipImage(obj,filepath)
            % A helper to load a raw file as a dipimage.  Used to load rawData and
            % calibration data, etc.  Can load .mat or .tif.  Uses the fixed
            % readmultitiff which avoids the bugs in readtimeseries
            %
            % [OUT]
            %   frames - gaurenteed to be non-empty on return or we throw
            %   an execpetion
            if ~exist(filepath,'file')
                error('SPData:loadRawDipImage','Unable to open file "%s"',filepath);
            end
            [~,~,ext]=fileparts(filepath);
            switch ext
                case '.mat'
                    d=load(filepath);
                    sequence=d.sequence;
                    if isempty(sequence)
                        error('SPData:loadRawDipImage','sequence is empty.');
                    end
                    if isa(sequence,'uint16')
                        frames=sequence;
                    elseif isa(sequence','dip_image')
                        if ~strcmp(datatype(sequence),'uint16')
                            error('SPData:loadRawDipImage','sequence is not uint16.  Unable to load.');
                        end
                        frames=uint16(sequence);
                    end
                    if isfield(d,'Params')
                        %This is a version 2 .mat produced by SPTRSRcollect must flip frames to load correctly
                        frames = permute(frames, [2,1,3]);
                    end
                case '.tif'
                    frames=readmultitif(filepath,obj.waitbarH);
                    if isempty(frames)
                        error('SPData:loadRawDipImage','Multi-page tif load is empty.');
                    end
                    if ~isa(frames,'uint16')
                        error('SPData:loadRawDipImage','sequence is not uint16.  Unable to load.');
                    end
                otherwise
                    error('SPData:loadRawDipImage','Unknown how to load ext "%s".',ext);
            end
            if isempty(frames)
                error('SPData:loadRawDipImage','Got an empty image');
            end
        end

        function name=makeROITitleString(obj,roi_in)
            %Make a nice title for the window
            roi = obj.getROI(roi_in);
            if isscalar(roi_in)
                name = sprintf('<%s>',obj.ROIname{roi_in});
            else
                name = 'ROI';
            end
            roi_cell = num2cell(roi);
            name = sprintf('%s {x:[%i,%i],y:[%i,%i],t:[%i,%i]}',name, roi_cell{:});
        end

        %% Abstract methods inherited from Pickle
        function val = getProtectedProperty(obj, name)
            %This is necessary for Pickle functionality to be able to access subclass protected variables
            val = obj.(name);
        end

        function modifyProtectedProperty(obj, name, newval)
            %This is necessary for Pickel functionality to be able to change subclass variables
            obj.(name)=newval;
        end
    end % protected methods

    methods (Static=true)
        function obj=loadobj(propstruct)
            %Cleanup the propstruct incase we run into old variables

            %Property 'Paths' was previously 'Filenames'
            if isfield(propstruct,'Filenames')
                propstruct.Paths=propstruct.Filenames;
                propstruct=rmfield(propstruct,'Filenames');
            end
            %Paths names changed from earlier versions, make sure to fix
            %these up.
            if isfield(propstruct,'Paths')
                paths=propstruct.Paths;
                if isfield(paths,'spdata')
                    paths.saveFile=paths.spdata;
                    paths=rmfield(paths,'spdata');
                end
                blankspd=SPData();
                newpaths=blankspd.Paths();
                newfields=fieldnames(newpaths);
                for i=1:length(newfields)
                    f=newfields{i};
                    if isfield(paths,f)
                        newpaths.(f)=paths.(f);
                    end
                end
                propstruct.Paths=newpaths;
            end
            %New property globalTBounds
            if ~isfield(propstruct,'globalTBounds')
                propstruct.globalTBounds=[1 propstruct.size(3)];
            end

            obj = Pickle.loadobj(propstruct); %Let pickle do the rest
        end

        function new_files=batchProcess(datapath, datafile_patterns, defaultParams, overwriteFlag) 
            % Make and save a .spdata for each of the data files matching a pattern
            %  found in the paths in datapath array. 
            % [IN]
            %  dataPath - A path to a directory where files are to be converted
            %  datafile_patterns - A cell-array of file patterns that can use wildcard '*'
            %                      to search for (a pattern can also just be a filename with no wildacard)
            %                      Ex: {'2015-01-01*.mat', '2015-01-03-condition2.mat'}
            %  defaultParams - An example SPData object or filename or preservedParams struct
            %  overwriteFlag - [optional] Integer 0=Do not overwrite; [Default] 
            %                                     1=Warn with selection dialog before overwrite;
            %                                     2=Force overwrite (caution!);
            % [OUT]
            %  new_files - Cell array of full-paths to all .spdata files corresponding to the given 
            %             data file patterns
            %             (This list inculdes the default and non-overwitten files, as we assum you will
            %             want to process all these files similarly in the next phase of batch
            %             processing).
            if nargin<4
                overwriteFlag = 0;
            end
            if ischar(defaultParams) && exist(defaultParams,'file')
                % We got a filename so, Open file make object and retrieve the properties struct
                defaultFilePath = defaultParams;
                spd=SPData(defaultFilePath);
                defaultParams = spd.getParamStruct();
            elseif isa(defaultParams,'SPData')
                % We got an object so call getParamsStruct to get it as a properties struct
                defaultFilePath = defaultParams.saveFilePath();
                defaultParams = defaultParams.getParamStruct();
            elseif isstruct(defaultParams)
                defaultFilePath = [];
            else
                error('SPData:batchProcess','Unknown default parameters format');
            end
            %Gather full paths for all files matching any of the patterns.
            filenames = cellmap(@(p) Pickle.listExistingFileNames(datapath,p), makecell(datafile_patterns));
            filenames = [filenames{:}];
            new_files = Pickle.changeFileExtensions(filenames, SPData.saveFileExt);
            try
                overwrite = (overwriteFlag==2) || SPData.confirmOverwriteDialog(new_files); % Determine if we need to overwrite
            catch EX
                switch EX.identifier
                    case 'Pickle:OverwriteCancel'
                        new_files = {};
                        return;
                    otherwise
                        rethrow(EX);
                end
            end
            Nfiles = numel(new_files);
            Nprocessed = 0;
            Nerror = 0;
            Nexisting = 0;
            H=waitbar(0,'Batch Processing ...');
            for n=1:Nfiles % Create and save the .spdata file for this data
                fprintf('\n*** [%i/%i] Batch Processing Datafile:\n  --->%s\n',n,Nfiles,filenames{n});
                waitbar(n/(Nfiles+1),H,sprintf('Batch Processing Datafile [%i/%i] ...', n,Nfiles));
                if exist(new_files{n},'file') && (strcmp(new_files{n},defaultFilePath) || ~overwrite)
                    Nexisting = Nexisting + 1;
                    continue %Don't overwrite the defaultParams file or any file if overwrite flag is off
                end 
                %Here is the actual creation and savin of the new .spdata file
                try
                    spd = SPData();
                    spd.disableWaitbar = true;
                    spd.setPreservedProperties(defaultParams); % Set the internal preserved properties based on the default params
                    spd.load(filenames{n}); % Load in the data with any known parameters
                catch err
                    fprintf('>>>(oops)<<< Batch Processing Error: "%s"\n',spd.saveFilePath);
                    disp(getReport(err))
                    new_files{n} = [];
                    Nerror = Nerror+1;
                    continue;
                end
                if ~strcmp('.mat',spd.rawDataFormat)
                    spd.convertToMat();
                end
                spd.save();
                Nprocessed = Nprocessed +1;
                new_files{n} = spd.saveFilePath;
            end
            close(H);
            fprintf('\n*** Batch Processing Complete. [Num Total: %i | Num Processed: %i | Num Existing: %i | Num Error: %i]\n',Nfiles,Nprocessed,Nexisting,Nerror);
        end
    end %public static methods
end % classdef
