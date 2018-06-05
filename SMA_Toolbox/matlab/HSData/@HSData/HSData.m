% HSData.m
% Mark J. Olah (mjo@cs.unm.edu)
% 07/21/14
%

classdef HSData < BaseData
    %HSData This is a streamlined replacement for the HSIData class, that focuses on
    % loading data for analysis and post processing.
    %
    % This class can load raw uint16 data from the .ics files, and also load
    % frames corrected for gain, background and baseline offset in single format.
    % Finally we can render frames in various formats including RGB dipimage format.
    %
    % All frames manipulated by this class are stored in column-major format wich are
    % indexed as [L Y X T]
    %
    %
    % Loading and Saving:
    %
    % The design of the class is simplified by assuming that all relevent files for this
    % experiment live in the same directory given by obj.workingDir;.  It is the users job
    % to put together the files.  An .hsd can be loaded without the calibration or bgMean
    % files because that data is stored in properties and saved like all other normal
    % properties.
    %
    % HSData objects can be initialized with a .hsi file from the old hsiData class, or
    % with an hsiData object.  Or they can be initialized with a .hsd data file from
    % this class or anouther HSData object.  These objects or files are given through the
    % constructor or the load command.


    properties (Constant=true)
        %% Inherited from Pickle
        saveFileExt='.hsdata'; % The extension used for saved files
        SaveableDataFormats={ '*.hsdata','HSData object files (.hsdata)'};
        LoadableDataFormats={'*.hsi;*.hsdata','All Loadable Sources (.hsi, .hsdata)';...
                             '*.hsdata','HSData object files (.hsdata)';...
                             '*.hsi','HSIData object files (.hsi)'};       
        
        RawDataFormats={'.hsi', '.hsdata'};

        NBaselineRows=3; % The number of columns blanked out with a high-pass filter.        
    end

    properties        
        %Fundemental aqcuisition parameters
        size;% [L Y X T] (size of actual frames after baseline correction)
        pixelSize; %spatial pixel size (um)
        psf; %[x y] in pixels
        sensorT; %time to aquire 1 sensor [L Y] image (s)
        rawLambda; %Wavelength of each row along lambda dimension of raw images (nm) (includes the baseline correction rows)

        %Gain calibration 
        CCDBackground; %mean dark background image Dims: [L, Y]
        CCDGain; % e-/ADU


        %Data visualization
        % Lambda mask
        currentLambdaMaskIdx;
        lambdaMask={}; %logical 1xsizeL 1=keep wavelength; 0=remove wavelength
        lambdaMaskName={}; %Name for masks

        colorMap;   % double [sizeL,3] maps lambda to RGB colors
        

        Paths=struct(... %All filenames are relative paths from the workingDir
             'saveFile',[],... %filename(w/extension) of the .hsdata file where this will save to
             'hsi',[],...  %filename(w/extension) where the hsi file (if any)
             'icsBaseName',[],... %relative path to location of .ics image files (w/o number or extension).
             'CCDBackground',[],... %filename(w/extension) where the bgMean information was read
             'CCDGainCal',[],... %filename(w/extension) where the gain and calibration information was read
             'HSRPT','HSRPT'...
             );

        %Reference aquisition information, extracted from .hsi file or object.
        %We don't make use of this directly, but wish to keep it recorded.
        AcquisitionParams=struct(...
            'acqDateTime', [],...       %Vector datetime of sequence acquisition
            'gainCalDateTime', [],...   %Vector datetime of gain calibration
            'wvCalDateTime', [],...     %Vector datetime of wavelength calibration
            'laserPower', [],...        %Power or intensity?? mW kW/cm^2 (?)
            'sensorTemperature', [],... %Sensor temperature deg C
            'cameraData', []...         %Camera struct
            );

        sumImageGray; % A grayscale image with lambda collapsed out based on the globalTBound
        sumImageRGB;  %
        spectra;   % size=[1, sizeL] intensity at each spectral value over the whole image
        spectralSeries; % size=[nFrames, sizeL] intensity summed over X,Y for each T,L pair
    end

    properties (Hidden=true)
        version=1; %For future file format version changes
        AnalysisPathName='HSAnalysis'; %subdirectory where analysis files are kept
    end

    properties (Dependent=true)
        nFrames; %number of frames after taking into account the globalTBounds.
        sizeX; %number of pixels in the X (scaned) dimension
        sizeY; %number of pixels in the Y (non-scanned) dimension
        sizeL; %number of pixels in the L dimension after baseline correction
        lambda; %Wavelength of each row along lambda dimension of frame (nm) (This is exactly lambdaRaw but without last NBaselineRows which are low-pass filtered).
        sensorSize; % [L Y]
        frameSize; % [L Y X]
        physicalSize; % [X Y] (micron)
        NSensorPixels; % sizeL*sizeY
        NFramePixels; % sizeL*sizeY*sizeX
        NExperimentPixels; % sizeL*sizeY*sizeX*sizeT
        
        %Timing information
        frameT;  %time to aquire 1 frame [L Y X] image (s)
        experimentT; %time to aquire entire experiment [L Y X T] (s)        

        icsBasePath;  %The full path without number or extension for the .ics files
        analysisPath; % Full path to analysis base for this data file (all files but .hsdata and rawData should be in here)
    end



    %These dependent properties are not shown in a property list, so it keeps the presentation
    %cleaner, but they are still public and should be used to get access to the files
    properties (Dependent=true, Hidden=true)
       %Raw sizes: These sizes include the NBaselineRows that are dropped from
        %the L dimension to correct for baseline offset.
        rawSize; %[L Y X T] The size of raw data for entire experiment
        rawSizeL; %number of pixels in the L dimension before baseline correction
        rawSensorSize; % [L Y]

        sequence;%The raw data as a dip_image object.  This conceptually mirrors the 'sequence' name in the .mat file.          
    end
    

    %% Private properties
    properties (Access=private, Constant=true)
        default_psf= [0.8 1.14]; % [psf_x psf_y]
        ifaceHandle= @HSData_Iface; %A handle to the iface class constructor
    end

    properties (Access=protected, Transient=true)      
        objectHandle=0; %Handle to HSData allocated c++ class in memory.  0=none
    end

    methods
        function obj=HSData( in )
            % This constructor calls the obj.load method with the arguments
            % given.  Empty input arguments give an empty object, otherwise the load
            % method determines the valid input arguments
            if nargin>0
                obj.load(in);
            end
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(obj)
            if obj.initialized && obj.objectHandle
                obj.call('@delete');
            end
        end

        function load(obj, in)
            % [in]
            %  (1) hsiData Object
            %  (2) .hsi filepath
            %  (3) .hsdata filepath
            if isa(in,'hsiData')
                % Case (1)
                obj.loadHSI(in.AcqParams, in.DataFile);
            elseif ischar(in)
                filepath=in;
                if ~exist(filepath,'file')
                    error('HSData:FileNotFound', 'unable to open file: %s', filepath);
                end
                [~, ~, ext]=fileparts(filepath);
                switch ext
                    case '.hsi' % Case (2)
                        hsi=load(filepath,'-mat');
                        obj.loadHSI(hsi.params, filepath);
                    case HSData.saveFileExt % Case (3)
                        obj.loadHSData(filepath);
                    otherwise
                        warning('HSData:HSData',['Unknown file type ' ext]);
                end
            else
                warning('HSData:HSData',['Unexpected input type ' class(in)]);
            end
        end

        function loadBackground(obj,filename)
            % The bgMean is stored in the save .hsd files and does not need to be reloaded
            % typically after it is read in the first time from a .hsi file load
            % [in] filename - string of full path to new background file.  If filename is
            %                 omitted, use current filename
            if nargin>1 %Change Filename.bg
                [path,file,ext]=fileparts(filename);
                if isempty(path)
                    path=obj.workingDir;
                elseif ~strcmp(path,obj.workingDir)
                    error('HSData:reloadBackground',['Cannot load file outside of dataDir: ' obj.workingDir]);
                end
                if ~exist(fullfile(path,[file ext]),'file')
                    error('HSData:reloadBackground',['File does not exist: ' filename]);
                end
                obj.Paths.CCDBackground=[file ext];
            end
            S=load(obj.getFilePath('CCDBackground'),'-mat');
            obj.CCDBackground=S.bgMean;
        end

       function loadCalibration(obj,filename)
            % The gain is stored in the save .hsd files and does not need to be reloaded
            % typically after it is read in the first time from a .hsi file load
            % [in] filename - string of full path to new calibration file.  If filename is
            %                 omitted, use current filename
            if nargin>1 %Change Filename.cal
                [path,file,ext]=fileparts(filename);
                if isempty(path)
                    path=obj.workingDir;
                elseif ~strcmp(path,obj.workingDir)
                    error('HSData:reloadCalibration',['Cannot load file outside of dataDir: ' obj.workingDir]);
                end
                if ~exist(fullfile(path,[file ext]),'file')
                    error('HSData:reloadCalibration',['File does not exist: ' filename]);
                end
                obj.Paths.CCDGainCal=[file ext];
            end
            S=load(obj.getFilePath('CCDGainCal'),'-mat');
            obj.CCDGain=single(S.gain(2));
            obj.AcquisitionParams.gainCalDateTime=S.time;
       end

        %% ROI Management
        function roi_size=getROIsize(obj, roi_in)
            % [out]
            %  roi_size: [L Y X T]
            roi=obj.getROI(roi_in);
            roi_size=[roi(6)-roi(5)+1, roi(4)-roi(3)+1, roi(2)-roi(1)+1, roi(8)-roi(7)+1];
        end
        
        function roi=getROI(obj, in)
            %
            % Get an roi based on user input.  Allows multiple ways of selecting a ROI
            %
            obj.assertInitialized();
            roi=in;
            if isempty(in)
                roi=[1 obj.sizeX 1 obj.sizeY 1 obj.sizeL obj.globalTBounds];
                return
            elseif isscalar(in) && round(in)==in %integer
                if in<=0 || in>length(obj.ROI)
                    error('HSData:getROI','Index out of bounds: %i',in);
                end
                roi=obj.ROI{in};
            end
            roi=obj.checkROI(roi);
        end
        
        %% Lambda mask methods
        function mask=getLambdaMask(obj)
            mask=obj.lambdaMask{obj.currentLambdaMaskIdx};
        end
        
        function setLambdaMask(obj, varargin)
            %
            % Form 1: No arguments.
            %   If obj.lambdaMask is empty, create a new default full-spectrum mask.
            %   Else set currentLambdaMaskaskIdex=1
            % Form 2:
            %   [in] mask_idx: integer index of lambdaMask to set currentLambdaMaskIdx
            % Form 3:
            %  [in] mask: 1xsizeL  true=remove wavelength false=keep wavelength 
            
            if nargin==1 || isempty(varargin{1})
                if isempty(obj.lambdaMask)
                    obj.currentLambdaMaskIdx=obj.addLambdaMask(false(1,obj.sizeL),'FullSpectrum'); % Make a default color mask                    
                else
                    obj.currentLambdaMaskIdx=1; % set default color mask
                end
            elseif isscalar(varargin{1})
                mask_idx=varargin{1};
                if mask_idx<1 || mask_idx>length(obj.lambdaMask)
                    error('HSData:setLambdaMask','Unable to set lambda mask to idx=%i',mask_idx);
                end
                obj.currentLambdaMaskIdx=mask_idx;
            else
                obj.currentLambdaIdx=obj.addLambdaMask(varargin{:});
            end
        end
            
        function idx=addLambdaMask(obj, mask, name)
            %  [in] mask: 1xsizeL  true=remove wavelength false=keep wavelength 
            % [out] mask: 1xsizeL logical is the new effective mask
            if nargin<3
                name=obj.nextUnusedName(obj.lambdaMaskName,'Mask%i');
            end
            mask=obj.checkLambdaMask(mask,name);
            obj.lambdaMask{end+1}=mask;
            obj.lambdaMaskName{end+1}=name;
            obj.dirty=true;
            idx=length(obj.lambdaMask);
        end
        
        function deleteLambdaMask(obj, mask_idx)
            if mask_idx>0 && mask_idx<=length(obj.lambdaMask)
                obj.lambdaMask(mask_idx)=[];
                obj.lambdaMaskName(mask_idx)=[];
                obj.dirty=true;
            end
        end
        
        function modifyLambdaMask(obj, mask_idx, new_mask, new_mask_name)
            if mask_idx>0 && mask_idx<=length(obj.lambdaMask)
                if strcmp(new_mask_name, obj.lambdaMaskName{mask_idx})
                    if all(new_mask==obj.lmabdaMask{mask_idx})
                        return %No changes
                    end
                    new_mask=obj.checkLambdaMask(new_mask);
                else
                    new_mask=obj.checkLambdaMask(new_mask,new_mask_name);
                end
                obj.lambdaMask{mask_idx}=new_mask;
                obj.lambdaMaskName{mask_idx}=new_mask_name;
                obj.dirty=true;
            end
        end

        function Lwv=lambdaPixelsToWavelength(obj, Lpx)
            %Given lambda in pixels return lambda in wavelength
            Lwv = interp1((1:obj.sizeL)-0.5, obj.lambda, Lpx, 'pchip');
        end

        function Lpx=lambdaWavelengthToPixels(obj, Lwv)
            %Given lambda in wavelength return lambda in pixels
            Lpx = interp1(obj.lambda,(1:obj.sizeL)-0.5, Lwv, 'pchip');
        end

        %% ColorMap
        function map=getColorMap(obj)
            map=obj.colorMap;
            map(obj.getLambdaMask(),:)=[0 0 0];
        end
        
        function setColorMap(obj, cmap)
            if nargin==1
                cmap=obj.hyperCM(obj.lambda, obj.getLambdaMask());
            end
            obj.colorMap=cmap;
            obj.dirty=true;
        end
        %% Tracking
        function rpt = trackHSRPT(obj, roi_idx)
            obj.assertInitialized();
            if nargin==1
                roi_idx = 1;
            end
            rpt = HSRPT(obj,roi_idx);
%             obj.appendOpenFigs(rpt.gui());
        end


        %% Frame access

        function imref = getImRef(obj, roi_in)
            if nargin==1; roi_in=[]; end
            roi=obj.getROI(roi_in);
            sizeX=roi(2)-roi(1)+1;
            sizeY=roi(4)-roi(3)+1;
            imref = imref2d([sizeY, sizeX], obj.pixelSize(1), obj.pixelSize(1));
        end

        function [frames,roi]=getFrames(obj, roi_in)
            % Get gain, background and baseline corrected frames as matlab single array.  
            % This uses the cached copy of the gain-corrected frames, so we only have to gain correct once.
            % In:
            %   roi_in - A 1x8 ROI or a index into the local ROI cell array
            % Out:
            %   frames - single array size: [L Y X T] 
            %   roi - 1x6 array: [xmin, xmax, ymin, ymax, Lmin, Lmax tmin, tmax]
            obj.assertInitialized()
            if nargin==1; roi_in=[]; end
            if ~obj.frames_loaded 
                obj.load_frames() %Try to load 
            end
            roi=obj.getROI(roi_in);
            frames=obj.frames_(roi(5):roi(6), roi(3):roi(4), roi(1):roi(2), roi(7):roi(8)); %L Y X T
        end
        
        function [frames,roi]=getFramesRGB(obj, roi_in)
            if nargin==1; roi_in=[]; end
            [frames,roi]=obj.getFrames(roi_in);
            if isempty(frames)
                return
            end
            %TODO
        end
        
        function [spectra,lambda,roi]=getSpectra(obj, roi_in)
            if nargin==1
                roi_in=[];
            end
            if isempty(roi_in) && ~isempty(obj.spectra) %Look for cached copy
                spectra=obj.spectra;
                roi=[1 obj.sizeX 1 obj.sizeY 1 obj.sizeL obj.globalTBounds];
                lambda=obj.lambda(roi(5):roi(6));
                return
            end
            [frames,roi]=obj.getFrames(roi_in);
            if isempty(roi) || isempty(frames) 
                error('HSData:getSpectra','Unable to load frames');
            end
            sz=size(frames);
            spectra=sum(reshape(frames,sz(1),prod(sz(2:end))),2);
            lambda=obj.lambda(roi(5):roi(6));
        end

        function [spectra,lambda,roi]=getSpectralSeries(obj, roi_in)
            if nargin==1
                roi_in=[];
            end
            if isempty(roi_in) && ~isempty(obj.spectralSeries) %Look for cached copy
                spectra=obj.spectralSeries;
                roi=[1 obj.sizeX 1 obj.sizeY 1 obj.sizeL obj.globalTBounds];
                lambda=obj.lambda(roi(5):roi(6));
                return
            end
                
                
            [frames,roi]=obj.getFrames(roi_in);
            if isempty(roi) || isempty(frames) 
                error('HSData:getSpectra','Unable to load frames');
            end
            sz=size(frames);
            spectra=squeeze(sum(reshape(frames,sz(1),prod(sz(2:3)),sz(4)),2));
            lambda=obj.lambda(roi(5):roi(6));
        end
        
        function [gray_im,roi]=getSumImageGray(obj,roi_in)
            if nargin==1
                roi_in=[];
            end
            if isempty(roi_in) && ~isempty(obj.sumImageGray) %Look for cached copy
                gray_im=obj.sumImageGray;
                roi=[1 obj.sizeX 1 obj.sizeY 1 obj.sizeL obj.globalTBounds];
                return
            end
            [frames,roi]=obj.getFrames(roi_in);
            if isempty(roi) || isempty(frames) 
                error('HSData:getSpectra','Unable to load frames');
            end
            max_hsim=squeeze(max(frames,[],4)); %collapse out time
            mean_hsim=squeeze(mean(frames,4)); %collapse out time
            max_gray_im=max(max_hsim,[],1);
            mean_gray_im=mean(mean_hsim,1);
            
            max_gray_im=max_gray_im./max(max_gray_im(:));
            mean_gray_im=mean_gray_im./max(mean_gray_im(:));
            gray_im=squeeze(max_gray_im+mean_gray_im);
        end
        
        function [rgb_im,roi]=getSumImageRGB(obj,roi_in)
            if nargin==1
                roi_in=[];
            end
            if isempty(roi_in) && ~isempty(obj.sumImageRGB) %Look for cached copy
                rgb_im=obj.sumImageRGB;
                roi=[1 obj.sizeX 1 obj.sizeY 1 obj.sizeL obj.globalTBounds];
                return
            end
            [frames,roi]=obj.getFrames(roi_in);
            cm = obj.colorMap(roi(5):roi(6),:);
            rgb_im = HSData.makeRGBSumImage(frames);
        end
        

        %% Visualization and Data Representations
        
         %Frame viewing
        function viewRawFrames(obj)
            frames=obj.getRawFrames();
            if isempty(frames)
                return
            end
            frames=permute(frames,[2,3,1,4]);
            dipframes=dip_image(frames);
            name=['Raw Frames - ' obj.makeROITitleString(roi_in)];
            obj.viewMaximizedDipFig(dipframes,'Name',name);
        end
        
        function viewFrames(obj, roi_in)
            if nargin==1
                roi_in=[];
            end
            frames=obj.getFrames(roi_in);
            frames=permute(frames,[2,3,1,4]);
            dipframes=dip_image(frames);
            if isempty(dipframes)
                return
            end
            name=['Calibrated Frames - ' obj.makeROITitleString(roi_in)];
            obj.viewMaximizedDipFig(dipframes,'Name',name);
        end      

        function viewFramesRGB(obj, roi_in)
            if nargin==1
                roi_in=[];
            end
            [frames,roi]=obj.getFrames(roi_in);
%             frames = permute(frames,[2,3,1,4]);
%             cm=obj.hyperCM(obj.lambda(roi(3):roi(4)));
%             cm=obj.hyperCM(obj.lambda(roi(5):roi(6)));
            cm=obj.colorMap(roi(5):roi(6),:);
            RGB=HSData.makeRGB(frames, cm );
            dipframes=joinchannels('RGB',RGB);
            if isempty(dipframes)
                return
            end
            name=['Calibrated Frames (RGB)- ' obj.makeROITitleString(roi_in)];
            obj.viewMaximizedDipFig(dipframes,'Name',name);
        end
        
        function f=viewSpectra(obj, roi_in)
            %view the results of plotSpectra in a new figure.
            if nargin==1
                roi_in=[];
            end
            f=figure();
            axH=axes();
            obj.plotSpectra(axH, roi_in);
            title('Spectral Distribution');
            whitebg(f);
        end

        function plotSpectra(obj, axH, roi_in)
            if nargin==2
                roi_in=[];
            end
            axes(axH);
            [I,L]=obj.getSpectra(roi_in);
            area(L,I);
            hold on;
            plot(L,I,'LineWidth',2);
            xlim([L(end) L(1)]);
            xlabel('Wavelength (nm)');
            ylabel('Intensity (photons)');
            hold off;
        end

        function f=viewSpectralSeries(obj, roi_in)
            %view the results of plotSpectralSeries in a new figure.
            if nargin==1
                roi_in=[];
            end
            f=figure();
            whitebg(f);
            obj.plotSpectralSeries(axes(), roi_in);
            title('Spectral-Temporal Distribution');
        end
        
        function plotSpectralSeries(obj, axH, roi_in)
            if nargin==2
                roi_in=[];
            end
            axes(axH);
            [I,L,roi]=obj.getSpectralSeries(roi_in);
            I(:)=max(0,I(:));
            ts=obj.frameT*((roi(7):roi(8))-1);
            surf(ts,L,double(I),'EdgeColor','none');
            ylim([L(end) L(1)]);
            xlim([ts(1) ts(end)]);
            zlim([0 max(I(:))]);
            colormap(jet);
            pbaspect([1 .33 .33])
            ylabel('Wavelength (nm)');
            xlabel('Time (s)');
            zlabel('Intensity (photons)');
            rotate3d(axH);
            box('on');
            set(axH,'Box','on','BoxStyle','full');
            set(axH,'XGrid','on','XMinorGrid','on','XMinorTick','on',...
                   'YGrid','on','YMinorGrid','on','YMinorTick','on',...
                   'ZGrid','on','ZMinorGrid','on','ZMinorTick','on');      
            hold('off');
        end
        
        function viewSumImageGray(obj,roi_in)
            if nargin==1
                roi_in=[];
            end
            [frames,roi]=obj.getSumImageGray(roi_in);
            dipframes=dip_image(frames);
            name=['Sum Image Gray - ' obj.makeROITitleString(roi)];
            obj.viewMaximizedDipFig(dipframes,'Name',name);
        end

        
        function viewSumImageRGB(obj,roi_in)
            if nargin==1
                roi_in=[];
            end
            [RGB,roi]=obj.getSumImageRGB(roi_in);
            
            dipframes=joinchannels('RGB',RGB);
            name=['Sum Image RGB - ' obj.makeROITitleString(roi)];
            obj.viewMaximizedDipFig(dipframes,'Name',name);
        end
        
        function f=viewSurfaceSliceIntensity(obj, frame_idx, roi_in)
            if nargin==1
                frame_idx=0;
            end
            if nargin<=2
                roi_in=[];
            end
            [im,roi]=obj.getPlottableFrame(frame_idx,roi_in);
            [xs,ys,Ls]=obj.pixelGridValues(roi);
            f=figure();
            HSData.surfaceSliceView(xs,ys,Ls,im, jet,'intensity');
            name=sprintf('Frame - %i',roi(7)-1+frame_idx);
            set(f,'Name',name);
            cbh=colorbar();
            set(get(cbh,'Label'),'String','Intensity (Photons)');
            obj.label3DFigure(f, cbh);
        end
        
        function f=viewSurfaceSliceWavelength(obj, frame_idx, roi_in)
            if nargin==1
                frame_idx=0;
            end
            if nargin<=2
                roi_in=[];
            end
            [im,roi]=obj.getPlottableFrame(frame_idx,roi_in);
            [xs,ys,Ls]=obj.pixelGridValues(roi);
            f=figure();
            HSData.surfaceSliceView(xs,ys,Ls,im, jet,'wavelength');
            name=sprintf('Frame - %i',roi(7)-1+frame_idx);
            set(f,'Name',name);
            cbh=colorbar();
            set(get(cbh,'Label'),'String','Wavelength $\lambda$ (nm)');
            obj.label3DFigure(f, cbh);
        end
        
        function f=viewVolumeIntensity(obj, frame_idx, roi_in)
            if nargin==1
                frame_idx=0;
            end
            if nargin<=2
                roi_in=[];
            end
            [im,roi]=obj.getPlottableFrame(frame_idx,roi_in);
            [xs,ys,Ls]=obj.pixelGridValues(roi);
            f=figure();
            ax=axes();
            HSData.volumeSliceView(xs,ys,Ls,im, jet,'intensity',ax);
            name=sprintf('Frame - %i',roi(7)-1+frame_idx);
            set(f,'Name',name);
            cbh=colorbar();
            set(get(cbh,'Label'),'String','Intensity (Photons)');
            obj.label3DFigure(f, cbh);
            view_angles=[171, 2];
            view(view_angles);
            set(gcf(),'Position',[10 10 600 450]);
            sliderStep=[1/(obj.nFrames-1) 10/(obj.nFrames-1)];
            function slider_CB(s,~)
                idx=round(s.Value);
                [im,roi]=obj.getPlottableFrame(idx,roi);
                [view_a1, view_a2]=view();
                cla(ax);
                HSData.volumeSliceView(xs,ys,Ls,im, jet,'intensity',ax);
                name=sprintf('Frame - %i',roi(7)-1+idx);
                set(f,'Name',name);
                set(get(cbh,'Label'),'String','Intensity (Photons)');
                obj.label3DFigure(f, cbh);
                view([view_a1 view_a2]);
            end
            if frame_idx>0
                total_frame=obj.getFrames(roi_in);
                mcount=max(total_frame(:));
                caxis([0 mcount*.75]);
                uicontrol('Style','slider','Min',obj.globalTBounds(1),'Max',obj.globalTBounds(2),'Value',frame_idx,...
                          'Position',[10,435,580,20],'SliderStep',sliderStep, 'Callback',@slider_CB);
            end
        end
        
        function f=viewVolumeWavelength(obj, frame_idx, roi_in)
            if nargin==1
                frame_idx=0;
            end
            if nargin<=2
                roi_in=[];
            end

            [im,roi]=obj.getPlottableFrame(frame_idx,roi_in);
            [xs,ys,Ls]=obj.pixelGridValues(roi);
            f=figure();
            ax=axes();
            if frame_idx>0
                fs=obj.getFrames(roi_in);
                max_val = 0.75*max(fs(:));
                im(:)=min(im(:),max_val);
            end
            HSData.volumeSliceView(xs,ys,Ls,im, jet,'wavelength',ax);
            name=sprintf('Frame - %i',roi(7)-1+frame_idx);
            set(f,'Name',name);
            cbh=colorbar();
            set(get(cbh,'Label'),'String','Wavelength $\lambda$ (nm)');
            obj.label3DFigure(f, cbh);
            view_angles=[171, 2];
            view(view_angles);
            set(gcf(),'Position',[10 10 600 450]);
            sliderStep=[1/(obj.nFrames-1) 10/(obj.nFrames-1)];
            function slider_CB(s,~)
                idx=round(s.Value);
                [im,roi]=obj.getPlottableFrame(idx,roi);
                [view_a1, view_a2]=view();
                cla(ax);
                im(:)=max(0,min(im(:),max_val)-0.02*max_val);
                HSData.volumeSliceView(xs,ys,Ls,im, jet,'wavelength',ax);
                name=sprintf('Frame - %i',roi(7)-1+idx);
                set(f,'Name',name);
                set(get(cbh,'Label'),'String','Wavelength $\lambda$ (nm)');
                obj.label3DFigure(f, cbh);
                view([view_a1 view_a2]);
            end
            if frame_idx>0
                total_frame=obj.getFrames(roi_in);
%                 mcount=max(total_frame(:));
%                 caxis([0 mcount*.75]);
                uicontrol('Style','slider','Min',obj.globalTBounds(1),'Max',obj.globalTBounds(2),'Value',frame_idx,...
                          'Position',[10,435,580,20],'SliderStep',sliderStep, 'Callback',@slider_CB);
            end
        end
        
        
        %% Dependent property get methods
        function val=get.lambda(obj)
            %lambda wavelengths for the spectral bins after baseline rows are removed
            if ~isempty(obj.rawLambda)
                val=obj.rawLambda(1:end-obj.NBaselineRows);
            else
                val=[];
            end
        end

        function val=get.nFrames(obj)
            %number of time points
            if ~isempty(obj.size)
                val=obj.globalTBounds(2)-obj.globalTBounds(1)+1;
            else
                val=[];
            end
        end

        function val=get.sizeX(obj)
            %%number of pixels in the X (scaned) dimension
            if ~isempty(obj.size)
                val=obj.size(3);
            else
                val=[];
            end
        end
        
        function val=get.sizeY(obj)
            %number of pixels in the Y (non-scanned) dimension
            if ~isempty(obj.size)
                val=obj.size(2);
            else
                val=[];
            end
        end
        
        function val=get.sizeL(obj)
            %number of pixels in the L dimension
            if ~isempty(obj.size)
                val=obj.size(1);
            else
                val=[];
            end
        end

        function val=get.sensorSize(obj)
            % effective sensor size after baseline correction is [L Y] 
            if ~isempty(obj.size)
                val=obj.size(1:2);
            else
                val=[];
            end
        end
        
        function val=get.frameSize(obj)
            % effective frame size after baseline correction is [L Y X] 
            if ~isempty(obj.size)
                val=obj.size(1:3);
            else
                val=[];
            end
        end
        
        function val=get.physicalSize(obj)
            % effective physcial frame size in microns [X Y] 
            if ~isempty(obj.size)
                val=obj.pixelSize*[obj.sizeX, obj.sizeY];
            else
                val=[];
            end
        end

        function val=get.NSensorPixels(obj)
            % effective # sensor pixels after baseline correction 
            if ~isempty(obj.size)
                val=prod(obj.size(1:2));
            else
                val=[];
            end
        end
        
        function val=get.NFramePixels(obj)
            % effective #pixels/frame after baseline correction 
            if ~isempty(obj.size)
                val=prod(obj.size(1:3));
            else
                val=[];
            end
        end

        function val=get.NExperimentPixels(obj)
            % effective #pixels for entire experiment after baseline correction 
            if ~isempty(obj.size)
                val=prod(obj.size);
            else
                val=[];
            end
        end

        function val=get.rawSize(obj)
            %The size of raw data for entire experiment [L Y X T]
            if ~isempty(obj.size)
                val=obj.size;
                val(1)=val(1)+obj.NBaselineRows;
            else
                val=[];
            end
        end

        function val=get.rawSizeL(obj)
            %Actual number of pixels in the L dimension in raw data
            if ~isempty(obj.size)
                val=obj.size(1)+obj.NBaselineRows;
            else
                val=[];
            end
        end

        function val=get.rawSensorSize(obj)
            %Actual sensor size [L Y]
            if ~isempty(obj.size)
                val=obj.size(1:2);
                val(1)=val(1)+obj.NBaselineRows;
            else
                val=[];
            end
        end
        
        function val=get.frameT(obj)
            %time(s) to aquire 1 [L Y X] frame
            if ~isempty(obj.size)
                val=obj.sensorT*obj.sizeX;
            else
                val=[];
            end
        end
        
        function val=get.experimentT(obj)
            %time(s) to aquire entire experiment [L Y X T] image
            if ~isempty(obj.size)
                val=obj.sensorT*obj.sizeX*obj.nFrames;
            else
                val=[];
            end
        end
        
        function val=get.icsBasePath(obj)
            %The fullpath without the final number and .ics for the ICS data files
            if ~isempty(obj.workingDir)
                val=fullfile(obj.workingDir,obj.Paths.icsBaseName);
            else
                val=[];
            end
        end
        

    end %Public Methods

    methods (Access=protected)
        function roi=checkROI(obj, roi_in, roi_name)
            if isempty(obj.size)
                error('HSData:checkROI', 'size not set.  Cannot add ROI.');
            end
            if nargin>2 && any(cellfun(@(s) strcmp(roi_name,s),obj.ROIname))
                error('HSData:checkROI', 'Already have ROI name "%s"',roi_name);
            end
            roi=roi_in(:)'; %make a row vector
            if length(roi)==4
                roi=[roi 1 obj.sizeL obj.globalTBounds];
            end
            if length(roi)==6
                roi=[roi obj.globalTBounds];
            end
            roi(1)=max(1,roi(1));
            roi(2)=min(obj.sizeX,roi(2));
            roi(3)=max(1,roi(3));
            roi(4)=min(obj.sizeY,roi(4));
            roi(5)=max(1,roi(5));
            roi(6)=min(obj.sizeL,roi(6));
            roi(7)=max(obj.globalTBounds(1),roi(7));
            roi(8)=min(obj.globalTBounds(2),roi(8));
            if length(roi)~=8 || any(round(roi)~=roi)
                error('HSData:checkROI','Bad ROI format');                 
            elseif roi(2)<roi(1)
                error('HSData:checkROI','Bad ROI X order');
            elseif roi(4)<roi(3)
                error('HSData:checkROI','Bad ROI Y order');
            elseif roi(6)<roi(5)
                error('HSData:checkROI','Bad ROI L order');
            elseif roi(8)<roi(7)
                error('HSData:checkROI','Bad ROI T order');
            end
        end

        function mask=checkLambdaMask(obj, mask, name)
            mask=mask(:)';
            
            if ~isvector(mask) 
                error('HSData:addLambdaMask','Mask is not a vector');
            end
            if length(mask)~=obj.sizeL 
                error('HSData:addLambdaMask','Mask should be length:%i',obj.sizeL);
            end
            if any(mask(:)~=logical(mask(:)))
                error('HSData:addLambdaMask','Mask is not 1/0');
            end
            if any(cellfun(@(m) all(m==mask),obj.lambdaMask))
                warning('HSData:addLambdaMask','Already have equivalent lambda mask pattern.');
            end
            if ~ischar(name)
                error('HSData:addLambdaMask','Name is invalid');
            end
            if any(cellfun(@(n) strcmp(n,name),obj.lambdaMaskName))
                error('HSData:addLambdaMask','Name %s already exists',name);
            end
            mask=logical(mask);
        end
           
        function loadHSData(obj, hsdatafile)
            obj.resetObject();
            [workingDir,~,~]=fileparts(hsdatafile);
            s=load(hsdatafile,'-mat');
            if obj.initialized
                obj.savePreservedProperties();
                pp=obj.preservedProperties;
                obj.resetObject();
                obj.preservedProperties=pp;
            end
            obj.copyobj(s.obj);
            obj.workingDir=workingDir;
            obj.initialized=true;
            obj.dirty=false;
            obj.openIface();
        end

        function loadHSI(obj, hsi, hsifilepath)
            %Common initialization for HSI from object or from .hsi file
            obj.resetObject();
            [path, file, ~]=fileparts(hsifilepath);
            obj.size=[hsi.sz hsi.nt];
            obj.size(1)=obj.size(1)-obj.NBaselineRows; % Account for baseline corrections
            obj.rawLambda=hsi.wv;
            obj.pixelSize=hsi.pixelSize/1000;
            obj.sensorT=hsi.t;
            obj.globalTBounds=[1 hsi.nt];
            obj.workingDir=path;
            acqInfo=strsplit(hsi.acqInfo);
            acqdatetime=datevec(strjoin(acqInfo(2:3),' '),'yyyy/mm/dd HH:MM:SS');
            obj.AcquisitionParams.acqDateTime=acqdatetime;
            obj.AcquisitionParams.wvCalDateTime=hsi.wvAcqTime;
            obj.AcquisitionParams.laserPower=hsi.LaserPower;
            obj.AcquisitionParams.sensorTemperature=hsi.temperature;
            obj.AcquisitionParams.cameraData=hsi.iXon;
            obj.psf=obj.default_psf;
            obj.Paths.saveFile=[file HSData.saveFileExt];
            obj.Paths.icsBaseName=[file '_'];
            obj.findBackgroundFile();
            obj.findCalibrationFile();
            obj.setLambdaMask();
            obj.setColorMap();           
            obj.initialized=true;
            obj.openIface();
            %Post initialization
            obj.sumImageGray=obj.getSumImageGray();
            obj.sumImageRGB=obj.getSumImageRGB();
            obj.spectra=obj.getSpectra();
            obj.spectralSeries=obj.getSpectralSeries();            
        end

        function findBackgroundFile(obj)
            bgFiles=dir(fullfile(obj.workingDir,'HSM_Background-*.mat'));
            if isempty(bgFiles)
                error('HSData:initializeBackground', 'Unable to locate background file.');
            elseif length(bgFiles)>1
                warning('HSData:initializeBackground', 'Located more than one background files.');
            end
            bgFile=fullfile(obj.workingDir, bgFiles(end).name); %Use last listed file
            obj.loadBackground(bgFile);
        end
        

        function findCalibrationFile(obj)
            calFiles=dir(fullfile(obj.workingDir,'HSM_Calibration-*.mat'));
            if isempty(calFiles)
                error('HSData:initializeCalibration', 'Unable to locate Calibration file.');
            elseif length(calFiles)>1
                warning('HSData:initializeCalibration', 'Located more than one Calibration files.');
            end
            calFile=fullfile(obj.workingDir, calFiles(end).name);
            obj.loadCalibration(calFile);
        end

        

        function success=openIface(obj)
            %Uset the obj.ifaceHandle to open the connection with the C++ class and return
            %the objectHandle to the allocated HSData instance.
            if ~obj.initialized
                error('HSData:openIface','Object not initialized!');
            end
            if ~obj.objectHandle
                %Check raw_size input
                int_raw_size=int32(obj.rawSize);
                assert(all(int_raw_size>0));
                assert(length(int_raw_size)==4);
                
                %Check gain input
                if ~isscalar(obj.CCDGain) || obj.CCDGain<=0
                    error('HSData:openIface','Invalid CCDGain: %f', obj.CCDGain);
                end
                
                %Check bgMean input
                single_bgMean=single(obj.CCDBackground);
                if ~all(size(single_bgMean)==obj.rawSensorSize) || ~all(single_bgMean(:)>0)
                    error('HSData:openIface','Invalid CCDBackground');
                end
                
                %Check basefilename input
                ics0=sprintf('%s%05i.ics',obj.icsBasePath,0);
                icsnT=sprintf('%s%05i.ics',obj.icsBasePath,obj.nFrames-1);
                assert(exist(ics0,'file')==2); %Check for first ics file
                assert(exist(icsnT,'file')==2); %Check for last ics file
                
                %Make the HSData object
                % (in) raw_size [1x4] int32 [LYXT]
                % (in) gain float scalar >0
                % (in) bgMean [raw_sizeL X sizeY] single >0
                % (in) basefilename string, the complete filename with path but without number and .ics
                obj.makeCObj(int_raw_size, single(obj.CCDGain), single_bgMean, obj.icsBasePath);
            end
            success=obj.objectHandle>0;
        end
    
        function load_raw_frames(obj)
            if ~obj.raw_frames_loaded
                obj.updateWaitbar(0,'Loading raw data');
                obj.raw_frames_=obj.call('loadRawFrames');
                obj.raw_frames_loaded= ~isempty(obj.raw_frames_) && all(size(obj.raw_frames_)==obj.rawSize);
                obj.updateWaitbar(1);
            end
        end

        function load_frames(obj)
            if ~obj.frames_loaded
                obj.updateWaitbar(0,'Loading gain corrected data');
                obj.frames_=obj.call('loadFrames');
                obj.frames_loaded= ~isempty(obj.frames_) && all(size(obj.frames_)==obj.size);
                obj.updateWaitbar(1);
            end            
        end

        function name=makeROITitleString(obj,roi_in)
            %Make a nice title for the window
            roi=obj.getROI(roi_in);
            if ~isempty(roi)
                roi_cell=arrayfun(@(a)a,roi_in,'Uniform',0);
                if isscalar(roi_in)
                    name=sprintf('<%s>',obj.ROIname{roi_in});
                else
                    name='ROI';
                end
                name=[name sprintf('{x:[%i,%i],y:[%i,%i],L:[%i,%i],t:[%i,%i]}',roi_cell{:})];
            else
                name=sprintf('Full Frame [X,Y,L,T]:[%i,%i,%i,%i]',obj.sizeX,obj.sizeY,obj.sizeL,obj.nFrames);
            end
        end

        function [xs,ys,Ls,ts]=pixelGridValues(obj,roi_in)
            %
            % Given an roi, return the physical grid for the frame in each dimension
            % In:
            %  roi - fully specified ROI using index values, 
            %        size:[1,8],[xmin, xmax, ymin, ymax, Lmin, Lmax, tmin, tmax]
            % Out:
            %  xs - [1,roi(2)-roi(1)+2] pixel boundary X values (micron)
            %  ys - [1,roi(4)-roi(3)+2] pixel boundary X values (micron)
            %  Ls - [1,roi(6)-roi(5)+2] pixel boundary X values (micron)
            %  ts - [1,roi(8)-roi(7)+2] pixel boundary X values (micron)
            if nargin<1
                roi_in=[];
            end
            roi=obj.getROI(roi_in);
            xs=((roi(1):roi(2)+1)-1)*obj.pixelSize;
            ys=((roi(3):roi(4)+1)-1)*obj.pixelSize;
            lambdas=obj.lambda(roi(5):roi(6));
            dl=diff(lambdas);
            Ls=[lambdas(1)-0.5*dl(1), lambdas(1:end-1)+0.5*dl, lambdas(end)+0.5*dl(end) ];
            ts=((roi(7):(roi(8)+1))-1)*obj.frameT;
        end
        
        function [frame,roi]=getPlottableFrame(obj, frame_idx, roi_in)
            roi=getROI(obj, roi_in);
            if frame_idx>0
                roi(7:8)=frame_idx;
                [frame,roi]=obj.getFrames(roi);
            else
                [frame,roi]=obj.getFrames(roi);
                frame=sum(frame,4);
            end
            frame=permute(frame,[2,3,1]);
        end
        
        function label3DFigure(obj, f, cbh)
            % f - figure handle
            % cbh - colorbar handle
            fig_bg_color=[0 0 0];
            ax_txt_color=[1 1 1];
            ax_font_size=12;
            xlabel('x ($\mu$m)','interpreter','latex','Color',ax_txt_color,'FontSize',ax_font_size);
            ylabel('y ($\mu$m)','interpreter','latex','Color',ax_txt_color,'FontSize',ax_font_size);
            zlabel('$\lambda$ (nm)','interpreter','latex','Color',ax_txt_color,'FontSize',ax_font_size);
            set(f,'Color',fig_bg_color);
            set(cbh,'Color',ax_txt_color);
            set(get(cbh,'Label'),'interpreter','latex','FontSize',ax_font_size);
        end
        
        %IFace Calling
        function varargout=call(obj, cmdstr, varargin)
            if ~obj.objectHandle && ~obj.open_iface()
                error('HSData:call','HSData_Iface could not be created.');
            end
            [varargout{1:nargout}]=obj.ifaceHandle(cmdstr,obj.objectHandle, varargin{:});
        end

        function makeCObj(obj, varargin)
            % Make a new C++ object and save the numeric handle to the allocated object
            obj.objectHandle = obj.ifaceHandle('@new', varargin{:});
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

    methods (Static)     
        
        RGB=hyperCM(lambda, mask) %in seperate file
        
        function RGB=makeRGB(image, cm)
            %
            % Tested at 20x speed of HSIData.makeRGB(), and with significantly cleaner
            % background.  Any RGB value less than 0 is trucated to 0.
            %
            % [in] cm - a Nx3 array of colors mapping the Nth wavelength to the Nth color.
            % [out] im - [y,x] or [y,x,t] RGB dip_image or dip_image sequence.
            %
            imsize=size(image);
            sizeL=imsize(1);
            if nargin<2
                cm=HSData.hyperCM(sizeL:-1:1);
            end
            newsize=imsize(2:end);
            npixels=prod(newsize);
            pixels=reshape(image, sizeL, npixels); % treat as 2D matrix of sizeL  col:Pixel
            RGB= (cm' * pixels)';
            RGB(:)=max(0,RGB(:));
            minval=min(RGB,[],1);
            maxval=max(RGB,[],1);

            %Normalize
            if maxval(1)==minval(1)
                RGB(:,1)=0;
            else
                RGB(:,1)=(RGB(:,1)-minval(1))/(maxval(1)-minval(1));
            end
            if maxval(2)==minval(2)
                RGB(:,2)=0;
            else
                RGB(:,2)=(RGB(:,2)-minval(2))/(maxval(2)-minval(2));
            end
            if maxval(3)==minval(3)
                RGB(:,3)=0;
            else
                RGB(:,3)=(RGB(:,3)-minval(3))/(maxval(3)-minval(3));
            end
            RGB=reshape(RGB, [newsize, 3]);
        end

        function rgb_im = makeRGBSumImage(frames, varargin)
            max_hsim=squeeze(max(frames,[],4)); %collapse out time
            mean_hsim=squeeze(mean(frames,4)); %collapse out time
            rgb_im = HSData.makeRGB(max_hsim+mean_hsim, varargin{:});
        end

        function viewRGB(varargin)
            dipframes=joinchannels('RGB', HSData.makeRGB(varargin{:}));
            if isempty(dipframes)
                return
            end
            name='RGB HSData';
            GUIBuilder.viewDip(dipframes,'Name',name);
        end

        function ax=surfaceSliceView(xs,ys,lambdas,im,cm,mode)
            if nargin==5
                mode='intensity';
            end
            max_alpha=0.7;
            lambda_aspect_ratio=2;
            sz=size(im);
            mxv=max(im(:));

            [X,Y,Z]=meshgrid(double(xs),double(ys),double(lambdas));
            A=zeros(sz+1);
            switch mode
                case 'intensity'
                    C=zeros(sz+1);
                    C(1:end-1,1:end-1,1:end-1)=im;                    
                    A=max(0,C)./mxv;%Normalized 0-1
                case 'wavelength'
                    C=repmat(reshape(lambdas,1,1,sz(3)+1),sz(1)+1,sz(2)+1,1);
                    A(1:end-1,1:end-1,1:end-1)=im;                    
                    A=max(0,A)./mxv;%Normalized 0-1                    
            end
            %Alpha data
            A(A<0.05)=0;
            A=double(A);
            ax=axes('SortMethod','childorder','XDir','reverse');
            hold('on');
            for L=1:sz(3)
                surface('XData',X(:,:,L),...
                        'YData',Y(:,:,L),...
                        'ZData',Z(:,:,L),...
                        'CData',C(:,:,L),...
                        'AlphaData',64*max_alpha*A(:,:,L),...
                        'FaceColor','flat','FaceAlpha','flat','EdgeColor','none',...
                        'AlphaDataMapping','direct');
            end
            colormap(cm);
            axis([xs(1) xs(end) ys(1) ys(end) min(lambdas) max(lambdas)]);
            x_stretch=(xs(end)-xs(1))/(length(xs)-1);
            y_stretch=(ys(end)-ys(1))/(length(ys)-1);
            lambda_stretch=abs(lambdas(1)-lambdas(end))/(length(lambdas)-1);
            daspect([x_stretch y_stretch lambda_stretch*lambda_aspect_ratio]);
            view(135,54);
            box('on');
            set(ax,'Box','on','BoxStyle','full');
            set(ax,'XGrid','on','XMinorGrid','on','XMinorTick','on',...
                   'YGrid','on','YMinorGrid','on','YMinorTick','on',...
                   'ZGrid','on','ZMinorGrid','on','ZMinorTick','on'); 
            ax_bg_color=[0 0 0];
            ax_txt_color=[1 1 1];            
            set(ax,'Color',ax_bg_color,'XColor',ax_txt_color,'YColor',ax_txt_color,'ZColor',ax_txt_color);
        end
        
        function ax=volumeSliceView(xs,ys,lambdas,im,cm,mode,ax)
            % xs - size=sizeX+1 -  The value of X on each pixel edge including the first and last.
            % xs - size=sizeY+1 -  The value of Y on each pixel edge including the first and last.
            % lambdas - size=sizeL+1 -  The value of Lambda on each pixel edge including the first and last.
            %
            % im - size(im)=[sizeY,sizeX,sizeL] - this is the format the
            %       matlab expects for surface plots
            %
            if nargin<7
                ax=axes();
            end
            if nargin<6
                mode='intensity';
            end
            max_alpha=0.7;
            lambda_aspect_ratio=3;
            sz=size(im);
            mxv=max(im(:));
            xs=xs-xs(1); %Remove these to allow different range
            ys=ys-ys(1);
            [X,Y,Z]=meshgrid(double(xs),double(ys),double(lambdas));
            A=zeros(sz+1);
            switch mode
                case 'intensity'
                    C=zeros(sz+1);
                    C(1:end-1,1:end-1,1:end-1)=im;                    
                    A=max(0,C)./mxv;%Normalized 0-1
                case 'wavelength'
                    C=repmat(reshape(lambdas,1,1,sz(3)+1),sz(1)+1,sz(2)+1,1);
                    A(1:end-1,1:end-1,1:end-1)=im;                    
                    A=max(0,A)./mxv;%Normalized 0-1                    
            end
            %Alpha data
%             A(A<0.1)=0;
            A=double(A);
            set(ax,'SortMethod','childorder','Clipping','off','Projection','perspective','XDir','reverse');
            hold('on');
            CX=C;
            CY=C;
            CZ=C;
            CX(:,2:end,:)=max(CX(:,1:end-1,:),CX(:,2:end,:));
            CY(2:end,:,:)=max(CY(1:end-1,:,:),CY(2:end,:,:));
            CZ(:,:,2:end)=max(CZ(:,:,1:end-1),CZ(:,:,2:end));
            AX=A;
            AY=A;
            AZ=A;
            AX(:,2:end,:)=max(AX(:,1:end-1,:),AX(:,2:end,:));
            AY(2:end,:,:)=max(AY(1:end-1,:,:),AY(2:end,:,:));
            AZ(:,:,2:end)=max(AZ(:,:,1:end-1),AZ(:,:,2:end));
            for x=1:sz(2)+1
                surface('XData',squeeze(X(:,x,:)),...
                       'YData',squeeze(Y(:,x,:)),...
                       'ZData',squeeze(Z(:,x,:)),...
                       'CData',squeeze(CX(:,x,:)),...
                       'FaceColor','flat','FaceAlpha','flat','EdgeColor','none',...
                       'AlphaDataMapping','direct',...
                       'AlphaData',64*max_alpha*squeeze(AX(:,x,:)));
            end
            for z=1:sz(3)+1
                surface('XData',squeeze(X(:,:,z)),...
                       'YData',squeeze(Y(:,:,z)),...
                       'ZData',squeeze(Z(:,:,z)),...
                       'CData',squeeze(CZ(:,:,z)),...
                       'FaceColor','flat','FaceAlpha','flat','EdgeColor','none',...
                       'AlphaDataMapping','direct',...
                       'AlphaData',64*max_alpha*squeeze(AZ(:,:,z)));
            end
            for y=1:sz(1)+1
                surface('XData',squeeze(X(y,:,:)),...
                       'YData',squeeze(Y(y,:,:)),...
                       'ZData',squeeze(Z(y,:,:)),...
                       'CData',squeeze(CY(y,:,:)),...
                       'FaceColor','flat','FaceAlpha','flat','EdgeColor','none',...
                       'AlphaDataMapping','direct',...
                       'AlphaData',64*max_alpha*squeeze(AY(y,:,:)));
            end
            
            
            
            colormap(cm);
            axis([xs(1) xs(end) ys(1) ys(end) min(lambdas) max(lambdas)]);
            x_stretch=(xs(end)-xs(1))/(length(xs)-1);
            y_stretch=(ys(end)-ys(1))/(length(ys)-1);
            lambda_stretch=abs(lambdas(1)-lambdas(end))/(length(lambdas)-1);
            daspect([x_stretch y_stretch lambda_stretch*lambda_aspect_ratio]);
            view(135,54);
            box('on');
            set(ax,'Box','on','BoxStyle','full');
            set(ax,'XGrid','on','XMinorGrid','on','XMinorTick','on',...
                   'YGrid','on','YMinorGrid','on','YMinorTick','on',...
                   'ZGrid','on','ZMinorGrid','on','ZMinorTick','on');
            ax_bg_color=[0 0 0];
            ax_txt_color=[1 1 1];            
            set(ax,'Color',ax_bg_color,'XColor',ax_txt_color,'YColor',ax_txt_color,'ZColor',ax_txt_color);
            hold off
        end

        
        function ax=pixelVolumeView(xs,ys,lambdas,im,cm,mode)
            if nargin==5
                mode='intensity';
            end
            max_alpha=0.5;
            lambda_aspect_ratio=2;
            sz=size(im);
            mxv=max(im(:));

            A=zeros(sz+1);
            switch mode
                case 'intensity'
                    C=zeros(sz+1);
                    C(1:end-1,1:end-1,1:end-1)=im;                    
                    A=max(0,C)./mxv;%Normalized 0-1
                case 'wavelength'
                    C=repmat(reshape(lambdas,1,1,sz(3)+1),sz(1)+1,sz(2)+1,1);
                    A(1:end-1,1:end-1,1:end-1)=im;                    
                    A=max(0,A)./mxv;%Normalized 0-1                    
            end
            %Alpha data
            ax=axes('SortMethod','childorder');
            hold('on');
            min_alpha=0.05;
            Pi=1;
            MaxP=prod(100000);
            Px=zeros(4,MaxP);
            Py=zeros(4,MaxP);
            Pz=zeros(4,MaxP);
            Pc=zeros(1,MaxP);
            Pa=zeros(1,MaxP);
            plotprops={'EdgeColor','none'};
            for Xi=1:sz(2)
                for Yi=1:sz(1)
                    for Zi=1:sz(3)
                        if A(Yi,Xi,Zi)<min_alpha
                            continue
                        end
                        if Pi>MaxP
                            warning('HSData:pixelVolumeView','Too many faces');
                            break
                        end
                        x0=xs(Xi);
                        y0=ys(Yi);
                        z0=lambdas(Zi);
                        x1=xs(Xi+1);
                        y1=ys(Yi+1);
                        z1=lambdas(Zi+1);
                        c=C(Yi,Xi,Zi);
                        a=A(Yi,Xi,Zi);
                        if Yi==1 || A(Yi-1,Xi,Zi)<a %front face
                            Px(:,Pi)=[x0;x1;x1;x0];
                            Py(:,Pi)=y0;
                            Pz(:,Pi)=[z0;z0;z1;z1];
                            Pc(Pi)=c;
                            Pa(Pi)=a;
                            Pi=Pi+1;
                        end
                        if Yi==sz(1) || A(Yi+1,Xi,Zi)<a %rear face
                            Px(:,Pi)=[x0;x1;x1;x0];
                            Py(:,Pi)=y1;
                            Pz(:,Pi)=[z0;z0;z1;z1];
                            Pc(Pi)=c;
                            Pa(Pi)=a;
                            Pi=Pi+1;
                        end
                        if Xi==1 || A(Yi,Xi-1,Zi)<a %front face
                            Px(:,Pi)=x0;
                            Py(:,Pi)=[y0;y1;y1;y0];
                            Pz(:,Pi)=[z0;z0;z1;z1];
                            Pc(Pi)=c;
                            Pa(Pi)=a;
                            Pi=Pi+1;
                        end
                        if Xi==sz(2) || A(Yi,Xi+1,Zi)<a %rear face
                            Px(:,Pi)=x1;
                            Py(:,Pi)=[y0;y1;y1;y0];
                            Pz(:,Pi)=[z0;z0;z1;z1];
                            Pc(Pi)=c;
                            Pa(Pi)=a;
                            Pi=Pi+1;
                        end
                        
                        if Zi==1 || A(Yi,Xi,Zi-1)<a %front face
                            Px(:,Pi)=[x0;x1;x1;x0];
                            Py(:,Pi)=[y0;y0;y1;y1];
                            Pz(:,Pi)=z0;
                            Pc(Pi)=c;
                            Pa(Pi)=a;
                            Pi=Pi+1;
                        end
                        if Zi==sz(3) || A(Yi,Xi,Zi+1)<a %rear face
                            Px(:,Pi)=[x0;x1;x1;x0];
                            Py(:,Pi)=[y0;y0;y1;y1];
                            Pz(:,Pi)=z1;
                            Pc(Pi)=c;
                            Pa(Pi)=a;
                            Pi=Pi+1;
                        end
                    end
                end
            end
            Px(:,Pi:end)=[];
            Py(:,Pi:end)=[];
            Pz(:,Pi:end)=[];
            Pc(:,Pi:end)=[];
            Pa(:,Pi:end)=[];
            hs=fill3(Px,Py,Pz,Pc,plotprops{:});
            for hi=1:length(hs)
                set(hs(hi),'FaceAlpha','flat','AlphaDataMapping','direct','FaceVertexAlphaData',64*max_alpha*Pa(hi));
            end
            colormap(cm);
            axis([xs(1) xs(end) ys(1) ys(end) min(lambdas) max(lambdas)]);
            x_stretch=(xs(end)-xs(1))/(length(xs)-1);
            y_stretch=(ys(end)-ys(1))/(length(ys)-1);
            lambda_stretch=abs(lambdas(1)-lambdas(end))/(length(lambdas)-1);
            daspect([x_stretch y_stretch lambda_stretch*lambda_aspect_ratio]);
            view(135,54);
            box('on');
            set(ax,'Box','on','BoxStyle','full');
            set(ax,'XGrid','on','XMinorGrid','on','XMinorTick','on',...
                   'YGrid','on','YMinorGrid','on','YMinorTick','on',...
                   'ZGrid','on','ZMinorGrid','on','ZMinorTick','on');      
            hold('off');
            ax_bg_color=[0 0 0];
            ax_txt_color=[1 1 1];            
            set(ax,'Color',ax_bg_color,'XColor',ax_txt_color,'YColor',ax_txt_color,'ZColor',ax_txt_color);
        end
        
    end % Public static methods
   
end

