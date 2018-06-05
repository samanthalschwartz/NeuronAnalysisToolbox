classdef ChannelRegistrationV1 < handle
    %ChannelRegistrationV1 Scanning Fiducial Bead Correction
    %   Aquire and use the postions found by scanning a bead that is
    %   observed in multiple (color) channels on the same camera.
    %   Example of four channels in a 2x2 layout and are indexed
    %   [1 3;2 4]
    %
    % ChannelRegistrationV1 Properties:
    %                    timeStamp - acquistion date and time as date number.
    %                     CameraID - camera identification. For future use?
    %                      StageID - stage Identification. Currently only
    %                               implemented for MCLNanoDrive
    %                    PixelSize - pixel size in microns. Default 16/150.
    %                      BoxSize - size of fitting box in pixels. Default 20.                      
    %                SearchBoxSize - size of search box for bead in pixels. Default 30. 
    %                     PSFsigma - point spread function sigma for each
    %                               channel in pixels. Default [1.3 1.3; 1.3 1.3]
    %                        xdist - spacing between fiducial sampling in x
    %                               dimension in microns. Default determined 
    %                               by Sampling for SamplingXY = [5 5]
    %                        ydist - spacing between fiducial sampling in y
    %                               dimension in microns. Default determined 
    %                               by Sampling for SamplingXY = [5 5]
    %     ChannelRegistrationArray - binary array indicating channel setup on
    %                               CCD and which channels to perform
    %                               registration on. Example [1 0;1 1] is 
    %                               four channel setup performing
    %                               registration on channels 1 2 4 where
    %                               channel indices are [1 3;2 4]
    %                MasterChannel - channel in ChannelRegistrationArray to
    %                               transform all other channels to
    %                    MCLhandle - handle for MCLNanoDrive GUI
    %                andor_handles - structured array with handles to
    %                               AndoriXon GUI
    %                   trnLogical - logical for dividing acquisition into
    %                               training and test sets. dimensions are 
    %                               equal to SamplingXY
    %                    ImageSize - size of image acquired on CCD in pixels
    %                   ImageStack - 4-D dip_image object with dimensions 
    %                               BoxSize x BoxSize x ImageNumber x Channel. 
    %                               Contains acquired roi
    %              CoordinateArray - Absolute localization coordinates in
    %                               pixels with dimensions ImageNumber x Channel x Dim
    %                BoxStartPixel - Absolute start position for each roi in
    %                               pixels ImageNumber x Channel x Dim
    %                      CCDgain - gain for CCD. Default 28.
    %                    CCDoffset - offset for CCD. Default 120.
    %                 BeadSelectXY - Absolute position in pixels of user click for
    %                               bead selection
    %                       BeadXY - Absolute position of localized bead identified 
    %                               from user click in microns.
    %                  MCLStartXYZ - Start position of stage at the begining
    %                               of acquisition
    %                   SamplingXY - number of calibration points in x and
    %                               y. Total number of calibration points
    %                               is prod(SamplingXY)
    %                        tform - structured array transformation
    %                               information for each channel into the
    %                               MasterChannel
    %                       tstObj - ChannelRegistration object for testing
    %                               the calibration
    %                            N - number of nearest neighbors to use for
    %                               locally weighted transform.
    %                               Default 5.
    %                            C - order of polynomial fit for transformation.
    %                               Polynomial function is 
    %                               1 + x + y + xy + x^2 + y^2
    %
    % ChannelRegistrationV1 Methods:
    %       ChannelRegistrationV1 - Constructor Finds Andor Camera and MCL GUI
    %       copyobj - copy object
    %       GetStageID - check for stage
    %       SetSamplingGrid - set xdist and ydist in microns
    %       StartBeadCapture - start fiducial acquisition
    %       UserSelectBead - Get ROI and Bead Position
    %       GetFitROI - acquire and fit bead in ROI
    %       GetFitROImaxAdj - ??
    %       CalcPositions - calculates CoordinateArray using ImageStack and BoxStartPixel
    %       CalcPixelSize - use the acquired data to determine pixelsize
    %       plotSubRegions - plot subregions for ChannelRegistration
    %       plot - plot localizations for ChannelRegistration
    %       plotTform - plot transformation results for ChannelRegistration
    %       calcTform - calculate transform into master channel
    %       tstTform - test transform into master channel
    %       locTform - transform localizations from ch into MasterChannel
    %       NoCalcTform - check that calcTform has been performed
    %       LWtform - Locally weighted transformation of coordinates using tform
    %       CalcLWtform - Calculate Locally Weighted tform
    %       GetStagePosition - ??
    %       SetStagePosition - ??
    %       
    %
    %Dependencies:
    %   cMakesubreions
    %   GPUgaussMLEv2
    %   64 bit matlab with Visual Studio Pro installed.
    %
    % $Rev: 502 $
    % $Date: 2013-06-21 11:56:44 -0600 (Fri, 21 Jun 2013) $
    % $Author: sheng $
    % $HeadURL: https://abbe.phys.unm.edu/svn/MATLAB/Instrumentation/ChannelRegistrationV1.m $
    %edited by Sam/Pat/Jason 110621 to have coefficient properties
    %
    % NOTE: Currently dependent on MCLNanoDrive GUI for acquisition
    %
    % See also CollectChannelRegistration
    
    properties
        % 0,0 is upper left for everything!
        % x = image cols, y = image rows
        timeStamp = []; %acquistion date and time as date number.
        CameraID=0; %camera identification. For future use?
        StageID=0; %stage Identification. Currently only implemented for MCLNanoDrive
        PixelSize=16/150;    %microns -needs to be defined for each instrument. For TIRF =16/150.
        BoxSize=20;     %size of fitting box in pixels. Default 20.
        SearchBoxSize=30;     %size of search box for bead in pixels. Default 30.
        PSFsigma=[1.3 1.3; 1.3 1.3];     %point spread function sigma for each channel in pixels. Default [1.3 1.3; 1.3 1.3]
        %xdist - spacing between fiducial sampling in x dimension in 
        %microns. Default determined by Sampling for SamplingXY = [5 5]
        xdist=0;
        %ydist - spacing between fiducial sampling in y dimension in 
        %microns. Default determined by Sampling for SamplingXY = [5 5]
        ydist=0;
        %     ChannelRegistrationArray - binary array indicating channel setup on
        %                               CCD and which channels to perform
        %                               registration on. Example [1 0;1 1] is
        %                               four channel setup performing
        %                               registration on channels 1 2 4 where
        %                               channel indices are [1 3;2 4]
        ChannelRegistrationArray=[0 0; 0 0]; 
        MasterChannel; % channel in ChannelRegistrationArray to transform all other channels to
        MCLhandle=0; %handle for MCLNanoDrive GUI
        andor_handles=0; %structured array with handles to AndoriXon GUI
        trnLogical = 0; %logical array with size SamplingXY indicating which point to use to train and testing the calibration
        ImageSize=0;    %size of image acquired on CCD in pixels
        ImageStack=0;   %4-D dip_image object with dimensions BoxSize x BoxSize x ImageNumber x Channel. Contains acquired roi
        CoordinateArray;  %Absolute localization coordinates in pixels with dimensions ImageNumber x Channel x Dim
        BoxStartPixel;  %Absolute start position for each roi in pixels ImageNumber x Channel x Dim
        CCDgain=28; %gain for CCD. Default 28.
        CCDoffset=120; %offset for CCD. Default 120.
        BeadSelectXY=0; %Absolute position in pixels of user click for bead selection
        BeadXY=0;  %Absolute position of localized bead identified from user click in microns.
        MCLStartXYZ=0; %Start position of stage at the begining of acquisition
        SamplingXY = []; %number of calibration points in x and y. Total number of calibration points is prod(SamplingXY)
        %tform - structured array transformation information for each 
        %channel into theMasterChannel
        tform = struct('Measured',[],'Sensed',[],'pX',[],'pY',[],'svec',[],'trn',[],'tst',[]);
        tstObj = 0; %ChannelRegistration object for testing the calibration
        N = 5; %number of nearest neighbors to use for locally weighted transform. Default 5.
        C = 4; %order of polynomial fit for transformation. Polynomial function is 1 + x + y + xy + x^2 + y^2
    end
    
    methods
        %Constructor: Finds Andor Camera and MCL GUI
        function obj=ChannelRegistrationV1(varargin)
            % Avoid being in dipimage menu
            if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'DIP_GetParamList')
                %obj = struct('menu','none');
                error('Go away Dipimage');
            end
            if nargin == 0 %return default settings
                 error('ChannelRegistrationV1:WrongNumberOfInputs',...
                    'Usage: ChannelRegistrationV1(PixelSize,ChannelRegistrationArray) or ChannelRegistrationV1(DefaultCRobj)');
                return
            end
            %error checking
            if nargin == 1 && ~isa(varargin{1},'ChannelRegistrationV1')
                error('ChannelRegistrationV1:WrongNumberOfInputs',...
                    'Usage: ChannelRegistrationV1(PixelSize,ChannelRegistrationArray) or ChannelRegistrationV1(DefaultCRobj)');
            end
            if nargin > 2
                error('ChannelRegistrationV1:WrongNumberOfInputs',...
                    'Usage: ChannelRegistrationV1(PixelSize,ChannelRegistrationArray) or ChannelRegistrationV1(DefaultCRobj)');
            end
            %set default values
            if nargin == 1
                copyobj(varargin{1},obj); 
            else
                obj.PixelSize=varargin{1};
                obj.ChannelRegistrationArray=varargin{2};
            end
            ChannelRegistrationV1.versn;
            andor=findall(0,'Name','AndoriXon');
            if isempty(andor)
                error('ChannelRegistrationV1:NoCamera','AndoriXon.m must be running')
            end
            set(andor,'handlevisibility','callback');
            obj.andor_handles = guidata(andor);
            %force andor to calculate image sizes and setup for 'capture'
            AndoriXon('ActionCaptureButton_Callback',obj.andor_handles.ActionCaptureButton,[],obj.andor_handles);
            close(findall(0,'name','capture'))
            obj.andor_handles = guidata(andor);
            obj.andor_handles.t = [];
            obj.ImageSize=obj.andor_handles.ImageSize;
            obj.timeStamp=now;
            
            %default to 5x5 grid
            %updated PJC 2011-08-29
            if ~obj.xdist && ~obj.ydist
                if isempty(obj.SamplingXY)
                    obj = SetSamplingGrid(obj,[5 5]);
                else
                    obj = SetSamplingGrid(obj,obj.SamplingXY);
                end
            end

        end
        
        function copyobj(thisObj,newObj)
            % Construct a new object based on a deep copy of the current
            % object of this class by copying properties over.
            props = properties(thisObj);
            for i = 1:length(props)
                % Use Dynamic Expressions to copy the required property.
                % For more info on usage of Dynamic Expressions, refer to
                % the section "Creating Field Names Dynamically" in:
                % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                if any(strcmp({'xdist' 'ydist'}, props{i})) && ~thisObj.(props{i})
                    continue; %enable copying of xdist, ydist properties with value of zero
                end 
                newObj.(props{i}) = thisObj.(props{i});
            end
        end
        
        function obj = GetCameraID(obj)
            
        end
        
        function obj = GetStageID(obj)
            %code that checks for MCL stage
            tmp=findall(0,'name','MCLNanoDrive');
            
            if ~isempty(tmp); obj.StageID='MCL';
                h=guidata(tmp);
                obj.MCLhandle = h.MCLhandle;
                return;
            end;
            %code that checks for Prior Stage
            if 0
                obj.StageID='Prior';
                return;
            end;
            error('ChannelRegistrationV1:NoStage','ChannelRegistrationV1: No stage found. MCL Stage must have MCLNanoDrive GUI Active.');
        end
        
        function obj = SetSamplingGrid(obj,SamplingXY)
            % SetSamplingGrid(x,y) set xdist and ydist in microns based on an x by y
            % grid.
            %updated PJC 2011-08-29
            [Ny Nx] = size(obj.ChannelRegistrationArray);
            ChannelWindowSize = obj.PixelSize*obj.ImageSize./([Nx Ny]); %micron
            %adds a buffer of BoxSize on each side fo sampling
            tmp = (ChannelWindowSize-3*obj.BoxSize*obj.PixelSize)./SamplingXY;
            obj.SamplingXY = SamplingXY;
            obj.xdist = tmp(1);
            obj.ydist = tmp(2);
            obj.trnLogical = false(SamplingXY);
            obj.trnLogical(1:2:end) = true;
        end
        
        function obj = StartBeadCapture(obj)
            %start fiducial acquisition
            
            close(findall(0,'name','Channel Registration SubRegions'))
            
            %Get StageID
            obj = GetStageID(obj);
            
            %Get Selected Bead Coordinates
            obj=UserSelectBead(obj); %sets obj.BeadXY in micron
            
            %number of registered channels
            Nreg=sum(sum(obj.ChannelRegistrationArray));
            
            %updated PJC 2011-08-29
            [Ny Nx] = size(obj.ChannelRegistrationArray);
            ChannelWindowSize = obj.PixelSize*obj.ImageSize./[Nx Ny]; %micron
            
            
            %Calculate the number of Relative Stage Positions Needed
            %updated PJC 2011-08-29
            %adds a buffer of BoxSize on each side fo sampling
            obj.SamplingXY = round((ChannelWindowSize-3*obj.BoxSize*obj.PixelSize)./[obj.xdist obj.ydist]);
            
            if any(size(obj.trnLogical) ~= obj.SamplingXY)
                obj.trnLogical = false(obj.SamplingXY);
                obj.trnLogical(1:2:end) = true;
            end
            
            %channel offsets from top left (microns)
            %updated PJC 2011-08-29
            XchOffsets = zeros(size(obj.ChannelRegistrationArray));
            YchOffsets = zeros(size(obj.ChannelRegistrationArray));
            if Ny == 1
                for ii = 1:Nx
                    XchOffsets(ii) = (ii-1)*ChannelWindowSize(1);
                end
            else if Nx == 1
                    for ii = 1:Ny
                        YchOffsets(ii) = (ii-1)*ChannelWindowSize(2);
                    end
                else
                    for ii = 1:Nx
                        XchOffsets(ii,:) = (ii-1)*ChannelWindowSize(1);
                    end
                    for ii = 1:Ny
                        YchOffsets(:,ii) = (ii-1)*ChannelWindowSize(2);
                    end
                end
            end
            XchOffsets = XchOffsets';
            YchOffsets = YchOffsets';
            XchOffsets = XchOffsets(:)';
            YchOffsets = YchOffsets(:)';
            
            %offsets from master channel (micron)
            MCoffsetX=XchOffsets-XchOffsets(obj.MasterChannel);
            MCoffsetY=YchOffsets-YchOffsets(obj.MasterChannel);
            
            %desired first location in master channel (micron)
            %updated PJC 2011-08-27 center sampling in window
            %adds a buffer of BoxSize on each side fo sampling
            MCstartX=XchOffsets(obj.MasterChannel)+(ChannelWindowSize(1) - obj.xdist*obj.SamplingXY(1))/2;
            MCstartY=YchOffsets(obj.MasterChannel)+(ChannelWindowSize(2) - obj.ydist*obj.SamplingXY(2))/2;
            
            %distance to move stage for first position (micron)
            deltaX=MCstartX-obj.BeadXY(1);
            deltaY=MCstartY-obj.BeadXY(2);
            
            %desired first location of bead in all channels (micron)
            StartX=MCstartX+MCoffsetX;
            StartY=MCstartY+MCoffsetY;
            
            %Stage start position
            StageStartX=obj.MCLStartXYZ(1)+deltaX;
            StageStartY=obj.MCLStartXYZ(2)+deltaY;
            StageStartZ=obj.MCLStartXYZ(3);
            
            %setup saved variables
            obj.ImageStack=newim(obj.BoxSize,obj.BoxSize,obj.SamplingXY(1)*obj.SamplingXY(2),Nreg);
            obj.CoordinateArray=zeros(obj.SamplingXY(1)*obj.SamplingXY(2),Nreg,2);
            
            % initialize figure for displaying subregions
            h = figure('name','Channel Registration SubRegions');
            
            %loop over the scan positions
            for ny=0:obj.SamplingXY(2)-1
                for nx=0:obj.SamplingXY(1)-1
                    
                    %move stage
                    X=StageStartX+nx*obj.xdist;
                    Y=StageStartY+ny*obj.ydist;
                    Z=StageStartZ;
                    obj.SetStagePosition(obj,X,Y,Z);
                    
                    %bead locations (center of ROI) in pixels!
                    tmp=StartX+nx*obj.xdist;
                    Xroi=tmp(obj.ChannelRegistrationArray>0)/obj.PixelSize;
                    tmp=StartY+ny*obj.ydist;
                    Yroi=tmp(obj.ChannelRegistrationArray>0)/obj.PixelSize;
                    
                    %image and calculate sub-pixel locations
                    %                     [Startout images]=GetFitROI(obj,Xroi,Yroi);
                    [Startout images]=GetFitROImaxAdj(obj,Xroi,Yroi);
                    
                    %write images, locations into object
                    obj.ImageStack(:,:,ny*obj.SamplingXY(1)+nx,:)=images;
                    obj.BoxStartPixel(ny*obj.SamplingXY(1)+nx+1,:,:)=Startout;
                    obj.plotSubRegions(h)
                    
                end
            end
            
            %calculate positions
            obj.CalcPositions();
            
        end
        
        

        function obj = UserSelectBead(obj)
            %Get ROI and Bead Position
            
            %figure(obj.MCLhandle.figure1)
            %PauseForCamera
            close(findall(0,'name','focus'))
            disp('Find a field of view containing a bead and adjust to best focus')
            AndoriXon('ActionFocusButton_Callback',obj.andor_handles.ActionFocusButton,[],obj.andor_handles);
            im = squeeze(evalin('base','focus'));
            close(findall(0,'name','focus'))
            h = dipshow(im);
            diptruesize(h,200);
            dipmapping('log');
            udata = get(h,'userdata');
            udata.imname = 'Click on bead in Master Channel to be used for calibration';
            set(h,'name','Click on bead in Master Channel to be used for calibration','numbertitle','off',...
                'userdata',udata);
            % get initial coordinates for bead
            disp('Click on bead in Master Channel to be used for calibration')
            obj.BeadSelectXY = dipgetcoords(h,1);
            [obj.MCLStartXYZ(1) obj.MCLStartXYZ(2) obj.MCLStartXYZ(3)] = obj.GetStagePosition(obj);
            close(h)
            
            %check current stage position
            [Ny Nx] = size(obj.ChannelRegistrationArray);
            ChannelWindowSize = obj.PixelSize*obj.ImageSize./[Nx Ny]; %micron
%            if any((100-obj.MCLStartXYZ([1 2]) < ChannelWindowSize) | (obj.MCLStartXYZ([1 2]) < ChannelWindowSize))
%                error('ChannelRegistrationV1:MCLStageTooCloseToLimits',...
%                    'ChannelRegistrationV1: stage close to limits center and restart');
%            end
                     
            %Fit Data to Get Sub-Pixel Coordinates
            xstart=round(obj.BeadSelectXY(1)-obj.SearchBoxSize/2);
            ystart=round(obj.BeadSelectXY(2)-obj.SearchBoxSize/2);
            xend=xstart+obj.SearchBoxSize-1;
            yend=ystart+obj.SearchBoxSize-1;
            roi=im(xstart:xend,ystart:yend);
            %adj roi to center on max pixel
            [val xy] = max(roi);
            xstart1 = floor(xstart-obj.BoxSize/2+xy(1));
            ystart1 = floor(ystart-obj.BoxSize/2+xy(2));
            xend1=xstart1+obj.BoxSize-1;
            yend1=ystart1+obj.BoxSize-1;
            roi1=im(xstart1:xend1,ystart1:yend1);
            
            %identify initial MasterChannel from bead selection
            %             obj.MasterChannel = sub2ind([Ny Nx],...
            %                 floor(obj.BeadSelectXY(1)/(obj.ImageSize(1)/Ny))+1,...
            %                 floor(obj.BeadSelectXY(2)/(obj.ImageSize(2)/Nx))+1);
            obj.MasterChannel = sub2ind([Ny Nx],...
                floor(obj.BeadSelectXY(2)/(obj.ImageSize(2)/Ny))+1,...
                floor(obj.BeadSelectXY(1)/(obj.ImageSize(1)/Nx))+1);
            
            if ~obj.ChannelRegistrationArray(obj.MasterChannel)
                error('ChannelRegistrationV1:MasterNotInRegistrationArray', ...
                    ['ChannelRegistrationV1: Bead selected in channel that is ' ...
                    '0 in the ChannelRegistrationArray']);
            end
            
            
            [P]=GPUgaussMLEv2(permute(single(roi1-min(roi1)),[2 1 3]),obj.PSFsigma(obj.MasterChannel),30,1);
            %absolute position of bead on the ccd, nothing to do with stage
            %coordinates!
            obj.BeadXY=(P(1,1:2)+[xstart ystart])*obj.PixelSize;
        end
        
        function [Startout ims]=GetFitROI(obj,X,Y)
            %acquire and fit bead in ROI
            %global sequence %this just gets currrnt
            %X,Y are vectors with length equal to number of channels
            close(findall(0,'name','capture'))
            AndoriXon('ActionCaptureButton_Callback',obj.andor_handles.ActionCaptureButton,[],obj.andor_handles);
            im = evalin('base','capture');
            h=figure(1111);
            im=squeeze(sum( (im-obj.CCDoffset)/obj.CCDgain,[],3));
            dipshow(h,im);
            diptruesize(h,200);
            pause(.1);
            
            Nroi=length(X);
            ims=newim(obj.BoxSize,obj.BoxSize,Nroi);
            %Fit Data to Get Sub-Pixel Coordinates
            for ii=0:Nroi-1;
                
                xstart=round(X(ii+1)-obj.BoxSize/2);
                if xstart<0;xstart=0;end;
                ystart=round(Y(ii+1)-obj.BoxSize/2);
                if ystart<0;ystart=0;end;
                
                xend=xstart+obj.BoxSize-1;
                if xend>(size(im,1)-1);xend=(size(im,1)-1);xstart=xend-obj.BoxSize+1;end;
                yend=ystart+obj.BoxSize-1;
                if yend>(size(im,2)-1);yend=(size(im,2)-1);ystart=yend-obj.BoxSize+1;end;
                
                Startout(ii+1,:)=[xstart ystart];
                roi=im(xstart:xend,ystart:yend);
                ims(:,:,ii)=roi-min(roi); %CHECK FIX
            end
        end
        
        function [Startout ims]=GetFitROImaxAdj(obj,X,Y)
            %global sequence %this just gets currrnt
            %X,Y are vectors with length equal to number of channels
            close(findall(0,'name','capture'))
            AndoriXon('ActionCaptureButton_Callback',obj.andor_handles.ActionCaptureButton,[],obj.andor_handles);
            im = evalin('base','capture');
            %close capture figure
            close(findall(0,'name','capture'))
            %             h=figure(1111);
            im=squeeze(sum( (im-obj.CCDoffset)/obj.CCDgain,[],3));
            %             dipshow(h,im);
            %             diptruesize(h,200);
            pause(.1);
            
            Nroi=length(X);
            ims=newim(obj.BoxSize,obj.BoxSize,Nroi);
            %Fit Data to Get Sub-Pixel Coordinates
            for ii=0:Nroi-1;
                
                xstart=round(X(ii+1)-obj.SearchBoxSize/2);
                if xstart<0;xstart=0;end;
                ystart=round(Y(ii+1)-obj.SearchBoxSize/2);
                if ystart<0;ystart=0;end;
                
                xend=xstart+obj.SearchBoxSize-1;
                if xend>(size(im,1)-1);xend=(size(im,1)-1);xstart=xend-obj.SearchBoxSize+1;end;
                yend=ystart+obj.SearchBoxSize-1;
                if yend>(size(im,2)-1);yend=(size(im,2)-1);ystart=yend-obj.SearchBoxSize+1;end;
                
                %initial roi
                roi=im(xstart:xend,ystart:yend);
                %adj roi to center on max pixel
                [val xy] = max(roi);
                xstart1 = floor(xstart-obj.BoxSize/2+xy(1));
                ystart1 = floor(ystart-obj.BoxSize/2+xy(2));
                xend1=xstart1+obj.BoxSize-1;
                yend1=ystart1+obj.BoxSize-1;
                if xend1>(size(im,1)-1);xend1=(size(im,1)-1);xstart1=xend1-obj.BoxSize+1;end;
                if yend1>(size(im,2)-1);yend1=(size(im,2)-1);ystart1=yend1-obj.BoxSize+1;end;
                Startout(ii+1,:)=[xstart1 ystart1]; %make output
                roi1=im(xstart1:xend1,ystart1:yend1);
                ims(:,:,ii)=roi1-min(roi1); %CHECK FIX
            end
        end
        
        
        function obj=CalcPositions(obj)
            %calculates CoordinateArray using ImageStack and BoxStartPixel
            
            FourDsize=size(obj.ImageStack);
            linsize=[FourDsize(1) FourDsize(2) FourDsize(3)*FourDsize(4)];
            %prepare stack for GPUgaussMLE
            stack=reshape(obj.ImageStack,linsize);
            
            %CHECK convergence of GPUgaussMLEv2 for high bg cases
            [P]=GPUgaussMLEv2(permute(single(stack),[2 1 3]),obj.PSFsigma,30,1);
            
            %convert back to channels
            X=reshape(P(:,1),[FourDsize(3) FourDsize(4) 1]);
            Y=reshape(P(:,2),[FourDsize(3) FourDsize(4) 1]);
            
            obj.CoordinateArray(:,:,1)=X+obj.BoxStartPixel(:,:,1);
            obj.CoordinateArray(:,:,2)=Y+obj.BoxStartPixel(:,:,2);
            
        end
        
        
        function PixelSize = CalcPixelSize(obj)
            % CalcPixelSize     use the acquired data to determine pixelsize
            %
            % OUTPUT
            %   PixelSize - 1x2 vector with values [x y] where x is pixel
            %               height and y is pixel width in microns
            %
            %updated by Pat Cutler August 2, 2011
            %
            % NOTE: still needs to be tested for anything besides 2x2!!!!
            
            if isempty(obj.CoordinateArray)
                warning('This function requires a data set');
                return
            end
            
            fprintf('This function will not update the object''s PixelSize property!\n')
            [Ny Nx] = size(obj.ChannelRegistrationArray); % PJC 2011-07-21
            
            %X direction
            %Use first channel
            Xfound=squeeze(obj.CoordinateArray(:,:,1));
            %microns per pixel
            PixelSize(1) = obj.xdist/mean(mean(diff(reshape(Xfound,(obj.SamplingXY.*[1 Nx*Ny])))));
            
            %Y direction
            %Use first channel
            Yfound=squeeze(obj.CoordinateArray(:,:,2));
            Yscan = [];
            for ii = 1:obj.SamplingXY(2)
                Yscan = [Yscan Yfound(ii:obj.SamplingXY(1):end,:)];
            end
            %microns per pixel
            PixelSize(2) = obj.ydist/mean(mean(diff(Yscan)));
            %             Yscan=repmat(obj.ydist*(1:Ny),[1 Nx]);
            %
            %             Py=polyfit(Yscan,Yfound,2);
            %             YPixelsPerMicron=Py(1);
            %
            %             PixelSize=[1/XPixelsPerMicron 1/YPixelsPerMicron];
            
            %             figure
            %             plot(Xscan,Xfound,'bo','linewidth',2)
            %             plot(Xscan,Xscan*Px(2)+Px(1),'b-','linewidth',2)
            %             plot(Yscan,Yfound,'ro','linewidth',2)
            %             plot(Yscan,Yscan*Py(2)+Py(1),'r-','linewidth',2)
            %             xlabel('Relative Stage Position (micron)','textsize',16)
            %             ylabel('Found Position (Pixels)','textsize',16)
            %             legend('X Found Position','X Fit','Y found Position','Y Fit')
            
        end
        
        function plotSubRegions(obj,h)
            %PLOTSUBREGIONS   plot subregions for ChannelRegistration
            %
            % USAGE:
            %   plotSubRegions(obj,h)
            %   obj.plotSubRegions
            %   obj.plotSubRegions(h)
            %
            % INPUTS
            %   obj - ChannelRegistration object
            %   h - figure handle in which to place plot
            %
            %created by Pat Cutler August 27, 2011
            
            if ~exist('h','var')
                % Create figure
                h = figure;
            end
            im = newim(obj.ImageSize);
            boxIdx = (1:obj.BoxSize)-1;
            for ii = 1:size(obj.BoxStartPixel,1)
                for jj = 1:size(obj.BoxStartPixel,2)
                    im(obj.BoxStartPixel(ii,jj,1)+boxIdx,...
                        obj.BoxStartPixel(ii,jj,2)+boxIdx) = obj.ImageStack(:,:,ii-1,jj-1);
                end
            end
            dipshow(h,im);
            diptruesize(h,200);
            dipmapping('log');
            udata = get(h,'userdata');
            udata.imname = 'Channel Registration SubRegions';
            set(h,'name','Channel Registration SubRegions','numbertitle','off',...
                'userdata',udata);
            drawnow
        end
            
            
        function plot(obj,ax)
            %PLOT   plot localizations for ChannelRegistration
            %
            % USAGE:
            %   plot(obj,ax)
            %   obj.plot
            %   obj.plot(ax)
            %
            % INPUTS
            %   obj - ChannelRegistration object
            %   ax - axes handle in which to place plot
            %
            %updated by Pat Cutler August 2, 2011
            %
            % NOTE: still needs to be tested for anything besides 2x2!!!!
            
            if ~exist('ax','var')
                % Create figure
                figure1 = figure('name',datestr(now,'dd-mmm-yyyy HH:MM:SS'));
                % Create axes
                ax = axes('Parent',figure1);
            end
            set(ax,'YDir','reverse','PlotBoxAspectRatio',[obj.ImageSize/max(obj.ImageSize) 1])
            % Uncomment the following line to preserve the X-limits of the axes
            xlim(ax,[0 obj.ImageSize(1)-1]);
            % Uncomment the following line to preserve the Y-limits of the axes
            ylim(ax,[0 obj.ImageSize(2)-1]);
            hold all
            plot(obj.CoordinateArray(obj.trnLogical,:,1),obj.CoordinateArray(obj.trnLogical,:,2),'bo')
            plot(obj.CoordinateArray(~obj.trnLogical,:,1),obj.CoordinateArray(~obj.trnLogical,:,2),'ro')
            nRowCols = size(obj.ChannelRegistrationArray);
            for ii = 1:nRowCols(1)-1 %loop through rows
                x = [0 obj.ImageSize(1)-1];
                y = ii*repmat(obj.ImageSize(2),[1 2])/nRowCols(1);
                plot(x,y,'color','black','linewidth',2)
            end
            for ii = 1:nRowCols(2)-1 %loop through cols
                y = [0 obj.ImageSize(2)-1];
                x = ii*repmat(obj.ImageSize(1),[1 2])/nRowCols(2);
                plot(x,y,'color','black','linewidth',2)
            end
            hold off
            xlabel('x (pixels)')
            ylabel('y (pixels)')
        end
        
        function h = plotTform(obj)
            %PLOTTFORM   plot transformation results for ChannelRegistration
            %
            % Makes a subPlot for each channel in the same arrangement as the
            % channels are spread on camera.
            % Plot master channel vs shifted channel vs transformed channel
            % for each transformed channel
            % For master channel subplot use plot method
            % Training points are indicated in blue and test points are indicated in red 
            %
            % USAGE:
            %   plot(obj)
            %   obj.plot
            %
            % INPUTS
            %   obj - ChannelRegistration object
            % OUPUTS
            %   h - figure handle
            %
            %updated by Pat Cutler August 3, 2011
            %
            % NOTE: still needs to be tested for anything besides 2x2!!!!
            
            %check that transform has been performed
            if obj.NoCalcTform
                error('ChannelRegistation:plotTformNoCalcTform', 'ChannelRegistation: calcTform must be performed prior to plotTform')
            end
            % Create figure
            if obj.timeStamp == obj.tstObj.timeStamp
                figName = ['ChannelRegistation(trn&tst' datestr(obj.timeStamp,'yyyy-mm-dd-HH-MM-SS') ')'];
            else
                figName = ['ChannelRegistation(trn' datestr(obj.timeStamp,'yyyy-mm-dd-HH-MM-SS')...
                ':tst' datestr(obj.tstObj.timeStamp,'yyyy-mm-dd-HH-MM-SS') ')'];
            end
            h = figure('name',figName);
            nRowCols = size(obj.ChannelRegistrationArray); %get number of rows and columns of channels
            subSize = obj.ImageSize./nRowCols([2 1]); %get image size per channel
            nChannels = prod(nRowCols); %get total number of channels
            %get the row and column of the master channel
            masterRowCol = [rem(obj.MasterChannel,nRowCols(1)) ceil(obj.MasterChannel/nRowCols(1))];
            masterRowCol(masterRowCol == 0) = nRowCols(1);
            for ii = 1:nChannels %loop through channels
                if obj.ChannelRegistrationArray(ii) %check that registration data is available
                    % Create axes
                    chRowCol = [rem(ii,nRowCols(1)) ceil(ii/nRowCols(1))];
                    chRowCol(chRowCol == 0) = nRowCols(1);
                    ax = subplot(nRowCols(1),nRowCols(2),chRowCol(2)+nRowCols(2)*(chRowCol(1)-1),'Parent',h);
                    if ii == obj.MasterChannel %plot all channels in master channel axes
                        obj.plot(ax)
                        title(ax,{['Channel ' num2str(ii) ' is the Master Channel'];...
                            'trn ''b'' tst ''r'''})
                    else
                        plot(ax,obj.tform(ii).tst.Measured(:,1),obj.tform(ii).tst.Measured(:,2),'ro',...
                            obj.tform(ii).tst.posShift(:,1),obj.tform(ii).tst.posShift(:,2),'rs',...
                            obj.tform(ii).tst.pos(:,1),obj.tform(ii).tst.pos(:,2),'rx')
                        hold(ax,'on')
                        plot(ax,obj.tform(ii).Measured(:,1),obj.tform(ii).Measured(:,2),'bo',...
                            obj.tform(ii).trn.posShift(:,1),obj.tform(ii).trn.posShift(:,2),'bs',...
                            obj.tform(ii).trn.pos(:,1),obj.tform(ii).trn.pos(:,2),'bx')
                        quiver(ax,obj.tform(ii).tst.Measured(:,1),obj.tform(ii).tst.Measured(:,2),...
                            obj.tform(ii).tst.xdiff,obj.tform(ii).tst.ydiff,.25,'k','linewidth',2)
                        hold(ax,'off')
                        % set axes properties
                        set(ax,'YDir','reverse','PlotBoxAspectRatio',[subSize/max(subSize) 1])
                        xlim(ax,subSize(1)*[masterRowCol(2)-1 masterRowCol(2)]);
                        ylim(ax,subSize(2)*[masterRowCol(1)-1 masterRowCol(1)]);
                        xlabel('x (pixels)')
                        ylabel('y (pixels)')
                        if any(~obj.trnLogical)
                            shiftError = obj.tform(ii).tst.shift;
                            fre = obj.tform(ii).tst.fre;
                        else
                            shiftError = obj.tform(ii).trn.shift;
                            fre = obj.tform(ii).trn.fre;
                        end
                        title({['Channel ' num2str(ii) ' transform'];...
                            'master ''o'' shift ''s'' tform ''x''';...
                            [num2str(size(obj.CoordinateArray,1)) ' control points & ' ...
                            num2str(obj.tform(ii).tst.dir.count) ' at ' ...
                            num2str(obj.tform(ii).tst.dir.min) '-'  ...
                            num2str(obj.tform(ii).tst.dir.max)];...
                            ['shift error pixels(' num2str(shiftError)...
                            ') nm(' num2str(shiftError*obj.PixelSize*1e3) ')'];...
                            ['tform error pixels(' num2str(fre)...
                            ') nm(' num2str(fre*obj.PixelSize*1e3) ')']},...
                            'fontsize',8)
                    end
                end
            end
        end
        
        
        function calcTform(obj)
            %CALCTFORM  calculate transform into master channel
            %
            % USAGE:
            %   calcTform(obj)
            %   obj.calcTform
            %
            % ---- tform property description ----
            % Structured array with size 1 by the number of channels
            % containing results from the image transformation.
            %   Measured - measured position for points in master channel
            %   Sensed - sensed position for points in transformed channel
            %   pX - matrix of coefficients for each local weighted fitting in X.
            %        Dimensions are number of sensed points by C.
            %   pY - matrix of coefficients for each local weighted fitting in Y.
            %        Dimensions are number of sensed points by C.
            %   svec - shift vector
            %   trn - structure containing training information.
            %   tst - structure containing testing information.
            %   
            % ---- fields of trn/tst ----
            %     pos - transformed position for trn/tst points to master channel position
            %     posShift - shifted position for trn/tst points to master channel position
            %     shift - standard error in pixels after simple shift
            %     fre - fiducial registration error in pixels. Standard error between
            %           measured positions in the master channel with
            %           sensed positions in the transformed channel.
            %     xdiff - x difference between transformed position and
            %             position in MasterChannel
            %     xdiff - y difference between transformed position and
            %             position in MasterChannel
            %     dir - structure with information about correlated
            %           directional error. Contains fields:
            %               count - max number of measurements in a single
            %                       directional bin. Total of 10 bins with
            %                       36 degrees in each bin.
            %               min - min value in degrees for bin with
            %                     corresponding to count field
            %               max - max value in degrees for bin with
            %                     corresponding to count field
            %
            %Created by Pat Cutler August 2011
            
            % identify registered channels
            registeredChannels = find(obj.ChannelRegistrationArray);
            registeredChannels = registeredChannels(:)';
            % adjust master channel index for Channel Registration Array
            adjMasterChannel = find(registeredChannels == obj.MasterChannel);
            M = squeeze(obj.CoordinateArray(:,adjMasterChannel,:));
            %             registeredChannels(adjMasterChannel) = []; %channels to register
            %             imageArrayIdx = obj.ChannelRegistrationArray; %channels to register
            %             imageArrayIdx(obj.MasterChannel) = 0;
            nPoints = sum(obj.trnLogical(:));
            for ii = registeredChannels %loop through channels
                if ii ~= adjMasterChannel
                    idx = find(registeredChannels == ii);
                    S = squeeze(obj.CoordinateArray(:,idx,:));
                    %                     tform = CalcLWtform(M,S,obj.N,obj.C);
                    tform = obj.CalcLWtform(M(obj.trnLogical(:),:),S(obj.trnLogical(:),:),obj.N,obj.C);
                    tform.trn.pos = obj.LWtform(S(obj.trnLogical(:),:),tform);
                    tform.trn.posShift = tform.Sensed+repmat(tform.svec,[nPoints 1]);
                    %tform.trn.shift = sqrt(sum(sum((tform.trn.posShift-tform.Measured).^2,2))/(2*nPoints));
                    tform.trn.shift = sqrt(sum(sum((tform.trn.posShift-tform.Measured).^2,2)/nPoints));
                    tform.trn.fre = sqrt(sum(sum((tform.trn.pos-tform.Measured).^2,2)/nPoints));
                    tform.trn.xdiff = tform.trn.pos(:,1)-tform.Measured(:,1);
                    tform.trn.ydiff = tform.trn.pos(:,2)-tform.Measured(:,2);
                    theta = cart2pol(tform.trn.xdiff,-tform.trn.ydiff);
                    [tout,rout] = rose(theta,10);
                    deg = tout(2:4:end)*180/pi;
                    tform.trn.dir.count = max(rout(2:4:end));
                    tform.trn.dir.min = deg(find(max(rout(2:4:end)) == rout(2:4:end),1));
                    tform.trn.dir.max = tform.trn.dir.min+36;
                    tform.tst = [];
                    obj.tform(ii) = tform;                  
                end
            end
        end
        
        
        function tstTform(obj,tstObj)
            %TSTTFORM  test transform into master channel
            %
            % USAGE:
            %   tstTform(obj)
            %   obj.tstTform
            %
            %
            %Created by Pat Cutler August 29, 2011
            
            %check that transform has been performed
            if obj.NoCalcTform
                error('ChannelRegistrationV1:tstTformNoCalcTform', 'ChannelRegistrationV1: calcTform must be performed prior to tstTform')
            end
            %use training set for 
            if nargin < 2 || isempty(tstObj)
                tstObj = obj;
            end
            obj.tstObj = ChannelRegistrationV1(obj);
            copyobj(tstObj,obj.tstObj);
            % identify registered channels
            registeredChannels = find(obj.ChannelRegistrationArray);
            registeredChannels = registeredChannels(:)';
            % adjust master channel index for Channel Registration Array
            adjMasterChannel = find(registeredChannels == obj.MasterChannel);
            M = squeeze(tstObj.CoordinateArray(:,adjMasterChannel,:));
            nPoints = sum(~tstObj.trnLogical(:));
            for ii = registeredChannels %loop through channels
                if ii ~= adjMasterChannel
                    tform = obj.tform(ii);
                    idx = find(registeredChannels == ii);
                    S = squeeze(tstObj.CoordinateArray(:,idx,:));
                    tform.tst.pos = obj.LWtform(S(~tstObj.trnLogical(:),:),tform);
                    tform.tst.Measured = M(~tstObj.trnLogical(:),:);
                    tform.tst.Sensed = S(~tstObj.trnLogical(:),:);
                    tform.tst.posShift = tform.tst.Sensed+repmat(tform.svec,[nPoints 1]);
                    tform.tst.shift = sqrt(sum(sum((tform.tst.posShift-tform.tst.Measured).^2,2)/nPoints));
                    tform.tst.fre = sqrt(sum(sum((tform.tst.pos-tform.tst.Measured).^2,2)/nPoints));
                    tform.tst.xdiff = tform.tst.pos(:,1)-tform.tst.Measured(:,1);
                    tform.tst.ydiff = tform.tst.pos(:,2)-tform.tst.Measured(:,2);
                    theta = cart2pol(tform.tst.xdiff,-tform.tst.ydiff);
                    [tout,rout] = rose(theta,10);
                    deg = tout(2:4:end)*180/pi;
                    tform.tst.dir.count = max(rout(2:4:end));
                    tform.tst.dir.min = deg(find(max(rout(2:4:end)) == rout(2:4:end),1));
                    tform.tst.dir.max = tform.tst.dir.min+36;
                    obj.tform(ii) = tform;
                end
            end
        end
        
        function locMaster = locTform(obj,ch,locCh)
            %LOCTFORM  transform localizations from ch into MasterChannel
            %
            % USAGE:
            %   locMaster = locTform(obj,ch,locCh)
            %   locMaster = obj.locTform(ch,locCh)
            %
            % INPUTS
            %   obj - ChannelRegistration object
            %   ch - channel in which localizations performed. Example of
            %        channel indexing for 2x2 image [1 3;2 4].
            %   locCh - localizations in absolute CCD pixels related to obj.ImageSize
            % OUTPUT
            %   locMaster - localizations in MasterChannel coordinates
            %
            %Created by Pat Cutler August 31, 2011
            
            %error checking
            if nargin ~= 3
                error('ChannelRegistrationV1:locTformNargin',...
                    'ChannelRegistrationV1: locTform requires ch and loc inputs')
            end 
            if ~obj.ChannelRegistrationArray(ch)
                error('ChannelRegistrationV1:locTformChNotChannelRegistrationArray',...
                    'ChannelRegistrationV1: ch must be in ChannelRegistrationArray')
            end
            %check that transform has been performed
            if obj.NoCalcTform
                error('ChannelRegistation:locTformNoCalcTform', 'ChannelRegistation: calcTform must be performed prior to locTform')
            end
            
            if ch == obj.MasterChannel
                locMaster = locCh; %return input if ch is same as MasterChannel
            else
                %transform into MasterChannel
                locMaster = obj.LWtform(locCh,obj.tform(ch));
            end
        end
        
        function check = NoCalcTform(obj)
            %CHECKCALCTFORM     check that calcTform has been performed
            tmp = obj.ChannelRegistrationArray;
            tmp(obj.MasterChannel) = 0;
            check = isempty(obj.tform(find(tmp,1)).Measured);
        end
        
        
        
       
        
                
        function set.ChannelRegistrationArray(obj,value)
            if nnz((value == 1) + (value == 0))~=numel(value)
                error('ChannelRegistrationV1:ChannelRegistrationArray','ChannelRegistrationArray must be an array of 0 or 1''s');
            end
            obj.ChannelRegistrationArray = value;
        end
        
        function set.xdist(obj,value)
            if (value <= 0)
                error ('ChannelRegistrationV1:xdist','xdist must be a positive value');
            end
            if (value > 100)
                error ('ChannelRegistrationV1:xdist','xdist must be less than the travel range of the stage');
            end
             obj.xdist = value;
        end
        
        function set.ydist(obj,value)
            if (value <= 0)
                error ('ChannelRegistrationV1:ydist','ydist must be a positive value');
            end
            if (value > 100)
                error ('ChannelRegistrationV1:ydist','ydist must be less than the travel range of the stage');
            end
            obj.ydist = value;
        end
        
    end
    
    
    
    methods(Static)
        
         function out = LWtform(X,tform,N)
            %LWtform Locally weighted transformation of coordinates using tform
            %   X: input cooridates [x y]
            %   tform: tform calculated from CalcLWtform.
            %   N: number of nearest neighbors to be used for transform. Must be
            %   greater than C. Default value is C+1.
            %   out: transformed coordinates
            %
            %See Also:
            %   CalcLWtform
            
            if nargin == 1 & ischar(X) & strcmp(X,'DIP_GetParamList')
                out = struct('menu','none');
                return
            end
            
            
            out=zeros(size(X));
            Npoints=size(X,1);
            pX=tform.pX;
            pY=tform.pY;
            S=tform.Sensed;
            
            C = size(pX,2);
            if ~exist('N','var') || isempty(N)
                N = C+1; % number of nearest neighbors
            end
            
            %             n=size(pX,2);
            
            for pp=1:Npoints
                if mod(pp,10000) == 0
                    fprintf('\n %d out of %d has been done\n',pp,Npoints);
                end
                %calculate the distance from the point to the n-1 nearest sensed points
                tmp=repmat(X(pp,:),[size(S,1) 1]);
                D=sqrt(sum( (S-tmp).^2,2));
                %sort these into distance and keep the IDs
                [Dsort IDs]=sort(D);
                %the IDs of the n closest points are then
                IDclose=IDs(1:N);
                R=Dsort(1:N)/Dsort(N); %R ranges from >0 to 1 (last value used to scale weighting)
                W=1-3*R.^2+2*R.^2; %last value of W is always 0
                %last point only used to scale weighting but weight is zero i.e. (1-3(1)^2+2*(1)^2 = 0
                
                %calculate the weighted correction
                x2=repmat(X(pp,1),[N 1]);
                y2=repmat(X(pp,2),[N 1]);
                DM=[x2*0+1 x2 y2 x2.*y2 x2.^2 y2.^2];
                DM=DM(:,1:C);
                tpX=pX(IDclose,:);
                tpY=pY(IDclose,:);
                out(pp,1)=sum(sum(DM.*tpX,2).*W)/sum(W)+x2(1);
                out(pp,2)=sum(sum(DM.*tpY,2).*W)/sum(W)+y2(1);
                
            end
            
         end
        
        function tform = CalcLWtform(M,S,N,C)
            %CalcLWtform Calculate Locally Weighted tform
            %   M: Set of measured points (channel 1)
            %   S: Set of sensed points   (channel 2)
            %   N: number of nearest neighbors to be used for transform. Must be
            %   greater than C. Default value is C+1.
            %   C: order of polynomial to fit. Allowable values 1-6. Default 4.
            %
            %   NOTES:
            %A polynomial is fit for every control point using n=6 local points
            %1 x y x*y x^2 xy y^2
            %The polynomial will give the correction required for the second channel.
            %There are two sets of polynomials, that for X and Y corrections.
            %Coordintes are absolute coordinates (0,0) at top left as given by GPUSPT
            %svec can be used for shifting images but is not used in the transform
            %
            % Created by Keith Lidke UNM
            % Updated by Pat Cutler 2011_07_21
            
            % Avoid being in dipimage menu
            if nargin == 1 & ischar(M) & strcmp(M,'DIP_GetParamList')
                tform = struct('menu','none');
                return
            end
            
            % define default parameters
            if nargin < 4 || isempty(N)
                C = 4; % order of polynomial
            end
            if nargin < 3 || isempty(N)
                N = C+1; % number of nearest neighbors
            end
            
            if N+1 < C
                error('CalcLWtform:NlessC','CalcLWtform: input N must be less than C')
            end
            
            % n=3; % number of nearest neighbors for fit and order of polynomial
            tform = struct('Measured',[],'Sensed',[],'pX',[],'pY',[],'svec',[],'trn',[],'tst',[]);
            Npoints=size(M,1);
            if N+1 > Npoints
                error('CalcLWtform:NgreaterNpoints','CalcLWtform: M must have at least as many rows as input N+1')
            end
            pX=zeros(Npoints,C);
            pY=zeros(Npoints,C);
            tform.svec=mean(M-S,1);
            
            
            for pp=1:Npoints
                %calculate the distance from the point to the n-1 nearest points
                tmp=repmat(M(pp,:),[Npoints 1]);
                D=sqrt(sum( (M-tmp).^2,2));
                %sort these into distance and keep the IDs
                [Dsort IDs]=sort(D);
                %the IDs of the n closest points are then
                IDclose=IDs(1:N+1);
                x2=S(IDclose,1);
                y2=S(IDclose,2);
                
                % x-corrections
                Xc=M(IDclose,1)-S(IDclose,1);
                %'design matrix'
                DM=[x2*0+1 x2 y2 x2.*y2 x2.^2 y2.^2];
                DM=DM(:,1:C);
                pX(pp,:)=inv(DM'*DM)*DM'*Xc;
                %test the coorections:
                %     tmp=(DM*pX(pp,:)'+S(IDclose,1))-M(IDclose,1);
                % y-corrections
                Yc=M(IDclose,2)-S(IDclose,2);
                %'design matrix'
                pY(pp,:)=inv(DM'*DM)*DM'*Yc;
                %test the coorections:
                %     tmp=(DM*pY(pp,:)'+S(IDclose,2))-M(IDclose,2);
            end
            
            tform.Measured=M;
            tform.Sensed=S;
            tform.pX=pX;
            tform.pY=pY;
            
            %svec can be used for shifting images, but is not used in the transform
            
            tform.svec=mean(M-S,1);
            
        end
        
        function [X Y Z] = GetStagePosition(obj)
            switch obj.StageID
                case 'Prior'
                    
                case 'MCL'
                    [X Y Z] = GetMCLNanoDrive();
            end
            
        end
        
        function SetStagePosition(obj,X,Y,Z)
            switch obj.StageID
                case 'Prior'
                    
                case 'MCL'
                    SetMCLNanoDrive(X,Y,Z);
            end
        end
        
        function versn()
            global VERSION
            r = '$Rev: 502 $';
            d = '$Date: 2013-06-21 11:56:44 -0600 (Fri, 21 Jun 2013) $';
            a = '$Author: sheng $';
            h = '$HeadURL: https://abbe.phys.unm.edu/svn/MATLAB/Instrumentation/ChannelRegistrationV1.m $';
            eval(['VERSION.' mfilename('class') '.Rev=''' r(7:end-2) ''';']);
            eval(['VERSION.' mfilename('class') '.Date=''' d(8:end-2) ''';']);
            eval(['VERSION.' mfilename('class') '.Author=''' a(10:end-2) ''';']);
            eval(['VERSION.' mfilename('class') '.HeadURL=''' h(11:end-2) ''';']);
        end
    end
end





