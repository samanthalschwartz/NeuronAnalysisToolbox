% TrackMovie.m
%
% A Class for generating and viewing track movies, which are tracks and localizations rendered on top
% of an image that can be animated and rendered as a movie or stepped through manually, etc.  This class
% is intended tro be future compatable with the new matlab HG2 graphics system and should work on 
% Matlab 2014b+.  
%
% Mark J. Olah (mjo@cs.unm.edu)
% 04 - 2015
%
% NOTE: This is inspired by the features of Pat Cutler's dipTrack, but we are using matlab's image plotting
% to remove any depenedency on dipimage and the undocumented features of dipimage that dipTrack relies
% on.  Also the class based interface makes it much easier to provide many interfaces and options in a
% more organized way.
%
%

classdef TrackMovie < handle
    properties (Constant=true)
        TRACK_COLOR_METHODS={'Sequence', 'Speed', 'Temporal', 'Wavelength'};
    end
    properties
        frames;
        sumImage;
        tracks; % Tracks are in physical aka. world coordinates.
        ROIphysical; % The real coordinate the image spans [xmin xmax ymin ymax]. In physical (aka. world coordinates)
        
        frameT; % This is computed from the tracks
        tBounds; % time bounds [minT maxT] this is computed from the tracks
        frameBounds; % time bounds [minFrameIdx maxFrameIdx] this is computed from the tracks must not be bigger than size(obj.frames,3)
        imRef; % an imref2d object representing this image coordiate systems


        FigureSize=[1080, 720];
        MovieSize=[540, 720];
        FrameLag=-1;
        TrackLineWidth=1;        
        TrackVertexMarker='none';
        TrackColorMethod = 'Sequence'; 
        TrackColorMap = @jet;
        TrackColorRange = 256;  %Number of different discrete colors in the tracks
        ImageColorMap = @gray;
        ImageColorRange = 256;  %Number of different discrete colors in the image
        ImageGlobalNormalize = false; % Set to falst to normalize each frame, set to true to normalize all frames globally.
        ImageAlpha = 1.0; 
        ColorbarOn = 0;
        ColorbarLabel = 'Sequence';
        FigureBoarderSize = [10 20 10 30]; %[left bottom right top]  set these to positive numbers to make a boarder between figure and axes
        TitleOn = 1;
        AxesFontSize = 14;
        TitleFontSize = 14;
        DrawLocalizations = 1; %Set to 0 to prevent drawing localizations
        DrawBackgroundFrames = 1; %Set to 0 to prevent drawing any background frames
        LocalizationSigmaMultiple = 3; % The size of the localization marker in sigmas (2sigma=95%)
        LocalizationSphereSize = 20; %This is the number of gradation levels draw in the sphere

        VideoWriterProfile = 'Archival';
        VideoWriterFrameRate = 25;

        trackCs; % Track Colors.  This is set by obj.setTrackColorMethod()
    end

    properties (Transient=true)
        data;
        video;
        figureH;
        axH;
        colorbarH;
    end

    properties (Access=protected)
%         frameMax; %maximum scaling level used in normalization
        trackMin;
        trackMax; 
        cMaxVal;
    end
    methods
        function obj = TrackMovie(varargin)
            % [case1]: Load from rpt
            %   tm = TrackMovie(rpt, tracks)
            %   rpt - RPT object or filename
            %   tracks - [optional] an RPT format tracks cell-array to use.  
            %            [default=use all tracks in RPT]
            % [case2]:
            %   tm = TrackMovie(tracks, frames, ROIphysical, frameBounds)
            % tracks - A RPT tracks cell-array formated set of tracks
            % frames - The image to display.  
            % ROIphysical - [optional] Physical region of space spaned by image [xmin xmax ymin ymax]
            %               Note: we use the matlab image coordinates where rows=y cols=x
            %               [default = [0 size(2) 0 size(1)]
            % frameBounds - [optional] The frame bounds for the tracks movie format: [minFrameIdx, maxFrameIdx]
            %               This should be set to make sure the added movie (frames) has the correct length to
            %               correspond to the tracks.  The tracks frameIdx column is always in absolute coordinates
            %               but the frames(:,:,1) frame will be displayed corresponding to the obj.framesBounds(1)
            %               index.
            %               Note: to retrieve from RPT, use rpt.ROI(5:6).
            %               [default = Use minimum and maximum frame index from tracks.]
            if nargin>0
                obj.initialize(varargin{:});
            end
        end

        function initialize(obj,varargin)
            %Arguments the same as for the constructor
            if nargin<1
                error('TrackMovie:initialize','No arguments given');
            elseif isa(varargin{1},'RPT')
                obj.loadRPT(varargin{:}); %Load from RPT object
            elseif ischar(varargin{1})
                fname=varargin{1}; %Load from .rpt file
                [~,~,ext]=fileparts(fname);
                switch ext
                    case '.rpt'
                        obj.loadRPT(RPT(fname), varargin{2:end});
                    otherwise
                        error('TrackMovie:initialize','Unknown file type: %s',ext);
                end
            elseif iscell(varargin{1})
                obj.loadTrackCellArray(varargin{:}); %Direct load from cell-array of tracks
            end
            %Common initializtion
            assert(obj.frameBounds(1)>=1 && obj.frameBounds(2)-obj.frameBounds(1)+1==size(obj.frames,3));
            obj.tBounds = (obj.frameBounds-1)*obj.frameT;
        end

         function setTrackColorMethod(obj, methodName, trackColorMap)
            if nargin == 3
                obj.TrackColorMap = trackColorMap;
            else 
                obj.TrackColorMap = @jet;
            end
            
            switch methodName
                case 'Sequence'
                    obj.trackCs = cellmap(@(i) i*ones(size(obj.tracks{i},1),1), 1:length(obj.tracks));
                    obj.ColorbarLabel = 'Sequence';
                    obj.ColorbarOn = length(obj.trackCs)>1;
                    obj.TrackColorRange = length(obj.tracks);
                case 'TrackLength'
                    obj.trackCs = cellmap(@(T) round((T(end,1)-T(1,1))/obj.frameT)*ones(size(T,1),1), obj.tracks);
                    obj.ColorbarLabel = 'Track Length (Frames)';
                    obj.ColorbarOn = length(obj.trackCs)>1;
                    obj.TrackColorRange = max(cellfun(@max, obj.trackCs));
                case 'Speed'
                    obj.trackCs = cellmap(@(T) TrackSegmentAnalysis.estimateSpeedWindow(T,6), obj.tracks); 
                    obj.ColorbarLabel = 'Speed ($\mu^2/\mathrm{s}$)';
                    obj.ColorbarOn = 1;
                    obj.TrackColorRange = 256;
                case 'Temporal'
                    obj.trackCs = cellmap(@(T) T(:,1), obj.tracks);
                    obj.ColorbarLabel = 'Time (s)';
                    obj.ColorbarOn = 1;
                    obj.TrackColorRange = 256;
                case 'CsInput'
                    obj.ColorbarLabel = 'Unique';
                otherwise
                    error('TrackMovie:setTrackColorMethod','Unknown coloring method: %s', methodName);
            end
            obj.TrackColorMethod = methodName;
            [obj.trackCs, obj.trackMin, obj.trackMax] = RPT.normalizeTrackColors(obj.trackCs, obj.TrackColorRange, obj.ImageColorRange);
        end
        
        function h=viewSequence(obj, frameIdx)
            % frame - [optional] initial frame
            if nargin==1
                frameIdx = obj.frameBounds(1);
            end
            obj.figureH=figure();
            obj.figureH.Units='Pixels';
            obj.figureH.Position=[10 10 obj.FigureSize];
            obj.axH=obj.initializeAxes();
            obj.axH.Units='Pixels';
            obj.figureH.SizeChangedFcn = @sizeChanged_CB;
            imH=obj.plotBackgroundFrame(frameIdx);
            trackHs=obj.plotTracksFrame(frameIdx);
            locHs=obj.plotLocsFrame(frameIdx);
            if obj.TitleOn
                title(sprintf('Frame:%i Time %.3f (s)',frameIdx,(frameIdx-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
            end
            sliderH = uicontrol('Style','slider','Min',obj.frameBounds(1), 'Max',obj.frameBounds(2),'Value',frameIdx,...
                                'SliderStep',[1/(diff(obj.frameBounds)) 10/(diff(obj.frameBounds))],...
                                'Units','normalized','Position',[0.05, 0, 0.8, 0.035], 'Callback',@slider_CB);
            whitebg(obj.figureH);
%             set(findobj(obj.figureH,'Type','text'),'Color',[0 1 0])
            drawnow();
            
            h = obj.figureH;
            function slider_CB(hObj,~)
                idx=round(hObj.Value);
                delete(imH);
                delete(trackHs);
                delete(locHs);
                imH=obj.plotBackgroundFrame(idx);
                if obj.FrameLag>0
                    trackHs=obj.plotTracksFrame(idx,idx-obj.FrameLag);
                else
                    trackHs=obj.plotTracksFrame(idx);
                end
                locHs=obj.plotLocsFrame(idx);
                if obj.TitleOn
                    title(sprintf('Frame:%i Time %.3f (s)',idx,(idx-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
                end
            end
            function sizeChanged_CB(~,~)
                f = gcbo;
                pos = get(f,'Position');
                pos(1:2) = 0;
                buf = [ 10 50 10 10];
                obj.axH.ActivePositionProperty='OuterPosition';
                GUIBuilder.positionAxes(obj.axH,pos,buf); 
            end
        end

        function vw = initializeVideoSequence(obj, filename)
            
            obj.video = VideoWriter(filename, obj.VideoWriterProfile);
            obj.video.FrameRate = obj.VideoWriterFrameRate;
            obj.video.open();
            obj.figureH=figure('Units','Pixels','Position',[0 0, obj.MovieSize]);
            obj.axH=obj.initializeAxes();
            obj.axH.Units='pixels';
            boarder= obj.FigureBoarderSize;
            if obj.ColorbarOn %adjust right figure boarder if color bar is used
                boarder(3) = boarder(3) + 70;
            end
            GUIBuilder.positionAxes(obj.axH, [0, 0, obj.MovieSize], boarder); 
            whitebg(obj.figureH);

            vw = obj.video;
        end

        function finalizeVideoSequence(obj)
            obj.video.close();
            close(obj.figureH);
        end

        function configureVideoZoom(obj)
            zoom();
            disp('Press any key to stop zoom configuration');
            while ~waitforbuttonpress()
            end
        end

        function renderActiveTracks(obj, frameIdxs)
            if nargin<2
                frameIdxs = obj.frameBounds(1):obj.frameBounds(2);
            end
            lag=obj.FrameLag;
            for i=frameIdxs
                imH=obj.plotBackgroundFrame(i);
                if lag>0
                    trackHs=obj.plotTracksFrame(i,i-lag);
                else
                    trackHs=obj.plotTracksFrame(i);
                end
                locHs=obj.plotLocsFrame(i);
                if obj.TitleOn
                    title(sprintf('Frame:%i Time %.3f (s)',i,(i-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
                end
                drawnow();
                F=getframe(obj.figureH);
                obj.video.writeVideo(F);
                delete(imH);
                delete(trackHs);
                delete(locHs);
            end
        end

        function renderStaticTracks(obj)
            obj.plotBackgroundFrame(obj.sumImage);
            obj.plotTracksFrame();
        end

        function renderStaticFlyover(obj, nframes)
            if nargin<2
                nframes=360;
            end
            obj.renderStaticTracks();
            nVPans=ceil(nframes/8);
            nHPans=ceil(nframes*5/8);
            view(0,90);
            for i=1:nVPans
                F=getframe(obj.figureH);
                obj.video.writeVideo(F);
            end
            for i=1:nVPans
                camorbit(0,-90/nVPans,'camera')
                drawnow
                F=getframe(obj.figureH);
                obj.video.writeVideo(F);
            end
            for i=1:nHPans
                camorbit(360/nHPans,0,'camera')
                drawnow
                F=getframe(obj.figureH);
                obj.video.writeVideo(F);
            end
            for i=1:nVPans
                camorbit(0,90/nVPans,'camera')
                drawnow
                F=getframe(obj.figureH);
                obj.video.writeVideo(F);
            end           
            view(0,90);
            for i=1:nVPans
                F=getframe(obj.figureH);
                obj.video.writeVideo(F);
            end
        end

        function initializeColorMap(obj)
            colormap([obj.ImageColorMap(obj.ImageColorRange); obj.TrackColorMap(obj.TrackColorRange)]);
        end

        function ax=initializeAxes(obj,ax)
            if nargin==1
                ax=gca();
            end
            cla(ax,'reset');

            %setup axes
            ax.SortMethod='ChildOrder';
            ax.YDir='reverse';
            ax.TickDir='out';
            ax.XMinorTick='on';
            ax.YMinorTick='on';
            ax.ZMinorTick='on';
            ax.XGrid='on';
            ax.YGrid='on';
            ax.ZGrid='on';
            ax.Box='on';
            ax.BoxStyle='full';
            ax.Projection='Orthographic';
            obj.initializeColorMap();
            caxis([0 obj.TrackColorRange]);
            if obj.ColorbarOn
                obj.colorbarH = RPT.configureTracksColorBar(obj.TrackColorRange, obj.ImageColorRange, [obj.trackMin, obj.trackMax], obj.ColorbarLabel);
            else
                colorbar('off');
            end            

            %Set aspect ratio of axes
            axis([obj.ROIphysical (obj.frameBounds-1)*obj.frameT]);
            asp(1) = (obj.ROIphysical(2)-obj.ROIphysical(1))/(obj.ROIphysical(4)-obj.ROIphysical(3));
            asp(2) = 1;
            asp(3) = mean([asp(1),asp(2)]);
            pbaspect(asp);
            %Fix grids and ticks
            view(0,90);

            %Labels
            zlabel('t (s)','interpreter','latex','fontsize',12);
            xlabel('x (um)','interpreter','latex','fontsize',12);
            ylabel('y (um)','interpreter','latex','fontsize',12);
        end

        function imageH = plotBackgroundFrame(obj, frame)
            %
            % IN:
            %   argin - either a frame index (scalar int) or an image to display on the background
            if obj.DrawBackgroundFrames %Only draw frames if DrawBackgroundFrames is set to true
               if isscalar(frame)
                    im = obj.frames(:,:,frame-(obj.frameBounds(1)-1)); %treat as index
                else
                    im = RPT.normalizeIndividualFrames(frame, obj.ImageColorRange);
                end
                BX = repmat(obj.ROIphysical(1:2),2,1);
                BY = repmat(obj.ROIphysical(3:4)',1,2);
                BZ = repmat(obj.tBounds(1),2,2); %draw bg image at Z=minTime
                imageH = surface(BX,BY,BZ,im,'FaceColor','texturemap','EdgeColor','none',...
                                 'FaceAlpha',obj.ImageAlpha,'CDataMapping','direct');
            else
                imageH = [];
            end
        end
        
        function trackHs = plotTracksFrame(obj, maxframe, minframe)
            % IN:
            %  maxframe - [optional] the maximum frame to show tracks for [default: last frame]
            %  minframe - [optional] the minimum frame to show tracks for [default: first frame]
            if nargin<3
                minframe = obj.frameBounds(1);
            end
            if nargin<2
                maxframe = obj.frameBounds(2);
            end
            hold('on');
            Ts = obj.tracks(cellfun(@(t) size(t,1)>1,obj.tracks));
            trackHs = cell(1,length(Ts));
            for i = length(Ts):-1:1 %go backwards to allocate trackHs array autmatically
                T=Ts{i};
                minIdx = find(T(:,end)>=minframe, 1, 'first');
                maxIdx = find(T(:,end)<=maxframe, 1, 'last');
                if isempty(minIdx) || isempty(maxIdx) || minIdx==maxIdx; continue; end %Track not visible
                trackHs{i} = BaseRPT.drawTrackSurface(T(minIdx:maxIdx,:), obj.trackCs{i}(minIdx:maxIdx));
            end
            hold('off');
            trackHs = [trackHs{:}];
        end

        function locHs = plotLocsFrame(obj, frameIdx)

            %Make a half-sphere of given size with alpha shading based on distance
            %This is the object we will use to represent localizations
            S = BaseRPT.makeLocalizationSurface();
            hold('on');
            locHs = cell(1,length(obj.tracks));
            for i = length(obj.tracks):-1:1
                T=obj.tracks{i};
                Lidx = find(T(:,end)==frameIdx, 1, 'first');
                if isempty(Lidx); continue; end %Track does not have a localization this frame
                locHs{i} = RPT.drawLocalizationSurface(T(Lidx,[2,3]), obj.LocalizationSigmaMultiple*T(Lidx,[7,8]), obj.tBounds(1)+eps, obj.trackCs{i}(Lidx), S);
            end
            hold('off');
            locHs = [locHs{:}];
        end
        
        function setFrames(obj, frames)
            if obj.ImageGlobalNormalize
                obj.frames = RPT.normalizeIndividualFramesGlobal(frames, obj.ImageColorRange);
            else
                obj.frames = RPT.normalizeIndividualFrames(frames, obj.ImageColorRange);
            end
        end
        
        function setSumImage(obj, sumim)
            obj.sumImage = RPT.normalizeIndividualFrames(sumim, obj.ImageColorRange);
        end
    end %public methods

    methods (Access=protected)        
        function loadRPT(obj, rpt, tracks)
            if nargin==2
                tracks = rpt.getTracks();
            elseif isnumeric(tracks) %tracks was given a list of IDs, and we'll use rpt.getTracks to get them
                tracks = rpt.getTracks(tracks);
            end
            RPT.checkTracks(tracks); %RPT knows how to check tracks for correctness.
            obj.tracks=tracks;
            obj.setFrames(rpt.getFrames());
            obj.setSumImage(rpt.sumImage);
            obj.ROIphysical = rpt.ROIPhysical;
            obj.frameT = rpt.data.frameT;
            obj.frameBounds = rpt.ROI(5:6);
            obj.imRef = rpt.data.getImRef(rpt.ROI);
        end
        
        function loadTrackCellArray(obj,tracks, frames, ROIphysical, frameBounds)
            if nargin>=3
                obj.ROIphysical=[0, size(frames,2), 0, size(frames,1)];
            else
                obj.ROIphysical = ROIphysical;
            end
            RPT.checkTracks(tracks); %RPT knows how to check tracks for correctness.
            obj.tracks = tracks;
            obj.setFrames(frames);
            obj.setSumImage(sumImage2D(frames));
            obj.frameT = (tracks{1}(2,1)-tracks{1}(1,1))/(tracks{1}(2,12)-tracks{1}(1,12)); % coputed as: (t2-t1)/(frame_idx2-frame_idx1)
            if nargin>=4
                obj.frameBounds = frameBounds;
            else
                obj.frameBounds = [min(cellfun(@(T) min(T(:,end)),obj.tracks)), max(cellfun(@(T) max(T(:,end)),obj.tracks))];
            end
            obj.imRef = imref2d(size(obj.frames(:,:,1)), obj.ROIphysical(1:2), obj.ROIphysical(3:4));
        end
    end %private methods   
end

