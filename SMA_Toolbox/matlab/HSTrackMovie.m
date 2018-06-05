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

classdef HSTrackMovie < TrackMovie
    properties
        lambdaColormap;
        framesColormap;
        sumImageColormap;
    end

    methods
        function obj = HSTrackMovie(varargin)
            % [case1]: Load from hsrpt
            %   tm = TrackMovie(rpt, tracks)
            % hsrpt - HSRPT object or filename
            % tracks - [optional] an HSRPT format tracks cell-array to use.  
            %          [default=use all tracks in HSRPT object]
            % [case2]:
            %   tm = TrackMovie(tracks, frames, ROIphysical)
            % tracks - A RPT tracks cell-array formated set of tracks
            % frames - The full 4D [L Y X T] movie for this ROI  
            % ROIphysical - [optional] Physical region of space spaned by image [xmin xmax ymin ymax Lmin Lmax]
            %               Note: we use the matlab image coordinates where rows=L cols=y slices=x
            %               [default = [0 size(3) 0 size(2) 0 size(1)]
            if nargin>0
                obj.initialize(varargin{:});
            end
        end

        function initialize(obj,varargin)
            %Arguments the same as for the constructor
            if nargin<1
                error('HSTrackMovie:initialize','No arguments given');
            elseif isa(varargin{1},'HSRPT')
                obj.loadHSRPT(varargin{:}); %Load from HSRPT object
            elseif ischar(varargin{1})
                fname=varargin{1}; %Load from .hsrpt file
                [~,~,ext]=fileparts(fname);
                switch ext
                    case '.hsrpt'
                        obj.loadHSRPT(HSRPT(fname));
                    otherwise
                        error('HSTrackMovie:initialize','Unknown file type: %s',ext);
                end
%             elseif iscell(varargin{1})
%                 obj.loadTrackCellArray(varargin{:}); %Direct load from cell-array of tracks
            end
            %Common initializtion
            obj.tBounds = (obj.frameBounds-1)*obj.frameT;
        end

        function setTrackColorMethod(obj, methodName)
            switch methodName
                case 'Wavelength'
                    obj.trackCs = cellmap(@(i) obj.tracks{i}(:,4), 1:length(obj.tracks));
                    obj.ColorbarLabel = 'Wavelength (nm)';
                    obj.ColorbarOn = true;
                    obj.TrackColorMap = obj.lambdaColormap;
%                     obj.TrackColorMap = 
                    obj.TrackColorRange=size(obj.TrackColorMap,1);
                otherwise
                    setTrackColorMethod@TrackMovie(obj, methodName);
                    return
            end
            obj.TrackColorMethod = methodName;
            [obj.trackCs, obj.trackMin, obj.trackMax] = RPT.normalizeTrackColors(obj.trackCs, obj.TrackColorRange, obj.ImageColorRange);
        end

        function initializeColorMap(obj)
            if ishandle(obj.TrackColorMap)
                colormap([obj.framesColormap; obj.TrackColorMap(obj.TrackColorRange)]);
            else
                colormap([obj.framesColormap; obj.TrackColorMap]);
            end
        end

        function [xErr,yErr]=getUncertainty(obj,trackIdx,locIdx)
             xErr=obj.LocalizationSigmaMultiple*obj.tracks{trackIdx}(locIdx,10);
             yErr=obj.LocalizationSigmaMultiple*obj.tracks{trackIdx}(locIdx,11);
        end


        function h=viewSequence(obj, frameIdx)
            % frame - [optional] initial frame
            if nargin==1
                frameIdx = obj.frameBounds(1);
            end
            frameLag=obj.FrameLag;
            obj.figureH=figure();
            obj.figureH.Units='Pixels';
            obj.figureH.Position=[10 10 obj.FigureSize];
            obj.axH=obj.initializeAxes();
            xsz=obj.ROIphysical(2)-obj.ROIphysical(1);
            ysz=obj.ROIphysical(4)-obj.ROIphysical(3);
            obj.axH.Units='Pixels';
            GUIBuilder.positionImageAxes(obj.axH,[xsz,ysz],[0, 0, obj.FigureSize],[20, 20, 40, 20]); 

            imH=obj.plotBackgroundFrame(frameIdx);
            trackHs=obj.plotTracksFrame(frameIdx);
            locHs=obj.plotLocsFrame(frameIdx);
            if obj.TitleOn
                title(sprintf('Frame:%i Time %.3f (s)',frameIdx,(frameIdx-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
            end
            drawnow();
            sliderH = uicontrol('Style','slider','Min',obj.frameBounds(1), 'Max',obj.frameBounds(2),'Value',frameIdx,...
                                'SliderStep',[1/(1+diff(obj.frameBounds)) 10/(diff(obj.frameBounds))],...
                                'Units','normalized','Position',[0.1, 0, 0.8, 0.055], 'Callback',@slider_CB);
            
            
            function slider_CB(hObj,~)
                idx=round(hObj.Value);
                delete(imH);
                delete(trackHs);
                delete(locHs);
                imH=obj.plotBackgroundFrame(idx);
                if frameLag<=0
                    trackHs=obj.plotTracksFrame(idx);
                else
                    trackHs=obj.plotTracksFrame(idx,idx-frameLag);
                end
                locHs=obj.plotLocsFrame(idx);
                if obj.TitleOn
                    title(sprintf('Frame:%i Time %.3f (s)',idx,(idx-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
                end
            end
            
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

            %Set aspect ratio of axes
            axis([obj.ROIphysical obj.tBounds]);
            asp=pbaspect();
            asp(1) = (obj.ROIphysical(2)-obj.ROIphysical(1))/(obj.ROIphysical(4)-obj.ROIphysical(3));
            asp(3) = mean([asp(1),asp(2)]);
            pbaspect(asp);
            %Fix grids and ticks
            view(0,90);

            %Labels
            zlabel('t (s)','interpreter','latex','fontsize',obj.AxesFontSize);
            xlabel('x (um)','interpreter','latex','fontsize',obj.AxesFontSize);
            ylabel('y (um)','interpreter','latex','fontsize',obj.AxesFontSize);            
        end

        function imageH = plotBackgroundFrame(obj, frame)
            %
            % IN:
            %   argin - either a frame index (scalar int) or an image to display on the background
            if isscalar(frame)
                im = obj.frames(:,:,frame-obj.frameBounds(1)+1); %treat as index
            else
                im = obj.normalizeImage(frame);
            end
            BX = repmat(obj.ROIphysical(1:2),2,1);
            BY = repmat(obj.ROIphysical(3:4)',1,2);
            BZ = repmat(obj.tBounds(1),2,2); %draw bg image at Z=minTime
            imageH = surface(BX,BY,BZ,im,'FaceColor','texturemap','EdgeColor','none',...
                             'FaceAlpha',obj.ImageAlpha,'CDataMapping','direct');            
        end
        
        function trackHs = plotTracksFrame(obj, maxframe, minframe)
            % IN:
            %  maxframe - [optional] the maximum frame to show tracks for [default: last frame]
            %  minframe - [optional] the minimum frame to show tracks for [default: first frame]
            if nargin<3
                minframe = obj.frameBounds(1);
            else
                minframe = max(minframe, obj.frameBounds(1));
            end
            if nargin<2
                maxframe = obj.frameBounds(2);
            end
            hold('on');
            Ts = obj.tracks(cellfun(@(t) size(t,1)>1,obj.tracks));
            trackHs = cell(1,length(Ts));
            nHs=0;
            for i = length(Ts):-1:1 %go backwards to allocate trackHs array autmatically
                T=Ts{i};
                minIdx = find(T(:,end)>=minframe, 1, 'first');
                maxIdx = find(T(:,end)<=maxframe, 1, 'last');
                if isempty(minIdx) || isempty(maxIdx); continue; end %Track not visible
                if minIdx==maxIdx; continue; end %Only a single localization.  Too short to draw
                ts = T(minIdx:maxIdx,1)';
                xs = T(minIdx:maxIdx,2)';
                ys = T(minIdx:maxIdx,3)';
                C = obj.trackCs{i}(minIdx:maxIdx)';
                nHs = nHs+1;
                trackHs{nHs} = surface([xs;xs],[ys;ys],[ts;ts],[C;C],'EdgeColor','interp',...
                                     'LineWidth',obj.TrackLineWidth,...
                                     'Marker',obj.TrackVertexMarker);
                %trackHs{nHs}.AlphaDataMapping='Direct';
                trackHs{nHs}.CDataMapping='Direct';
            end
            hold('off');
            trackHs = [trackHs{1:nHs}];
        end

%        function f=renderVolumeWavelength(obj)
% 
%             for i=frameIdxs
%                 [im,roi]=obj.getPlottableFrame(frame_idx,roi_in);
%             [xs,ys,Ls]=obj.pixelGridValues(roi);
%             if frame_idx>0
%                 fs=obj.getFrames(roi_in);
%                 max_val = 0.75*max(fs(:));
%                 im(:)=min(im(:),max_val);
%             end
%             HSData.volumeSliceView(xs,ys,Ls,im, jet,'wavelength',ax);
%             name=sprintf('Frame - %i',roi(7)-1+frame_idx);
%             set(f,'Name',name);
%             cbh=colorbar();
%             set(get(cbh,'Label'),'String','Wavelength $\lambda$ (nm)');
%             obj.label3DFigure(f, cbh);
%             view_angles=[171, 2];
%             view(view_angles);
%             set(gcf(),'Position',[10 10 600 450]);
%             sliderStep=[1/(obj.nFrames-1) 10/(obj.nFrames-1)];
%             function slider_CB(s,~)
%                 idx=round(s.Value);
%                 [im,roi]=obj.getPlottableFrame(idx,roi);
%                 [view_a1, view_a2]=view();
%                 cla(ax);
%                 im(:)=max(0,min(im(:),max_val)-0.02*max_val);
%                 HSData.volumeSliceView(xs,ys,Ls,im, jet,'wavelength',ax);
%                 name=sprintf('Frame - %i',roi(7)-1+idx);
%                 set(f,'Name',name);
%                 set(get(cbh,'Label'),'String','Wavelength $\lambda$ (nm)');
%                 obj.label3DFigure(f, cbh);
%                 view([view_a1 view_a2]);
%             end
%             if frame_idx>0
%                 total_frame=obj.getFrames(roi_in);
% %                 mcount=max(total_frame(:));
% %                 caxis([0 mcount*.75]);
%                 uicontrol('Style','slider','Min',obj.globalTBounds(1),'Max',obj.globalTBounds(2),'Value',frame_idx,...
%                           'Position',[10,435,580,20],'SliderStep',sliderStep, 'Callback',@slider_CB);
%             end
%         end

        function setFrames(obj, frames)
            if size(frames,ndims(frames))==3 
                RGB=frames;
            else
                RGB = HSData.makeRGB(frames, obj.lambdaColormap);
            end
            sz=size(RGB);
            RGB = reshape(RGB, [], sz(3), 3); %resize to 2Dx3color image and scale all at once
            [obj.frames, obj.framesColormap] = rgb2ind(RGB, obj.ImageColorRange);
            obj.ImageColorRange=size(obj.framesColormap,1);
            obj.frames = reshape(obj.frames, sz(1:3));
        end
        function setSumImage(obj, sumim)
            assert(size(sumim,ndims(sumim))==3); %check we are an rgb
            [obj.sumImage, obj.sumImageColormap ] = rgb2ind(sumim, obj.ImageColorRange);
        end
    end %Public methods

    methods (Access=protected)
        function loadHSRPT(obj, hsrpt, tracks)
            if nargin==2 || isempty(tracks)
                tracks = hsrpt.getTracks();
            elseif isnumeric(tracks) %tracks was given a list of IDs, and we'll use rpt.getTracks to get them
                tracks = hsrpt.getTracks(tracks);
            end
            HSRPT.checkTracks(tracks); %RPT knows how to check tracks for correctness.
            obj.ImageColorRange = 1024;
            obj.TrackColorMethod = 'Wavelength';
            obj.tracks = tracks;
            ROI = hsrpt.ROI;
            obj.lambdaColormap = flip(hsrpt.data.colorMap(ROI(5):ROI(6),:));
            obj.setFrames(hsrpt.getFramesRGB());
            obj.setSumImage(hsrpt.sumImage);
            obj.imRef = hsrpt.data.getImRef();
%             roi_phys = hsrpt.ROIPhysical;
            obj.ROIphysical = hsrpt.ROIPhysical(1:4);
            obj.frameT = hsrpt.data.frameT;
            obj.frameBounds = ROI(7:8);

            obj.FigureSize=[780, 1024];
        end
    end % protected methods

    methods (Static=true)
    end %static methods
end
