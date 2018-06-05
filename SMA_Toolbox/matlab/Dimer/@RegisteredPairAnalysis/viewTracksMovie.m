function  f = viewTracksMovie(obj, Tidxs, method, pBounds, tBounds)
    % This makes an interactive movie of the two tracks along with overlaid frames
    if nargin<3
        method = 'Wavelength';
    end
    if nargin<2 || isempty(Tidxs);
        Tidxs = 1:obj.nTracks;
        pBounds = obj.physicalBounds;
        tBounds = obj.timeBounds;
        fBounds = obj.frameBounds;
    else
        if nargin<4 || isempty(pBounds)
            pBounds = RPT.trackPhysicalBounds(obj.tracks(Tidxs));
        end
        if nargin<5 || isempty(tBounds)
            tBounds = RPT.trackTimeBounds(obj.tracks(Tidxs));
            fBounds = RPT.trackFrameBounds(obj.tracks(Tidxs));
        end
    end
    f=figure('Units','pixels','Position',[5,5,obj.FigureSize]);
    f.SizeChangedFcn = @sizeChanged_CB;
    whitebg(f);
    f.Color=[0,0,0];
    axH=axes();
    obj.RenderTracks = obj.tracks(Tidxs);
    obj.RenderTrackCs = obj.setTrackColors(Tidxs,method);
    obj.initializeAxes(axH, pBounds, tBounds);
    [nim1, nim2] = obj.normalizeFrameColors(obj.rpt_channel1.getFrames(), obj.rpt_channel2.getFrames());
    obj.RenderNormIms={nim1,nim2};
    obj.RenderImCoords = obj.makeImageSurfaceCoords(tBounds(1));
    lag =0;
    imHs = obj.plotBackgroundFrames(fBounds(1));
    trackHs = obj.plotTracksFrame(fBounds(1));
    locHs = obj.plotLocsFrame(fBounds(1));
    if obj.TitleOn
        title(sprintf('Frame:%i Time %.3f (s)',fBounds(1),(fBounds(1)-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
    end
    drawnow();

    uicontrol('Style','slider','Min',fBounds(1), 'Max',fBounds(2),'Value',fBounds(1),...
              'SliderStep',[1/(diff(fBounds)) 10/(diff(fBounds))],...
              'Units','normalized','Position',[0.05, 0, 0.8, 0.035], 'Callback',@slider_CB);
            
    function slider_CB(hObj,~)
        idx=round(hObj.Value);
        delete(imHs);
        delete(trackHs);
        delete(locHs);
        hold('on');
        imHs = obj.plotBackgroundFrames(idx);
        if lag>0
            trackHs = obj.plotTracksFrame(idx,idx-lag);
        else
            trackHs = obj.plotTracksFrame(idx);
        end
        locHs = obj.plotLocsFrame(idx);
        if obj.TitleOn
            title(sprintf('Frame:%i Time %.3f (s)',idx,(idx-1)*obj.frameT),'interpreter','latex','fontsize',obj.TitleFontSize);
        end
    end

    function sizeChanged_CB(~,~)
        pos = f.Position;
        pos(1:2) = 0;
        buf = [ 10 50 10 10];
        axH.ActivePositionProperty='OuterPosition';
        GUIBuilder.positionAxes(axH,pos,buf); 
    end
end


