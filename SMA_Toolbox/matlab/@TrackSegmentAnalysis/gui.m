% TrackSegmentAnalysis.gui()

function guiFig = gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end

    %% Constants and Global Initialization
    gui@GUIBuilder(obj); %Call GUIBuilder initialization
    gui_name='Track Segment Analysis';
    
    lastPath=obj.workingDir; %save this for file open dialogs after a reset
    LocsTableColFormat=struct(); %Store the relevent column parameters for the Localizations Table
    SegsTableColFormat=struct(); %Store the relevent column parameters for the Segments Table
    TracksTableColFormat=struct();%Store the relevent column parameters for the Tracks Table
    uH = GUIBuilder.default_unitHeight; % unit height for elements
    boarder = GUIBuilder.default_boarder;%Boarder width around the outside of the gui
    sp = GUIBuilder.default_spacing; %spacing between elements.
    but_sz = GUIBuilder.default_buttonSize; %Button size
    savedPlotInfo.trackId=[];
    fig_sz=[1120 800]; %figure size

    fullw_sz=[fig_sz(1)-2*boarder-2*sp uH]; %size of a full-width component that stays within the boarders.
    halfw=512;
    plot_height = 500;
    halfw_sz=[halfw-2*sp uH];
    if isunix
        ExportExt = '.csv';
        ExportFormats = obj.CSVFormatsLinux;
    else
        ExportExt = '.xls';
        ExportFormats = obj.CSVFormatsWin;
    end
    %% Create figure and controls
    % Make figure
    guiFig = figure('Units','pixels','Position',[10 0 fig_sz],'Resize','off',...
                    'MenuBar','none','ToolBar','none','NumberTitle','off',...
                    'Name',gui_name,'Visible','on',...
                    'CloseRequestFcn',@close_CB);
    obj.guiFig=guiFig;

    

    tabs_pos=[boarder boarder fullw_sz(1) fig_sz(2)-plot_height];
    halfPan_pos=[boarder boarder+halfw halfw_sz];
    handles.tabG = uitabgroup('Parent',guiFig,'Units','Pixels','Position',tabs_pos,'TabLocation','top');
    handles.tabs.Locs=uitab('Parent',handles.tabG ,'Title','Localizations');
    handles.tabs.file=uitab('Parent',handles.tabG ,'Title','File');
    handles.panels.Segments=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Segments');
    handles.panels.Tracks=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Tracks');
    
    createMenus()
    populateFilePanel(handles.tabs.file);    
    populateLocsPanel(handles.tabs.Locs);
    populateSegmentsPanel(handles.panels.Segments);
    populateTracksPanel(handles.panels.Tracks);
    
    handles.axes= axes('Units', 'Pixels','Visible','off');
    alignPanels();
    positionAxes(handles.axes);
    setTracksTableData();
    setSegsTableData();
    setLocsTableData();
    updateControls();
    pause(0.05); %See if this helps with the redraw maybe JavaEDT issues?
    drawnow();

    function populateFilePanel(panH)
        % File Panel
        % Bottom row of buttons
        bot_row.pos=[sp, sp, halfw_sz];
        bot_row.names={'Load','Save','SaveAs'};
        bot_row.CBs={@load_CB,@save_CB,@saveAs_CB};
        GUIBuilder.buttonRow(panH, bot_row.pos, but_sz, bot_row.names, bot_row.CBs);
        %Quit is out to the side so we don't include it.
        handles.quitButton = uicontrol('Parent',panH,'Style','pushbutton','String','Quit',...
            'Position',[fullw_sz(1)-but_sz(1) sp but_sz],...
            'Callback',@close_CB);
        col_names={'Associated Files','Full Path'};
        table_pos = [sp, 2*sp+uH, fullw_sz(1)-2*sp 120];
        col_width = {200, table_pos(3)-200-2*boarder-2*sp};
        handles.fileTable = uitable(panH, 'ColumnName', col_names, 'SelectionHighlight','off',...
                                    'Position', table_pos, 'ColumnWidth',col_width,'RowName',[]);
        setFileTableData();
    end

   

    function populateTracksPanel(panH)
        % File Panel
        % Bottom row of buttons
        table_pos=[sp sp halfw_sz(1) plot_height/2-6*boarder];
        col_names={'Track<br />ID','Num<br />Segs','Segment<br />Class','Time(s)','Num<br />Locs',...
                   'Total<br />dist(um)','Net<br />Dist(um)','Max<br />Dist(um)',...
                   'Confinement<br />Ratio','D_MLE<br />(um^2/s)',...
                   'Mean<br />Spd(um/s)','Net<br />Spd(um/s)','Linearity<br />Ratio'};
        col_names=cellmap(@(t) sprintf('<html><div align="center">%s</div></html>',t) ,col_names);
        col_formats={'', '', 'char', '', '', '', '', '','','','','',''};
        col_editable=false(1,length(col_names));
        col_editable(3)=true;
        data=table2cell(obj.trackStatsTable());
        col_widths={40, 40 , 'auto','auto', 40,'auto','auto','auto',85,'auto','auto','auto','auto'};        
        TracksTableColFormat.names = col_names;
        TracksTableColFormat.formats = col_formats;
        TracksTableColFormat.widths = col_widths;
        TracksTableColFormat.editable = col_editable;
        [handles.TracksTable, handles.TracksTablePanel] = GUIBuilder.makeTreeTable(panH,...
                                            TracksTableColFormat, data, table_pos, 'Groupable', false);
        handles.TracksTable.MousePressedCallback = {@GUIBuilder.javaTableMouse_CB,handles.contextmenus.tracks};
    end
    

    function populateSegmentsPanel(panH)
        % File Panel
        % Bottom row of buttons
        table_pos = [sp sp halfw_sz(1) plot_height/2-6*boarder];
        col_names={'Segment<br />Id','Segment<br />Class','Time(s)','Num<br />Locs',...
                   'Total<br />dist(um)','Net<br />Dist(um)','Max<br />Dist(um)',...
                   'Confinement<br />Ratio','D_MLE<br />(um^2/s)', ...
                   'Mean<br />Spd(um/s)','Net<br />Spd(um/s)','Linearity<br />Ratio'};
        col_names=cellmap(@(t) sprintf('<html><div align="center">%s</div></html>',t) ,col_names);
        col_names=['TrackID', col_names]; %First column can't be html
        col_formats={'', '', 'char', '', '', '', '', '','','','','',''};
        data=[];
        col_widths={100, 30, 30,30,30,30,30,30,30,30,30,30};
        col_editable=false(1,length(col_widths));
        col_editable(3)=true;
        SegsTableColFormat.names=col_names;
        SegsTableColFormat.formats=col_formats;
        SegsTableColFormat.widths=col_widths;
        SegsTableColFormat.editable = col_editable;
        [handles.SegsTable, handles.SegsTablePanel] = GUIBuilder.makeTreeTable(panH,...
                                        SegsTableColFormat, data, table_pos, 'Groupable', true);
        handles.SegsTable.MousePressedCallback = {@GUIBuilder.javaTableMouse_CB,handles.contextmenus.segs};
    end

    function populateLocsPanel(panH)
        % File Panel
        % Bottom row of buttons
        table_pos = [sp sp fullw_sz(1)-boarder fig_sz(2)-plot_height-6*boarder];
        col_names={'Localization<br />ID','Frame<br />Idx',...
                   'Time<br />(s)','x<br />(um)','y<br />(um)','I<br />(photons)','Background<br />(photons/px)','Sigma',...
                   'Displacement(um)','Total<br />CumDist(um)','Net<br />CumDist(um)','Inst. Speed<br />(um/s)',...
                   'Direction<br />Deg','Angle<br />Deg'};
        col_names=cellmap(@(t) sprintf('<html><div align="center">%s</div></html>',t) ,col_names);
        col_names=['TrackId', 'SegId', col_names]; %First 2-column can't be html
        ncols = length(col_names);
        data = [];
        [col_widths{1:ncols}] = deal(80);
        [col_formats{1:ncols}] = deal('');
        col_editable=false(1,ncols-1);
        LocsTableColFormat.names=col_names;
        LocsTableColFormat.formats=col_formats;
        LocsTableColFormat.widths=col_widths;
        LocsTableColFormat.editable = col_editable;
        [handles.LocsTable, handles.LocsTablePanel] = GUIBuilder.makeTreeTable(panH,...
                                        LocsTableColFormat, data, table_pos, 'Groupable', true);
        set(handles.LocsTable,'MousePressedCallback',{@GUIBuilder.javaTableMouse_CB,handles.contextmenus.locs});
    end

    function createMenus()
        %File Menu

        labels = {'Load ...','Save','Save As ...','Reset Object',[],'Quit'};
        CBs = {@load_CB, @save_CB, @saveAs_CB, @resetObject_CB,[],@close_CB};
        GUIBuilder.makeFigureMenu(guiFig,'File',labels,CBs);

        labels = {'Plot All Segments Sequence','Plot All Tracks Sequence','Plot All Tracks Speed','Plot All Temporal',[],...
                  'View All Segments Sequence...', 'View All Tracks Sequence...','View All Tracks Speed...','View All Temporal...',[],...
                  'View All Tracks SquaredDisplacementCDF...'};
        CBs = {@(~,~) plotMainAxes(@obj.plot3DSegmentSequence),@(~,~) plotMainAxes(@obj.plot3DTrackSequence), @(~,~) plotMainAxes(@obj.plot3DTrackSpeed),@(~,~) plotMainAxes(@obj.plot3DTrackTemporal),[],...
               @(~,~) obj.view3DSegmentSequence, @(~,~) obj.view3DTrackSequence(), @(~,~) obj.view3DTrackSpeed(), @(~,~) obj.view3DTrackTemporal(),[],...
               @(~,~) obj.viewSquaredDisplacementCDFOverall()};
        GUIBuilder.makeFigureMenu(guiFig,'View',labels,CBs); 
        
        labels = {'Export All...','Export All Track Stats...', 'Export All Segment Stats...', 'Export All Localization Stats...'};
        CBs = {@exportAllStats_CB, @exportAllTrackStats_CB, @exportAllSegStats_CB, @exportAllLocStats_CB};
        GUIBuilder.makeFigureMenu(guiFig,'Export',labels,CBs);                

        labels = {'List Segments','List Localizations',[],...
                  'Plot Selected Tracks Sequence','Plot Selected Tracks Speed','Plot Selected Tracks Temporal',[],...
                  'View Selected Tracks Sequence...','View Selected Tracks Speed...','View Selected Tracks Temporal...',[],...
                  'View Selected Tracks Sq.Disp. Overall CDF...','View Selected Tracks Sq.Disp. Individual CDF...','View Selected Tracks MSD Temporal...',[],...
                  'ExportSelectedTrackStats...', 'Export Selected Localization Stats...', [],...
                  'Join Tracks', 'Delete Tracks'};
        CBs = {@listSegments_CB,@listTrackLocalizations_CB,[],...
               @plotTrackSequence_CB, @plotTrackSpeed_CB, @plotTrackTemporal_CB, [],...
               @viewTrackSequence_CB, @viewTrackSpeed_CB, @viewTrackTemporal_CB, [],...
               @viewTrackSqDispCDFOverall_CB,@viewTrackSqDispCDFIndividual_CB,@viewTrackMSDTemporal_CB, [],...
               @exportTrackStats_CB,@exportTrackLocStats_CB,[],@joinTracks_CB, @deleteTracks_CB};
        handles.contextmenus.tracks = GUIBuilder.makeJavaContextMenu(labels, CBs);

        labels = {'List Localizations',[],...
                  'Plot Selected Segments Sequence','Plot Selected Segments Speed','Plot Selected Segments Temporal',[],...
                  'View Selected Segments Sequence...','View Selected Segments Speed...','View Selected Segments Temporal...',[],...
                  'View Selected Segments Sq.Disp. Overall CDF...','View Selected Segments SqDisp Individual CDF...','View Selected Segments MSD Temporal...',[],...
                  'Export Selected Segments Stats...','Export Selected Localization Stats...',[],...
                  'JoinSegments','Delete Segments'};
        CBs = {@listSegLocalizations_CB,[],...
               @plotSegSequence_CB, @plotSegSpeed_CB, @plotSegTemporal_CB, [],...
               @viewSegSequence_CB, @viewSegSpeed_CB, @viewSegTemporal_CB, [],...
               @viewSegmentSqDispCDFOverall_CB,@viewSegmentSqDispCDFIndividual_CB,@viewSegmentMSDTemporal_CB, [],...
               @exportSegStats_CB,@exportSegLocStats_CB,[],@joinSegments_CB, @deleteSegments_CB};
        handles.contextmenus.segs = GUIBuilder.makeJavaContextMenu(labels, CBs);

        labels = {'Export Localization Stats...',[],'Divide Segment', 'Divide Track',[], 'Delete Localization'};
        CBs = {@exportLocStats_CB,[],@divideSegment_CB,@divideTrack_CB,[],@deleteLocalization_CB};
        handles.contextmenus.locs = GUIBuilder.makeJavaContextMenu(labels, CBs);
        
        
        
        %DataCursorMode UIControlModification
        labels = {'List Segment', 'List Localizations',[],'View Segment...',[],...
                  'Divide Segment', 'Divide Track',[],...
                  'View Sq.Disp. CDF...','View MSD Temporal...',[],...
                  'Delete Localization','Delete Track'};
        CBs = {@listSegmentCursor_CB, @listLocalizationsCursor_CB,[],...
               @viewSegmentCursor_CB,[],...
               @divideSegmentCursor_CB, @divideTrackCursor_CB, [],...
               @viewSqDispCDFCursor_CB,@viewMSDTemporalCursor_CB, [],...
               @deleteLocaliztionCursor_CB, @deleteTrackCursor_CB};
        handles.contextmenus.datacursor = GUIBuilder.makeContextMenu(labels, CBs);
        dcm = datacursormode(gcf());
        dcm.UIContextMenu = handles.contextmenus.datacursor;
    end

    

    function panels_pos=alignPanels()
        pans={handles.tabG, handles.panels.Segments, handles.panels.Tracks};
        for n=1:length(pans)
            GUIBuilder.autoSizePanel(pans{n});
        end
        panels_pos=GUIBuilder.align(pans,'Left','Fixed',sp);
    end

    function positionAxes(axesH)
        pos = round(handles.panels.Segments.Position);
        left_pos = pos(1)+pos(3)+2*sp;
        axes_pos = [left_pos, pos(2)+uH+sp, fig_sz(1)-left_pos-boarder, fig_sz(2)-pos(2)-boarder-uH-sp];
        labelmargin = [20 20 20 0];
        GUIBuilder.positionAxes(axesH, axes_pos, labelmargin);
        axesH.Clipping='on';
        
        %These are the buttons to change the plot mode that emulate the toolbar
        c_pos = [left_pos, pos(2), uH, uH];
        but_pos = [left_pos, pos(2), 70, uH];
        H{1} = uicontrol('Parent',guiFig,'Style','togglebutton', 'CData',GUIBuilder.readMatlabIcon('tool_hand.png'),...
                         'Position',c_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('pan',a,b));
        H{2} = uicontrol('Parent',guiFig,'Style','togglebutton', 'CData',GUIBuilder.readMatlabIcon('tool_rotate_3d.png'),...
                         'Position',c_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('rotate3d',a,b));
        H{3} = uicontrol('Parent',guiFig,'Style','togglebutton', 'CData',GUIBuilder.readMatlabIcon('tool_zoom_in.png'),...
                         'Position',c_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('zoom',a,b));
        H{4} = uicontrol('Parent',guiFig,'Style','togglebutton', 'CData',GUIBuilder.readMatlabIcon('tool_data_cursor.png'),...
                         'Position',c_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('datacursormode',a,b));
        H{5} = uicontrol('Parent',guiFig,'Style','pushbutton', 'String','SaveView',...\
                         'Position',but_pos,'HorizontalAlignment','center','Callback',@(~,~) saveAllPlotView());
        H{6} = uicontrol('Parent',guiFig,'Style','pushbutton', 'String','Replot',...
                         'Position',but_pos,'HorizontalAlignment','center','Callback',@restorePlotView_CB);
        GUIBuilder.align(H,'Fixed',sp,'Bottom');
        handles.axes_buttons.pan = H{1};
        handles.axes_buttons.rotate3d = H{2};
        handles.axes_buttons.zoom = H{3};
        handles.axes_buttons.datacursormode = H{4};
        plotMainAxes(@obj.plot3DSegmentSequence,[]);
    end

    function axesMode_CB(mode,~,~)
        modes = {'rotate3d','pan','zoom','datacursormode'};
        settings = {'off','on'};
        for i = 1:length(modes)
            m = modes{i};
            func = str2func(m);
            if strcmp(mode,m)
                state = handles.axes_buttons.(m).Value;
                func(settings{state+1});
            else
                handles.axes_buttons.(m).Value = 0;
                func(settings{1});
            end
        end
        
    end
    %%Table data maniputlation
    
    function setFileTableData()
        data={'RPT (.rpt)', obj.getFilePath('RPT');...
              'SPT (.spt)', obj.getFilePath('SPT');...
              'SPData (.spdata)', obj.getFilePath('data');...
              'CSV Export', obj.getFilePath('CSV');...
              'TrackSegmentAnalysis (.tsa) ', obj.saveFilePath};
        handles.fileTable.Data = data;
    end

    function setTracksTableData()       
        data = table2cell(obj.trackStatsTable());
        GUIBuilder.setTreeTableData(handles.TracksTable, TracksTableColFormat, data);
    end

    function setSegsTableData(varargin)
        data = table2cell(obj.segmentStatsTable(varargin{:}));
        GUIBuilder.setTreeTableData(handles.SegsTable, SegsTableColFormat, data, 1);
    end

    function setLocsTableData(varargin)
        data = table2cell(obj.trackLocalizationsTable(varargin{:}));
        GUIBuilder.setTreeTableData(handles.LocsTable, LocsTableColFormat, data, [1,2]);
    end

    function trackIds=getSelectedTracks()
        %Get the selected trackIds
        obj.assertInitialized()
        tabH=handles.TracksTable;
        trackIds = arrayfun(@(r) tabH.getValueAt(r,0), tabH.getSelectedRows());
    end

    function [trackIds,segIds]=getSelectedSegs()
        %Get the selected trackIds
        obj.assertInitialized()
        tabH=handles.SegsTable;
        treeModelH=tabH.getModel().getActualModel();
        tableModelH=treeModelH.getActualModel(); %This model does not have the tree columns collapsed
        actualRows = arrayfun(@(r) treeModelH.getActualRowAt(r), tabH.getSelectedRows());
        actualRows = actualRows(actualRows>=0);
        trackIds = arrayfun(@(r) tableModelH.getValueAt(r,0),actualRows);
        segIds = arrayfun(@(r) tableModelH.getValueAt(r,1), actualRows);
    end

    function [trackIds,segIds,locIds]=getSelectedLocs()
        %Get the selected trackIds
        obj.assertInitialized()
        tabH=handles.LocsTable;
        treeModelH=tabH.getModel().getActualModel();
        tableModelH=treeModelH.getActualModel(); %This model does not have the tree columns collapsed
        actualRows = arrayfun(@(r) treeModelH.getActualRowAt(r), tabH.getSelectedRows());
        actualRows = actualRows(actualRows>=0);
        trackIds = arrayfun(@(r) tableModelH.getValueAt(r,0),actualRows);
        segIds = arrayfun(@(r) tableModelH.getValueAt(r,1), actualRows);
        locIds = arrayfun(@(r) tableModelH.getValueAt(r,2), actualRows);
    end

    %% Updating
    function updateControls()
        axes(handles.axes);        
    end

    function reloadAllControls()
        setFileTableData();
        setTracksTableData();
        setSegsTableData([]);
        setLocsTableData([],[]);
    end

    %% File loading and saving callbacks
    function load_CB(~, ~)
        if isempty(lastPath)
            lastPath=obj.workingDir;
        end
        [filename, pathname]=uigetfile(obj.LoadableDataFormats,'Select data source', lastPath);
        if ~filename; return; end
        obj.load(fullfile(pathname,filename))
        setTracksTableData();
        positionAxes(handles.axes);
        reloadAllControls();
    end

    function save_CB(~,~)
        obj.save();
        updateControls();
    end

    function saveAs_CB(~,~)       
        obj.saveas();
        updateControls();
    end

    function resetObject_CB(~,~)
        if obj.initialized
            lastPath=obj.workingDir; %Only change the lastPath if we actually have a valid workingDir
        end
        obj.resetObject();
        cla(handles.axes);
        setTracksTableData();
        positionAxes(handles.axes);
        reloadAllControls();
    end

    function close_CB(~,~)
        try
            set(handles.TracksTable,'MousePressedCallback',[]);
            set(handles.SegsTable,'MousePressedCallback',[]);
            set(handles.LocsTable,'MousePressedCallback',[]);
            GUIBuilder.clearJavaContextMenu(handles.contextmenus.tracks);
            GUIBuilder.clearJavaContextMenu(handles.contextmenus.segs);
            GUIBuilder.clearJavaContextMenu(handles.contextmenus.locs);
        catch err
            disp(getReport(err));
        end
        obj.closeGUI();
    end

    function deleteTracks_CB(~,~)
        obj.deleteTrack(getSelectedTracks());
        reloadAllControls();
        plotMainAxes(@obj.plot3DTrackSequence)
    end

    %% Viewing and Plotting
    function plotTrackSequence_CB(~,~)
        plotMainAxes(@obj.plot3DTrackSequence, getSelectedTracks());
    end

    function plotTrackSpeed_CB(~,~)
        plotMainAxes(@obj.plot3DTrackSpeed, getSelectedTracks());
    end

    function plotTrackTemporal_CB(~,~)
        plotMainAxes(@obj.plot3DTrackTemporal, getSelectedTracks());
    end

    function viewTrackSequence_CB(~,~)
        obj.view3DTrackSequence(getSelectedTracks());
    end

    function viewTrackSpeed_CB(~,~)
        obj.view3DTrackSpeed(getSelectedTracks());
    end

    function viewTrackTemporal_CB(~,~)
        obj.view3DTrackTemporal(getSelectedTracks());
    end
    
    function plotSegSequence_CB(~,~)
        [tid,sid]=getSelectedSegs();
        plotMainAxes(@obj.plot3DSegmentSequence, tid, sid);
    end

    function plotSegSpeed_CB(~,~)
        [tid,sid]=getSelectedSegs();
        plotMainAxes(@obj.plot3DSegmentSpeed, tid, sid);
    end

    function plotSegTemporal_CB(~,~)
        [tid,sid]=getSelectedSegs();
        plotMainAxes(@obj.plot3DSegmentTemporal, tid, sid);
    end

    function viewSegSequence_CB(~,~)
        [tid,sid]=getSelectedSegs();
        obj.view3DSegmentSequence(tid, sid);
    end

    function viewSegSpeed_CB(~,~)
        [tid,sid]=getSelectedSegs();
        obj.view3DSegmentSpeed(tid, sid);
    end

    function viewSegTemporal_CB(~,~)
        [tid,sid]=getSelectedSegs();
        obj.view3DSegmentTemporal(tid, sid);
    end

    function viewTrackSqDispCDFOverall_CB(~,~)
        obj.viewSquaredDisplacementCDFOverall(getSelectedTracks());
    end
    function viewTrackSqDispCDFIndividual_CB(~,~)
        obj.viewSquaredDisplacementCDFIndividual(getSelectedTracks());
    end
    function viewTrackMSDTemporal_CB(~,~)
        obj.viewMSDTemporal(getSelectedTracks());
    end

    function viewSegmentSqDispCDFOverall_CB(~,~)
        [tid,sid] = getSelectedSegs();
        obj.viewSquaredDisplacementCDFOverall(tid,sid);
    end
    function viewSegmentSqDispCDFIndividual_CB(~,~)
        [tid,sid] = getSelectedSegs();
        obj.viewSquaredDisplacementCDFIndividual(tid,sid);
    end
    function viewSegmentMSDTemporal_CB(~,~)
        [tid,sid]=getSelectedSegs();
        obj.viewMSDTemporal(tid,sid);
    end
    
    function viewSqDispCDFCursor_CB(~,~)
        [tid,sid] = getDataCursorIndexPosition();
        obj.viewSquaredDisplacementCDFIndividual(tid, sid);
    end
    function viewMSDTemporalCursor_CB(~,~)
        [tid,sid] = getDataCursorIndexPosition();
        obj.viewMSDTemporal(tid, sid);
    end
    %% DataCursor Mode UIContextMenu callbacks
    function [trackId,segId,locId] = getDataCursorIndexPosition()
        dcm=datacursormode(guiFig);
        info=dcm.getCursorInfo();
        [trackId,segId,locId] = deal(3);
        if ~isempty(info)
            trackH = info.Target;
            pos = info.Position;
            data = trackH.UserData;
            if ~isempty(data)
                stats=data{1};
                Ts=data{2};
                trackId = stats.TrackID;
                segId = stats.SegmentID;        
                locId = binarysearch(Ts(:,1),pos(3));
            end
        end
    end
    
    function listSegmentCursor_CB(~,~)
        [trackId, segId] = getDataCursorIndexPosition();
        setSegsTableData(trackId,segId);       
    end
    
    function listLocalizationsCursor_CB(~,~)
        [trackId, segId] = getDataCursorIndexPosition();
        setLocsTableData(trackId,segId);
    end
    function viewSegmentCursor_CB(~,~)
        [tid,sid] = getDataCursorIndexPosition();
        obj.view3DSegmentSequence(tid, sid);
    end
    function divideSegmentCursor_CB(~,~)
        [trackId, segId, locId] = getDataCursorIndexPosition();
        divideSegmentUpdate(trackId,segId,locId);
    end
  
    function divideTrackCursor_CB  (~,~)
        [trackId, segId, locId] = getDataCursorIndexPosition();
        divideTrackUpdate(trackId,segId,locId);
    end

    function deleteTrackCursor_CB(~,~)
        [trackId, ~, ~] = getDataCursorIndexPosition();
        obj.deleteTrack(trackId);
        reloadAllControls();
        saveAllPlotView();
        t=timer('StartDelay', 0.01, 'ExecutionMode','singleShot','TimerFcn',@restorePlotView_CB);
        start(t);
    end

    function divideSegment_CB(~,~)
        [trackId, segId, locId] = getSelectedLocs();
        divideSegmentUpdate(trackId,segId,locId);
    end

    function divideTrack_CB(~,~)
        [trackId, segId, locId] = getSelectedLocs();
        divideTrackUpdate(trackId,segId,locId);
    end
 
    function divideSegmentUpdate(trackId,segId,locId)
        if isempty(trackId)
            return
        end
        obj.divideSegment(trackId,segId,locId);
        setTracksTableData();
        setSegsTableData(trackId,[]);
        setLocsTableData(trackId,[]);
        saveAllPlotView()
        savedPlotInfo.trackId=trackId;
        t=timer('StartDelay', 0.01, 'ExecutionMode','singleShot','TimerFcn',@restorePlotView_CB);
        start(t);
    end
    
    function divideTrackUpdate(trackId,segId,locId)
        if isempty(trackId)
            return
        end
        newTrackIds = obj.divideTrack(trackId,segId,locId);
        setTracksTableData();
        setSegsTableData(newTrackIds,[]);
        setLocsTableData(newTrackIds,[]);
        saveAllPlotView()
        savedPlotInfo.trackId = newTrackIds;
        t=timer('StartDelay', 0.01, 'ExecutionMode','singleShot','TimerFcn',@restorePlotView_CB);
        start(t);
    end

    function joinSegments_CB(~,~)
        [trackId, segId] = getSelectedSegs();
        if ~all(trackId==trackId(1))
            error('Can only join segments from the same track');
        end
        obj.joinSegment(trackId(1),segId);
        setTracksTableData();
        setSegsTableData(trackId,[]);
        setLocsTableData(trackId,[]);
        saveAllPlotView()
        savedPlotInfo.trackId=trackId;
        t=timer('StartDelay', 0.01, 'ExecutionMode','singleShot','TimerFcn',@restorePlotView_CB);
        start(t);
    end



    function plotMainAxes(plotFun, varargin)
        %plot something into the main axes
        obj.assertInitialized()
        axes(handles.axes);
        cla(handles.axes,'reset');
        plotFun(varargin{:});
        savedPlotInfo.plotFun = plotFun;
        handles.axes.Clipping='on';
        axesMode_CB('datacursormode',[],[]); %let datacursormode be the default
    end
    
    function saveAllPlotView()
        savedPlotInfo.trackId=[];
        savedPlotInfo.CameraPosition = handles.axes.CameraPosition;
        savedPlotInfo.CameraTarget = handles.axes.CameraTarget;
        savedPlotInfo.CameraViewAngle = handles.axes.CameraViewAngle;
    end
    
    function restorePlotView_CB(~,~)
        plotMainAxes(savedPlotInfo.plotFun,savedPlotInfo.trackId);
        if isfield(savedPlotInfo,'CameraPosition')
            handles.axes.CameraPosition = savedPlotInfo.CameraPosition;
            handles.axes.CameraViewAngle = savedPlotInfo.CameraViewAngle;
            handles.axes.CameraTarget = savedPlotInfo.CameraTarget;
        end
        
    end

    
    %% List Callbacks allow the list of displayed Segs or Locs to be set
    function listSegments_CB(~,~)
        setSegsTableData(getSelectedTracks());
    end

    function listTrackLocalizations_CB(~,~)
        setLocsTableData(getSelectedTracks(),[]);
    end

    function listSegLocalizations_CB(~,~)
        [trackIds,segIds] = getSelectedSegs();
        setLocsTableData(trackIds,segIds);
    end

    %% Exporting
    function exportWrapper(appendix,title,exportFunc,varargin)
        %This is the generic function for all the export callbacks
        exppath = obj.getFilePath('Export');
        filepattern = [obj.saveFileBaseName appendix ExportExt];
        filepath = Pickle.selectUnusedFileName(exppath, filepattern,ExportFormats, title);
        if isempty(filepath); return; end
        obj.updateWaitbar(0,'Exporting...');
        try
            exportFunc(filepath,varargin{:});
        catch err
            obj.updateWaitbar(1);
            throw(err);
        end
        obj.updateWaitbar(1);
    end
    
    function exportAllStats_CB(~,~)
        exportWrapper('.AllStats', 'Export All Stats', @obj.exportAllStats);
    end
    
    function exportAllTrackStats_CB(~,~)
        exportWrapper('.AllTrackStats', 'Export All Track Stats', @obj.exportTrackStats);
    end

    function exportAllSegStats_CB(~,~)
        exportWrapper('.AllSegStats', 'Export All Segment Stats', @obj.exportSegmentStats);
    end

    function exportAllLocStats_CB(~,~)
        exportWrapper('.AllLocStats', 'Export All Localization Stats', @obj.exportLocalizationStats);
    end

    function exportTrackStats_CB(~,~)
        exportWrapper('.TrackStats', 'Export Selected Track Stats', @obj.exportTrackStats, getSelectedTracks());
    end

    function exportSegStats_CB(~,~)
        [trackIds,segIds] = getSelectedSegs();
        exportWrapper('.SegStats', 'Export Selected Segment Stats', @obj.exportSegmentStats, trackIds, segIds);
    end

    function exportTrackLocStats_CB(~,~)
        exportWrapper('.LocStats', 'Export Selected Track Localization Stats', @obj.exportLocalizationStats, getSelectedTracks());
    end

    function exportSegLocStats_CB(~,~)
        [trackIds,segIds] = getSelectedSegs();
        exportWrapper('.LocStats', 'Export Selected Segment Localization Stats', @obj.exportLocalizationStats, trackIds, segIds);
    end

    function exportLocStats_CB(~,~)
        [trackIds,segIds]=getSelectedLocs();
        exportWrapper('.LocStats', 'Export Selected Localization Stats', @obj.exportLocalizationStats, trackIds, segIds);
    end
end
