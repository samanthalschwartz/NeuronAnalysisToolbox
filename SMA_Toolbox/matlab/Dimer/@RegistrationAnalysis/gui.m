% RegistrationAnalysis.gui

function guiFig = gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end

    %% Constants and Global Initialization
    gui@GUIBuilder(obj); %Call GUIBuilder initialization
    gui_name = '[RegistrationAnalysis] 2-channel Registration Analysis';
    
    lastPath = obj.workingDir; %save this for file open dialogs after a reset
    DataTableFormat = struct(); % Format definition of data table
    MapsTableFormat = struct(); % Format definition of Map table
    uH = GUIBuilder.default_unitHeight; % unit height for elements
    boarder = GUIBuilder.default_boarder;%Boarder width around the outside of the gui
    sp = GUIBuilder.default_spacing; %spacing between elements.
    but_sz = GUIBuilder.default_buttonSize; %Button size
    
    fig_sz = [1200 820]; %figure size
    halfw = 650;
    dataPan_H = 250; % Height of data panel
    advancedParams_H = 200; % Height of advance params
    fullw_sz = [fig_sz(1)-2*boarder uH]; %size of a full-width component that stays within the boarders.
    halfw_sz = [halfw-sp-boarder uH];
    axes_pos = [halfw+sp boarder fig_sz(1)-halfw-boarder-sp fig_sz(2)-2*boarder]; % To be updated later

    halfPan_pos=[boarder boarder+halfw halfw_sz];
    fullPan_pos=[boarder boarder+halfw fullw_sz];

   
    %% Create figure and controls
    % Make figure
    guiFig = figure('Units','pixels','Position',[10 0 fig_sz],'Resize','off',...
                    'MenuBar','none','ToolBar','figure','NumberTitle','off',...
                    'Name',gui_name,'Visible','on',...
                    'CloseRequestFcn',@close_CB);
    obj.guiFig = guiFig;    
    whitebg(guiFig);

    handles.tabG = uitabgroup('Parent',guiFig,'Units','Pixels','Position',halfPan_pos-[0 0 boarder 0],'TabLocation','top');
    handles.tabs.general = uitab('Parent',handles.tabG ,'Title','General Params');
    handles.tabs.advanced = uitab('Parent',handles.tabG ,'Title','Advanced Params');
    handles.panels.maps = uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Maps');
    
    handles.panels.data = uipanel('Parent',guiFig,'Units','Pixels','Position',fullPan_pos,'Title','Data');
    handles.panels.file = uipanel('Parent',guiFig,'Units','Pixels','Position',fullPan_pos,'Title','File');
    
    handles.axes = axes('Parent',guiFig,'Units', 'Pixels','Visible','off','Position',axes_pos);
 
    
    createMenus();
    populateFilePanel();
    populateDataPanel();
    populateMapsPanel();
    populateGeneralTab();
    populateAdvancedTab();
    alignPanels();
    drawnow();
    positionAxes();
    
    function createMenus()
        labels = {'Load ...', 'Save','Save As ...','Clear All Data','Reset Object',[],'Quit'};
        CBs = {@load_CB, @save_CB, @saveAs_CB,@clearData_CB,@resetObject_CB,[],@close_CB};
        GUIBuilder.makeFigureMenu(guiFig,'File',labels,CBs);

        labels = {'View Sequence ...','View RAW Sequence ...','View Sequence Sturation ...',[],'View Selected Boxes...','View Frame Pairing...',[],...
                  'Plot Intensity','Plot Sigma','Plot LLH','Plot Positions','Plot Position SE',[],...
                  'Reprocess Sequence','Reprocess ALL Sequences',[], 'Delete Sequence', 'Clear ALL Sequences'};
        CBs = {@viewSequence_CB,@viewRawSequence_CB,@viewSequenceSaturation_CB,[],@viewBoxes_CB,@viewFramePairing_CB,[],...
               @plotIntensity_CB,@plotSigma_CB,@plotLLH_CB,@plotPosition_CB,@plotPositionSE_CB,[],...
               @reprocessData_CB,@reprocessAllData_CB,[], @deleteData_CB, @clearData_CB};
        handles.contextmenus.data = GUIBuilder.makeJavaContextMenu(labels, CBs);
        
        labels = {'Plot Map Error', 'View Map Error ...', 'View Map Displacement ...', 'View Map Vectors ...',[], ...
                  'Set Training Frames', 'Set Testing Frames', 'Export map function', [],...
                  'Reprocess Map(s)', 'Reprocess ALL Maps', [], 'Delete Map', 'Clear All Maps'};
              
        CBs = {@plotMapError_CB, @viewMapError_CB, @viewMapDisplacement_CB, @viewMapVectors_CB, [],...
               @setMapTrainingData_CB, @setMapTestingData_CB, @exportMapFunction_CB, [], ...
               @reprocessMaps_CB, @reprocessAllMaps_CB, [], @deleteMap_CB, @clearMaps_CB};
       handles.contextmenus.maps = GUIBuilder.makeJavaContextMenu(labels, CBs);
        
        labels = {'Reprocess All Data'};
        CBs = {@reprocessAllData_CB};
        handles.contextmenues.advancedParams = GUIBuilder.makeJavaContextMenu(labels, CBs);
    end

    function populateFilePanel()
        % File Panel
        % Bottom row of buttons
        bot_row.pos = [sp, sp, fullw_sz-2*boarder-2*sp];
        bot_row.names = {'Load', 'Save','SaveAs','Clear Data','Set Defaults'};
        bot_row.CBs = {@load_CB, @save_CB,@saveAs_CB,@clearData_CB, @setDefaultProps_CB};
        GUIBuilder.buttonRow(handles.panels.file, bot_row.pos, but_sz+[20 0], bot_row.names, bot_row.CBs);
        %Quit is out to the side so we don't include it.
        handles.quitButton = uicontrol('Parent',handles.panels.file,'Style','pushbutton','String','Quit',...
            'Position',[fullw_sz(1)-3*sp-but_sz(1) sp but_sz],...
            'Callback',@close_CB);
        edits.pos = bot_row.pos + [0, but_sz(2)+3*sp, 0, 2*sp];
        edits.hNames = {'saveFile'}; %names in the handles structure returned
        edits.labels = {'Save File:'};
        edits.values = {obj.saveFilePath};
        edits.CBs = {@saveAs_CB};
        handles.file = GUIBuilder.labeledHEdits(handles.panels.file, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
    end

    function populateDataPanel()
        table_pos = [sp sp fullw_sz(1) dataPan_H-20];
        col_names={'Frame<br />Id','Date','Valid','#Ref<br />Pairs', '#LeftCh.<br />Emitters','#RightCh.<br />Emitters',...
                   '#Saturated<br />Pixels', 'RMS Disp.<br />(um)','Max Disp.<br />(um)',...
                   'SourceName'};
        DataTableFormat.names =['Sequence', cellmap(@(t) sprintf('<html><div align="center">%s</div></html>',t) ,col_names)]; %First column can't be html
        DataTableFormat.formats={'','', 'char','', '', '', '', '', '','',''};
        data = obj.getDataCellArray();
        [handles.dataTable, handles.dataTablePanel] = ...
                GUIBuilder.makeTreeTable(handles.panels.data,DataTableFormat, data, table_pos, 'Groupable', true);
        handles.dataTable.MousePressedCallback = {@GUIBuilder.javaTableMouse_CB,handles.contextmenus.data};
    end

    function populateMapsPanel()
        
        table_pos = [sp sp halfw_sz(1) uH];
        col_names={'Name','Algorithm','Valid','RMSE<br />(nm)','MaxErr<br />(nm)','#Train<br />Points','#Test<br />Points',...
                    '#Train<br />Frames','#Test<br />Frames','Train<br />Data Idxs','Test<br />Data Idxs'};
        MapsTableFormat.names = cellmap(@(t) sprintf('<html><div align="center">%s</div></html>',t) ,col_names);
        MapsTableFormat.formats={'char','char', '', '', '', '', '', '', '', 'char', 'char'};
        data = obj.getMapsCellArray();
        [handles.mapsTable, handles.mapsTablePanel] = ...
                GUIBuilder.makeTreeTable(handles.panels.maps, MapsTableFormat, data, table_pos, 'Groupable', false);
        handles.mapsTable.MousePressedCallback = {@GUIBuilder.javaTableMouse_CB,handles.contextmenus.maps};
    end

    function populateGeneralTab()
        % Bottom row of buttons
        bot_row.pos = [sp, sp, halfw_sz(1)-boarder, uH];
        bot_row.names = {'Calibrate'};
        bot_row.CBs = {@gainCal_CB};
        GUIBuilder.buttonRow(handles.tabs.general, bot_row.pos, but_sz, bot_row.names, bot_row.CBs);
        
        edits.pos = [sp, sp+but_sz(2)+sp, halfw_sz(1)-boarder-2*sp, uH];
        edits.hNames = {'FrameSize','PixelSize', 'CCDGain','CCDBackground'}; %names in the handles structure returned
        edits.labels = {'Sensor Size [X,Y] (pixels):','Pixel Size [X,Y] (um):','CCD Gain (e-/ADU):', 'CCD Background (ADU):'};
        edits.values = {flip(obj.frameSize), obj.pixelSize, obj.CCDGain, obj.CCDBackground};
        edits.CBs = {[], @setPixelSize_CB, @setCCDGain_CB, @setCCDBackground_CB};
        handles.general = GUIBuilder.labeledHEdits(handles.tabs.general, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
    end

    function populateAdvancedTab()
        params_pos = [sp, sp, halfw-2*boarder-2*sp, advancedParams_H];
        params.BoxSize.disp = 'Box Size';
        params.BoxSize.desc = 'Size of emitter boxes [X, Y]';
        params.MinBoxes.disp = 'Min Number of Boxes to Find';
        params.MinBoxes.desc = 'Try to find at least this many boxes in each frame (use if not enough boxes found) [0=disable].';
        params.NeighborhoodSize.disp='Neighborhood Size';
        params.NeighborhoodSize.desc='Size of filtering Neighborhood (>=3)';
        params.NeighborhoodSize.range =  {int32(3), int32(5),int32(7),int32(9),int32(11),int32(13),int32(15)};
        params.PSFSigma.disp='PSF Sigma';
        params.PSFSigma.desc='PSF Gaussian Sigma (pixels) [X Y]';
        params.EstimationMethod.disp = 'Estimation Method';
        params.EstimationMethod.desc = 'Gaussian fitting maximization method for MLE estimation.';
        params.EstimationMethod.range = Gauss2DsMAP.EstimationMethods;
        params.BoundaryWidth.disp  = 'Boundary Width';
        params.BoundaryWidth.desc  = 'Width in pixels around each channel to automatically reject fits (Set to 0 for no effect)';
        params.PixelSaturationThreshold.disp = 'Pixel Saturation Threshold';
        params.PixelSaturationThreshold.desc = 'Threshold in ADUs above which the pixel has been saturated';
        params.MinIntensity.disp = 'Min Intensity';
        params.MinIntensity.desc = 'Minimum Intensity to accept for fits';
        params.MinSigma.disp = 'Min Sigma';
        params.MinSigma.desc = 'Minimum Gaussian Sigma (pixels) to accept.';
        params.MaxSigma.disp = 'Max Sigma';
        params.MaxSigma.desc = 'Maximum Gaussian Sigma (pixels) to accept.';
        params.MaxPositionSE.disp = 'Max Position SE';
        params.MaxPositionSE.desc = 'Maximum Estimated error in fit for either X or Y position to accept (pixels).';
        params.MinLocalizationDist.disp = 'Min Localization Dist';
        params.MinLocalizationDist.desc = 'Minimum acceptable localization distance to allow (pixels).';
        params.MaxRefDisplacement.disp = 'Max Reference Displacement';
        params.MaxRefDisplacement.desc = 'Maximum allowable displacement between overlayed images (pixels).  This prevents pairing the wrong reference points.';
        handles.advancedParamsTable = GUIBuilder.makePropertyGrid(handles.tabs.advanced,...
                                        'Advanced Parameters',obj.Params,params, params_pos);
        % Set up the java menu
        GUIBuilder.setPropertyGridMenu(handles.advancedParamsTable,...
                    {@GUIBuilder.javaTableMouse_CB, handles.contextmenues.advancedParams});
    end
    function seqIds = getSelectedSequence()
        tabH = handles.dataTable;
        treeModelH = tabH.getModel().getActualModel();
        tableModelH = treeModelH.getActualModel(); %This model does not have the tree columns collapsed
        selectedRows = tabH.getSelectedRows();
        actualRows = arrayfun(@(r) treeModelH.getActualRowAt(r), selectedRows);
        headerRows = selectedRows(actualRows<0);
        actualRows = actualRows(actualRows>=0);
        headerIdxs = cellfun(@(v) str2double(v{2}), cellmap(@(r) strsplit(char(treeModelH.getValueAt(r,0)),':'), headerRows));
        seqIds = [headerIdxs; arrayfun(@(r) tableModelH.getValueAt(r,0),actualRows)];
        seqIds = unique(seqIds);
    end

    function frameIds = getSelectedFrames()
        tabH = handles.dataTable;
        treeModelH = tabH.getModel().getActualModel();
        tableModelH = treeModelH.getActualModel(); %This model does not have the tree columns collapsed
        selectedRows = tabH.getSelectedRows();
        actualRows = arrayfun(@(r) treeModelH.getActualRowAt(r), selectedRows);
        actualRows = actualRows(actualRows>=0);
        frameIds = [arrayfun(@(r) tableModelH.getValueAt(r,0),actualRows), arrayfun(@(r) tableModelH.getValueAt(r,1),actualRows)];
    end

    function mapIds = getSelectedMaps()
        tabH = handles.mapsTable;
        treeModelH = tabH.getModel().getActualModel();
        selectedRows = tabH.getSelectedRows();
        actualRows = arrayfun(@(r) treeModelH.getActualRowAt(r), selectedRows);
        actualRows = actualRows(actualRows>=0);
        mapIds = actualRows+1;
    end

    function setSequenceTableData()
        data=obj.getDataCellArray();
        GUIBuilder.setTreeTableData(handles.dataTable, DataTableFormat, data,1);
    end

   function setMapsTableData()        
        data = obj.getMapsCellArray();
        GUIBuilder.setTreeTableData(handles.mapsTable, MapsTableFormat, data);
    end

    function positionAxes()
        labelmargin = [30 20 0 20];
        dataPan_pos = handles.panels.data.Position;
        axes_pos = [halfw+sp, sum(dataPan_pos([2,4]))+sp, 0, 0];
        axes_pos(3:4) = fig_sz - boarder - axes_pos(1:2);
        GUIBuilder.positionAxes(handles.axes, axes_pos, labelmargin);
        axes.Clipping = 'on';
        resetAxes();
    end

    function plotMainAxes(plotFun, varargin)
        %plot something into the main axes
        obj.assertInitialized();
        axes(handles.axes);
        cla(handles.axes,'reset');
        plotFun(varargin{:});
        handles.axes.Clipping='on';
    end

    function resetAxes()
        if obj.NData<1
            cla(handles.axes,'reset');
            handles.axes.Visible='off';
        elseif obj.NMaps>1 && any([obj.maps.valid])
            plotMainAxes(@obj.plotMapError, handles.axes);
        else
            plotMainAxes(@obj.plotEmitterIntensity, handles.axes);
        end
    end

    function panels_pos = alignPanels()
        filePan_pos = GUIBuilder.autoSizePanel(handles.panels.file);
        filePan_pos(2) = boarder;
        handles.panels.file.Position = filePan_pos;

        tabG_pos = GUIBuilder.autoSizeTabGroup(handles.tabG);
        dataPan_pos = fullPan_pos;
        dataPan_pos(4) = dataPan_H;
        handles.panels.data.Position = dataPan_pos;
        mapsPan_H = fig_sz(2)-2*boarder-2*sp - dataPan_H - tabG_pos(4) -filePan_pos(4);
        mapsPan_pos = halfPan_pos;
        mapsPan_pos(4) = mapsPan_H;
        handles.panels.maps.Position = mapsPan_pos;
        handles.mapsTablePanel.Position = [sp, sp, mapsPan_pos(3:4)-2*sp] - [0 0 0 16];
        pans = {handles.panels.file, handles.panels.data, handles.panels.maps, handles.tabG};
        panels_pos=GUIBuilder.align(pans,'Left','Fixed',sp);
    end

    function updateAdvancedParams()
        %Update the values in the advanced params PropertyGrid object using the obj.Params struct
        H = handles.advancedParamsTable;
        pgf = H.Properties;
        fns = fieldnames(obj.Params);
        for i = 1:length(fns)
            fn = fns{i};
            typ = pgf.FindByName(fn).Type.PrimitiveType;
            switch typ
                case 'int32'
                    pgf.FindByName(fn).Value = int32(obj.Params.(fn));
                otherwise
                    pgf.FindByName(fn).Value = obj.Params.(fn);
            end
        end
        H.Properties = pgf;
    end

    function reloadAdvancedParams()
        new_params = handles.advancedParamsTable.GetPropertyValues();
        nfs = fieldnames(new_params);
        for i = 1:length(nfs)
            obj.Params.(nfs{i}) = new_params.(nfs{i});
        end
    end

    function reloadAllControls()
        %This is a full-reload of all controls.  Calls updateControls() wich is the 'softer' update
        if isempty(obj.workingDir) || ~exist(obj.workingDir,'dir')
            lastPath = pwd;
        else
            lastPath = obj.workingDir; %save this for file open dialogs after a reset
        end
        updateAdvancedParams();
        resetDataControls()
    end
    
    function resetDataControls()
        setSequenceTableData();
        setMapsTableData();
        updateControls();
        resetAxes();
    end

%     function recalibrateData(gain, bg)
% 
%     end

    function updateControls()
        handles.file.saveFile.String = obj.saveFilePath;
        if obj.dirty
            handles.file.saveFile.BackgroundColor = obj.color_unsavedBG;
        else
            handles.file.saveFile.BackgroundColor = obj.color_editBG;
        end
        handles.general.FrameSize.String = arr2str(obj.frameSize);
        obj.setDefaultableControl(handles.general.PixelSize,'pixelSize');
        obj.setDefaultableControl(handles.general.CCDBackground,'CCDBackground');
        obj.setDefaultableControl(handles.general.CCDGain,'CCDGain');
    end

    function load_CB(~, ~)
        if ~isempty(obj.workingDir)
            lastPath = obj.workingDir;
        end
        filepaths = Pickle.selectExistingFileName(lastPath, [], obj.LoadableDataFormats, 'Load RegisteredPairsAnalysis','MultiSelect','on');
        if isempty(filepaths); return; end
        
        obj.load(filepaths);
        obj.guiFig = guiFig;
        obj.inGui = true;
        reloadAllControls();
    end

    function save_CB(~,~)
        reloadAdvancedParams();
        obj.updateWaitbar(0,'Saving...');
        obj.save();
        updateControls();
        obj.updateWaitbar(1);
    end

    function saveAs_CB(~,~)
        reloadAdvancedParams();
        F = handles.file.saveFile.String;
        title = 'Save RegistrationAnalysis As ...';
        [file,path] = uiputfile(obj.SaveableDataFormats,title,F);
        if ~file
            return;
        end
        obj.saveas(fullfile(path,file));
        reloadAllControls();
    end

    function resetObject_CB(~,~)
        if obj.initialized
            lastPath=obj.workingDir; %Only change the lastPath if we actually have a valid workingDir
        end
        obj.resetObject();        
        obj.guiFig = guiFig;
        obj.inGui = true;
        reloadAllControls();
    end

    function clearData_CB(~,~)
        obj.clearData();
        reloadAllControls();
    end

    function deleteData_CB(~,~)
        obj.deleteData(getSelectedSequence());
        reloadAllControls();
    end

    function reprocessData_CB(~,~)
        reloadAdvancedParams();
        obj.reprocessData(getSelectedSequence());
        reloadAllControls();
    end

    function reprocessAllData_CB(~,~)
        reloadAdvancedParams();
        obj.reprocessData();
        reloadAllControls();
    end

    function deleteMap_CB(~,~)
        obj.deleteMap(getSelectedMaps());
        reloadAllControls();
    end

    function clearMaps_CB(~,~)
        obj.clearMaps();
        reloadAllControls();
    end

    function reprocessMaps_CB(~,~)
        obj.reprocessMaps(getSelectedMaps());
        reloadAllControls();
    end

    function reprocessAllMaps_CB(~,~)
        reloadAdvancedParams();
        obj.reprocessMaps();
        reloadAllControls();
    end

    function close_CB(~,~)
        try
            GUIBuilder.clearPropertyGridMenus(handles.advancedParamsTable);
            handles.dataTable.MousePressedCallback=[];
            handles.mapsTable.MousePressedCallback=[];
            GUIBuilder.clearJavaContextMenu(handles.contextmenus.data);
            GUIBuilder.clearJavaContextMenu(handles.contextmenus.maps);
            GUIBuilder.clearJavaContextMenu(handles.contextmenues.advancedParams);
        catch err
            disp(getReport(err));
        end
        obj.closeGUI();
    end

    function setPixelSize_CB(H,~)
        pixelsize = str2num(H.String); %#ok<ST2NM>
        if (isscalar(pixelsize) || length(pixelsize)==2) && all(isfinite(pixelsize)) && all(pixelsize>0)
            obj.pixelSize=pixelsize([1,end]); %make it a [X, Y]
            obj.dirty=true;
        else
            H.String = arr2str(obj.pixelSize); %Reset the string to the original
            warning('SPData:gui','Invalid pixel size: %s',H.String);
        end
        obj.setDefaultableControl(handles.general.PixelSize,'pixelSize');
    end

    function setCCDBackground_CB(H,~)
        bg = str2double(H.String);
        if ~isfinite(bg) || bg<0
            H.String = obj.CCDBackground;
            error('RegistrationAnalysis:gui','Invalid CCD Background: %s', H.String);
        elseif ~isempty(obj.CCDGain) % If Gain is also set we are ready to calibrate
            obj.calibrateCCD(obj.CCDGain, bg);
            resetDataControls();
        else
            obj.CCDBackground = bg; %Just set the background because gain is not set yet.
        end
    end

    function setCCDGain_CB(H,~)
        gain = str2double(H.String);
        if ~isfinite(gain) || gain<=0
            H.String = obj.CCDGain;
            error('RegistrationAnalysis:gui','Invalid CCD Gain: %s', H.String);
        elseif ~isempty(obj.CCDBackground) % If BG is also set we are ready to calibrate
            obj.calibrateCCD(gain, obj.CCDBackground);
            resetDataControls();
        else
            obj.CCDGain = gain; %Just set the gain because background is not set yet.
        end
    end

    function gainCal_CB(~,~)
        obj.gainCalGui(@resetDataControls);
    end

    function viewSequence_CB(~,~)
        seqIdx = getSelectedSequence();
        if ~isempty(seqIdx)
            obj.viewSequence(seqIdx(1));
        end
    end

    function viewRawSequence_CB(~,~)
        seqIdx = getSelectedSequence();
        if ~isempty(seqIdx)
            obj.viewRawSequence(seqIdx(1));
        end
    end

    function viewSequenceSaturation_CB(~,~)
        seqIdx = getSelectedSequence();
        if ~isempty(seqIdx)
            obj.viewSequenceSaturation(seqIdx(1));
        end
    end

    function viewBoxes_CB(~,~)
        seqIdx = getSelectedSequence();
        if ~isempty(seqIdx)
            obj.viewSelectedBoxes(seqIdx(1));
        end
    end

    function viewFramePairing_CB(~,~)
        idx = getSelectedFrames();
        if ~isempty(idx)
            obj.viewFramePairing(idx(1,1),idx(1,2));
        end
    end

    function plotIntensity_CB(~,~)
        plotMainAxes( @obj.plotEmitterIntensity, handles.axes, getSelectedSequence());
    end

    function plotSigma_CB(~,~)
        plotMainAxes( @obj.plotEmitterSigma, handles.axes, getSelectedSequence());
    end

    function plotLLH_CB(~,~)
        plotMainAxes( @obj.plotEmitterLLH, handles.axes, getSelectedSequence());
    end

    function plotPosition_CB(~,~)
        plotMainAxes( @obj.plotEmitterPosition, handles.axes, getSelectedSequence());
    end

    function plotPositionSE_CB(~,~)
        plotMainAxes( @obj.plotEmitterPositionSE, handles.axes, getSelectedSequence());
    end

    function plotMapError_CB(~,~)
        plotMainAxes(@obj.plotMapError, handles.axes, getSelectedMaps());
    end

    function viewMapError_CB(~,~)
        idx = getSelectedMaps();
        if ~isempty(idx)
            obj.plotMapError([],idx(1));
        end
    end

    function viewMapDisplacement_CB(~,~)
        idx = getSelectedMaps();
        if ~isempty(idx)
            obj.viewMapDisplacement(idx);
        end
    end

    function viewMapVectors_CB(~,~)
        idx = getSelectedMaps();
        if ~isempty(idx)
            obj.viewMapVectors(idx);
        end
    end

end

