% PairAnalysis.gui

function guiFig = gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end

    %% Constants and Global Initialization
    gui@GUIBuilder(obj); %Call GUIBuilder initialization
    gui_name = '[PairAnalysis] Pair Interaction Analysis Class';
    
    lastPath = obj.workingDir; %save this for file open dialogs after a reset
    
    uH = GUIBuilder.default_unitHeight; % unit height for elements
    boarder = GUIBuilder.default_boarder;%Boarder width around the outside of the gui
    sp = GUIBuilder.default_spacing; %spacing between elements.
    but_sz = GUIBuilder.default_buttonSize; %Button size
    
    fig_sz = [1280 800]; %figure size
    halfw = 560;
    fullw_sz = [fig_sz(1)-2*boarder-2*sp uH]; %size of a full-width component that stays within the boarders.
    halfw_sz = [halfw-sp-boarder uH];
    axes_pos = [halfw boarder fig_sz(1)-halfw fig_sz(2)-2*boarder]; % To be updated later

    %% Create figure and controls
    % Make figure
    guiFig = figure('Units','pixels','Position',[10 0 fig_sz],'Resize','off',...
                    'MenuBar','none','ToolBar','figure','NumberTitle','off',...
                    'Name',gui_name,'Visible','on',...
                    'CloseRequestFcn',@close_CB);
    obj.guiFig = guiFig;    
    whitebg(guiFig);
    
    
    fullPan_pos = [boarder boarder fullw_sz];
    halfPan_pos=[boarder boarder+halfw halfw_sz];
    tabs_pos=[boarder boarder fullw_sz(1) fig_sz(2)-axes_pos(4)-axes_pos(2)];
    handles.tabG = uitabgroup('Parent',guiFig,'Units','Pixels','Position',tabs_pos,'TabLocation','top');
    handles.tabs.file=uitab('Parent',handles.tabG ,'Title','File Associations');
    handles.tabs.locs=uitab('Parent',handles.tabG ,'Title','Pair Localizations');
    handles.axes= axes('Units', 'Pixels','Visible','off');
    handles.panels.pairs = uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Pairs');
    
    PairsTableColFormat = struct();

    savedPlotInfo = struct();
    
    createMenus();
    populatePairPanel(handles.panels.pairs);
    populateFilePanel(handles.tabs.file);
    alignPanels();
    drawnow();
    positionAxes(handles.axes);
    reloadAllControls();
   
    function populatePairPanel(panH)
        table_pos = [sp sp halfw_sz(1)-2*sp 400];
        col_names={'Pair<br />Id','Track1<br />Id','Track2<br />Id', 'Num<br />Locs',...
                   'Duration<br />(frames)', 'Duration<br />(s)', 'Min<br />Dist(um)',...
                   'Dmle1<br />(um^2/s)','Dmle2<br />(um^2/s)'};
        col_names=cellmap(@(t) sprintf('<html><div align="center">%s</div></html>',t) ,col_names);
        col_formats=cell(1,numel(col_names));
        col_editable=false(1,numel(col_names));
        data = cell(1,numel(col_names));
        col_widths={20, 20, 20, 40, 40,'auto','auto','auto','auto'};
        PairsTableColFormat.names = col_names;
        PairsTableColFormat.formats = col_formats;
        PairsTableColFormat.widths = col_widths;
        PairsTableColFormat.editable = col_editable;
        [handles.pairsTable, handles.pairsTablePanel] = ...
                GUIBuilder.makeTreeTable(panH,PairsTableColFormat, data, table_pos, 'Groupable', false);
        handles.pairsTable.MousePressedCallback = {@GUIBuilder.javaTableMouse_CB, handles.contextmenus.pairs};
        setPairsTableData();
    end
    
    function populateFilePanel(panH)
        % File Panel
        % Bottom row of buttons
        bot_row.pos = [sp, sp, fullw_sz];
        bot_row.names = {'Load', 'Save','SaveAs','BatchProcess'};
        bot_row.CBs = {@load_CB, @save_CB,@saveAs_CB,  @batchProcess_CB};
        GUIBuilder.buttonRow(panH, bot_row.pos, but_sz+[20 0], bot_row.names, bot_row.CBs);
        %Quit is out to the side so we don't include it.
        handles.quitButton = uicontrol('Parent',panH,'Style','pushbutton','String','Quit',...
            'Position',[fullw_sz(1)-but_sz(1) sp but_sz],...
            'Callback',@close_CB);
        edits.pos = bot_row.pos + [0, but_sz(2)+sp, 0, 2*sp];
        edits.hNames = {'spdata','rpt_channel1', 'rpt_channel2', 'regAnalysis','saveFile'}; %names in the handles structure returned
        edits.labels = {'SPData (.spdata) Path:',  'RPT Channel1 (.rpt) Path:', 'RPT Channel2 (.rpt) Path:', 'RegistrationAnalysis (.reganalysis) Path:', 'Save File (.regpairs) Path:'};
        edits.values = {obj.getFilePath('data'),obj.getFilePath('rpt_channel1'),obj.getFilePath('rpt_channel2'),obj.getFilePath('regAnalysis') ,obj.saveFilePath };
        edits.CBs = {};
        handles.file = GUIBuilder.labeledVEdits(panH, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
        
        open_CB_list = {@openSPData_CB, @openRPTChannel1_CB, @openRPTChannel2_CB, @openRegAnalysis_CB};
        open_string_list = {'Open .spdata', 'Open .rpt', 'Open .rpt', 'Open .reganalysis'};
        load_CB_list = {@load_SPData_CB, [], [], @change_registration_CB};
        load_string_list = {'Load .spdata', [], [], 'Load .reganalysis'};
        for i=1:numel(edits.hNames);
            name = edits.hNames{i};
            editH = handles.file.(name);
            pos = editH.Position;
            load_pos = pos;
            pos(1) = pos(1) + 2*(but_sz(1)+sp);
            pos(3) = pos(3) - 2*(but_sz(1)+sp);
            editH.Position = pos;
            if i==numel(edits.hNames)
                break;
            end
            load_pos(3) = but_sz(1);
            open_pos = load_pos;
            open_pos(1) = open_pos(1) + but_sz(1) + sp;
            if ~isempty(load_CB_list{i})
                handles.fileButtons.(name) = uicontrol('Parent',panH,'Units','Pixels','Position',load_pos,'Style','pushbutton',...
                        'String',load_string_list{i}, 'Callback', load_CB_list{i});
            end
            handles.fileButtons.(name) = uicontrol('Parent',panH,'Units','Pixels','Position',open_pos,'Style','pushbutton',...
                    'String',open_string_list{i}, 'Callback', open_CB_list{i});
        end
    end

    function createMenus()
        %File Menu

        labels = {'Load ...','Load from SPData ...', 'Save','Save As ...','Reload Ch. Reg ...','Reset Object',[],'Quit'};
        CBs = {@load_CB, @load_SPData_CB, @save_CB, @saveAs_CB, @change_registration_CB, @resetObject_CB,[],@close_CB};
        GUIBuilder.makeFigureMenu(guiFig,'File',labels,CBs);

        labels = {'Plot ALL tracks Wavelength', 'Plot ALL tracks [Speed]', 'Plot ALL tracks [Temporal]', 'Plot ALL tracks [Sequence]',[],...
                  'View ALL tracks [Wavelength] ...', 'View ALL tracks [Speed] ...', 'View ALL tracks [Temporal] ...', 'View ALL tracks [Sequence] ...',[],...
                  'View ALL Track Movie [Wavelength] ...', 'View ALL Track Movie [Speed] ...','View ALL Track Movie [Temporal] ...','View ALL Track Movie [Sequence] ...'};
        CBs = {@plotAllTracksWavelength_CB, @plotAllTracksSpeed_CB, @plotAllTracksTemporal_CB, @plotAllTracksSequence_CB, [], ...
               @viewAllTracksWavelength_CB, @viewAllTracksSpeed_CB, @viewAllTracksTemporal_CB, @viewAllTracksSequence_CB, [], ...
               @viewAllTrackMovieWavelength_CB, @viewAllTrackMovieSpeed_CB, @viewAllTrackMovieTemporal_CB, @viewAllTrackMovieSequence_CB};
        GUIBuilder.makeFigureMenu(guiFig,'View',labels,CBs);
        
        %Pairs context menu
        labels = {'List Localizations',[],...
                  'Plot Selected Pairs [Wavelength]', 'Plot Selected Pairs [Speed]', 'Plot Selected Pairs [Temporal]', 'Plot Selected Pairs [Sequence]',[],...
                  'View Selected Pairs [Wavelength] ...', 'View Selected Pairs [Speed] ...', 'View Selected Pairs [Temporal] ...','View Selected Pairs [Sequence] ...',[],...
                  'View Selected Pairs Movie [Wavelength] ...','View Selected Pairs Movie [Speed] ...','View Selected Pairs Movie [Temporal] ...','View Selected Pairs Movie [Sequence] ...',[],...
                  'Delete Pair'};
        CBs = {@listPairLocalizations_CB,[],...
               @plotPairWavelength_CB, @plotPairSpeed_CB, @plotPairTemporal_CB, @plotPairSequence_CB, [],...
               @viewPairWavelength_CB, @viewPairSpeed_CB, @viewPairTemporal_CB, @viewPairSequence_CB, [],...
               @viewPairMovieWavelength_CB, @viewPairMovieSpeed_CB, @viewPairMovieTemporal_CB, @viewPairMovieSequence_CB, [],...
               @deletePair_CB};
        handles.contextmenus.pairs = GUIBuilder.makeJavaContextMenu(labels, CBs);
    end

    function positionAxes(axesH)
        pos = handles.tabG.Position;
        axes_pos(2) = pos(2)+pos(4)+sp;
        axes_pos(4) = fig_sz(2)-axes_pos(2)-boarder;
        labelmargin = [20 20 20 0];
        GUIBuilder.positionAxes(axesH, axes_pos, labelmargin);
        axesH.Clipping = 'on';
    end

    function alignPanels()
%         pans = fieldnames(handles.panels);
%         for n = 1:length(pans)
%             GUIBuilder.autoSizePanel(handles.panels.(pans{n}));
%         end
%         panels_pos = GUIBuilder.align(struct2cell(handles.panels),'Left','Fixed',sp);
        GUIBuilder.autoSizeTabGroup(handles.tabG);
        pos = handles.tabG.Position;
        tabh = pos(2)+pos(4)+sp;
        pan_pos = [boarder, tabh, halfw_sz(1)-4*sp, fig_sz(2)-tabh-2*sp];
        handles.panels.pairs.Position = pan_pos;
        handles.pairsTablePanel.Position = [sp, sp, pan_pos(3)-3*sp, pan_pos(4)-2*sp-20]; 
    end

    function setPairsTableData()       
        if obj.nPairs==0
            data = cell(1,numel(PairsTableColFormat.names));
        else
            data = num2cell([(1:obj.nPairs)', obj.pairIds(:,1), obj.pairIds(:,2), [obj.pairStats(:).length]',...
                    [obj.pairStats(:).duration_frames]', [obj.pairStats(:).duration_secs]', [obj.pairStats(:).min_distance]',...
                    [obj.pairStats(:).Dmle1]',[obj.pairStats(:).Dmle2]']);
        end
        GUIBuilder.setTreeTableData(handles.pairsTable, PairsTableColFormat, data);
    end
    
    function pairIdxs = getSelectedPairs()
        %Get the selected pairIdxs
        obj.assertInitialized()
        tabH=handles.pairsTable;
        pairIdxs = arrayfun(@(r) tabH.getValueAt(r,0), tabH.getSelectedRows());
    end

    function plotMainAxes(plotFun, varargin)
        %plot something into the main axes
        obj.assertInitialized()
        axes(handles.axes);
        cla(handles.axes,'reset');
        plotFun(varargin{:}); % call arbitrary plotting function given
        savedPlotInfo.plotFun = plotFun;
        handles.axes.Clipping='on';
    end

    function plotAllPairsWavelength_CB(~,~)
        plotMainAxes(@(~) obj.plotAllPairsWavelength, handles.axes);
    end

    function plotAllTracksWavelength_CB(~,~)
        plotMainAxes(@obj.plotTracks, handles.axes, [], 'Wavelength');
    end

    function plotAllTracksSpeed_CB(~,~)
        plotMainAxes(@obj.plotTracks, handles.axes, [], 'Speed');
    end

    function plotAllTracksTemporal_CB(~,~)
        plotMainAxes(@obj.plotTracks, handles.axes, [], 'Temporal');
    end

    function plotAllTracksSequence_CB(~,~)
        plotMainAxes(@obj.plotTracks, handles.axes, [], 'Sequence');
    end

    function viewAllTracksWavelength_CB(~,~)
        obj.plotTracks([], [], 'Wavelength');
    end

    function viewAllTracksSpeed_CB(~,~)
        obj.plotTracks([], [], 'Speed');
    end

    function viewAllTracksTemporal_CB(~,~)
        obj.plotTracks([], [], 'Temporal');
    end

    function viewAllTracksSequence_CB(~,~)
        obj.plotTracks([], [], 'Sequence');
    end

    function viewAllTrackMovieWavelength_CB(~,~)
        obj.viewTracksMovie([], 'Wavelength');
    end

    function viewAllTrackMovieSpeed_CB(~,~)
        obj.viewTracksMovie([], 'Speed');
    end

    function viewAllTrackMovieTemporal_CB(~,~)
        obj.viewTracksMovie([], 'Temporal');
    end

    function viewAllTrackMovieSequence_CB(~,~)
        obj.viewTracksMovie([], 'Sequence');
    end

    function plotPairWavelength_CB(~,~)
        plotMainAxes(@obj.plotPairs, handles.axes, getSelectedPairs(), 'Wavelength');
    end

    function plotPairSequence_CB(~,~)
        plotMainAxes(@obj.plotPairs, handles.axes, getSelectedPairs(), 'Sequence');
    end

    function plotPairSpeed_CB(~,~)
        plotMainAxes(@obj.plotPairs, handles.axes, getSelectedPairs(), 'Speed');
    end

    function plotPairTemporal_CB(~,~)
        plotMainAxes(@obj.plotPairs, handles.axes, getSelectedPairs(), 'Temporal');
    end

    function viewPairWavelength_CB(~,~)
        obj.plotPairs([], getSelectedPairs(), 'Wavelength');
    end

    function viewPairSequence_CB(~,~)
       obj.plotPairs([], getSelectedPairs(), 'Sequence');
    end

    function viewPairSpeed_CB(~,~)
        obj.plotPairs([], getSelectedPairs(), 'Speed');
    end

    function viewPairTemporal_CB(~,~)
        obj.plotPairs([], getSelectedPairs(), 'Temporal');
    end

    function viewPairMovieWavelength_CB(~,~)
        obj.viewPairMovie(getSelectedPairs(), 'Wavelength');
    end

    function viewPairMovieSpeed_CB(~,~)
        obj.viewPairMovie(getSelectedPairs(), 'Speed');
    end

    function viewPairMovieTemporal_CB(~,~)
        obj.viewPairMovie(getSelectedPairs(), 'Temporal');
    end

    function viewPairMovieSequence_CB(~,~)
        obj.viewPairMovie(getSelectedPairs(), 'Sequence');
    end


%     function saveAllPlotView()
%         savedPlotInfo.trackId=[];
%         savedPlotInfo.CameraPosition = handles.axes.CameraPosition;
%         savedPlotInfo.CameraTarget = handles.axes.CameraTarget;
%         savedPlotInfo.CameraViewAngle = handles.axes.CameraViewAngle;
%     end
%     
%     function restorePlotView_CB(~,~)
%         plotMainAxes(savedPlotInfo.plotFun,savedPlotInfo.trackId);
%         if isfield(savedPlotInfo,'CameraPosition')
%             handles.axes.CameraPosition = savedPlotInfo.CameraPosition;
%             handles.axes.CameraViewAngle = savedPlotInfo.CameraViewAngle;
%             handles.axes.CameraTarget = savedPlotInfo.CameraTarget;
%         end        
%     end

    function openSPData_CB(~,~)
        data = obj.getData();
        if ~isempty(data)
            data.gui();
        end
    end

    function openRPTChannel1_CB(~,~)
        if ~isempty(obj.rpt_channel1)
            obj.rpt_channel1.gui();
        end
    end

    function openRPTChannel2_CB(~,~)
        if ~isempty(obj.rpt_channel2)
           obj.rpt_channel2.gui();
        end
    end
    function openRegAnalysis_CB(~,~)
        if ~isempty(obj.regAnalysis)
           obj.regAnalysis.gui();
        end
    end
    function reloadAllControls()
        if isempty(obj.workingDir) || ~exist(obj.workingDir,'dir')
            lastPath = pwd;
        else
            lastPath = obj.workingDir; %save this for file open dialogs after a reset
        end
        updateControls();
        setPairsTableData();
        if obj.initialized
            plotMainAxes(@obj.plotAllPairsWavelength, handles.axes);
        end
    end

    function resetAxes()
        cla(handles.axes,'reset');
        if obj.initialized
            handles.axes.Visible = 'on';
        else
            handles.axes.Visible = 'off';
        end
    end

    function updateControls()
        handles.file.saveFile.String = obj.saveFilePath;
        handles.file.spdata.String = obj.getFilePath('data');
        handles.file.rpt_channel1.String = obj.getFilePath('rpt_channel1');
        handles.file.rpt_channel2.String = obj.getFilePath('rpt_channel2');
        handles.file.regAnalysis.String = obj.getFilePath('regAnalysis');
        if obj.dirty
            handles.file.saveFile.BackgroundColor = obj.color_unsavedBG;
        else
            handles.file.saveFile.BackgroundColor = obj.color_editBG;
        end
    end

    function load_CB(~, ~)
        filepath = Pickle.selectExistingFileName(lastPath, [], obj.LoadableDataFormats, 'Load RegisteredPairsAnalysis');
        if ~isempty(filepath)
            if ~isempty(obj.workingDir) && exist(obj.workingDir,'dir')
                lastPath=obj.workingDir;
            end
            [~,~,ext] = fileparts(filepath);
            switch ext
                case SPData.saveFileExt
                    load_SPData(filepath);
                case obj.saveFileExt
                    obj.updateWaitbar(0,'Loading from .regpairs ...');
                    obj.load(filepath);
            
            end
            reloadAllControls();
            obj.updateWaitbar(1);
        end
    end

    function change_registration_CB(~,~)
        if obj.initialized
            path = obj.regAnalysis.workingDir;
        else
            path = lastPath;
        end
        filepath = Pickle.selectExistingFileName(path, RegistrationAnalysis.saveFileExt,...
                             RegistrationAnalysis.LoadableDataFormats, 'Load RegistrationAnalysis');
        if ~isempty(filepath)
            obj.updateWaitbar(0,'Re Loading with new .reganalysis ...');
            obj.load(obj.rpt_channel1, obj.rpt_channel2, filepath, obj.channelWavelengths);
            obj.updateWaitbar(0.9);
            reloadAllControls();
            obj.updateWaitbar(1);
        end
    end

    function load_SPData_CB(~,~)
        if obj.initialized
            path = obj.workingDir;
        else
            path = lastPath;
        end
        filepath = Pickle.selectExistingFileName(path,SPData.saveFileExt,...
                                                 SPData.SaveableDataFormats, 'Load SPData');
        load_SPData(filepath);
    end

    function load_SPData(spdata_file)
        % Helper that does the actual loading from SPData which includes checking to see if we need
        % a RegAnalysis
        if isempty(spdata_file)
            return
        end
        [path,~,~] = fileparts(spdata_file);
        if ~obj.initialized && ~isa(obj.regAnalysis, 'RegistrationAnalysis')
            regAnalysis = Pickle.selectExistingFileName(path, RegistrationAnalysis.saveFileExt,...
                             RegistrationAnalysis.LoadableDataFormats, 'Load RegistrationAnalysis');
        else
            regAnalysis = obj.regAnalysis;
        end
        if ~isempty(regAnalysis)
            obj.updateWaitbar(0,'Loading from SPData...');
            obj.load(spdata_file, regAnalysis, obj.channelWavelengths);
            obj.updateWaitbar(0.9);
            resetAxes();
            reloadAllControls();
            obj.updateWaitbar(1);
        end
    end


    function save_CB(~,~)
        obj.updateWaitbar(0,'Saving...');
        obj.save();
        obj.updateWaitbar(0.95);
        updateControls();
        obj.updateWaitbar(1);
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
        obj.guiFig = guiFig;
        resetAxes();
        reloadAllControls();
    end

    function batchProcess_CB(~,~)
        if obj.dirty
            error('RegisteredPairAnalysis:gui:batchProcess_CB','Object is dirty. Try saving first.');
        end
        title = 'Select SPData Files To Regsiter';
        if ~isempty(obj.workingDir) && exist(obj.workingDir,'dir')
            batchDir = fullfile(obj.workingDir,'..');
        else
            batchDir = pwd();
        end
        if isempty(obj.regAnalysis)
            reganalysis = Pickle.selectExistingFileName(batchDir,RegistrationAnalysis.saveFileExt,...
                                        RegistrationAnalysis.SaveableDataFormats,...
                                        'Select RegistrationAnalysis to use in registration');
        else
            reganalysis = obj.regAnalysis;
        end
        filepaths = Pickle.selectBatchProccesingFileNames(batchDir,'*.spdata',SPData.SaveableDataFormats,title);        
        if isempty(filepaths); return; end
        overwriteFlag = 1; % ask for overwrite        
        rpairs_files = RegisteredPairAnalysis.batchProcess([], filepaths, reganalysis, [], overwriteFlag);
        if isempty(rpairs_files); return; end
        count = sum(cellfun(@ischar, rpairs_files));
        msgbox(sprintf('Succesfully processed %i/%i .regpairs files.',count,numel(rpairs_files)));
    end

    function close_CB(~,~)
        try
            handles.pairsTable.MousePressedCallback=[];
            GUIBuilder.clearJavaContextMenu(handles.contextmenus.pairs);
        catch err
            disp(getReport(err));
        end
        obj.closeGUI();
    end
end
