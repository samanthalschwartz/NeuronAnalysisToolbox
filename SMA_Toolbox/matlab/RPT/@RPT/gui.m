% RPT.gui

function guiFig = gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end

    %% Constants and Global Initialization
    gui@GUIBuilder(obj); %Call GUIBuilder initialization
    gui_name = '[RPT] Robust Particle Tracking';
    lastPath = obj.workingDir; %save this for file open dialogs after a reset
    
    uH = GUIBuilder.default_unitHeight; % unit height for elements
    boarder = GUIBuilder.default_boarder;%Boarder width around the outside of the gui
    sp = GUIBuilder.default_spacing; %spacing between elements.
    but_sz = GUIBuilder.default_buttonSize; %Button size
    
    fig_sz = [1080 600]; %figure size
    halfw = 530;
    reserveH = 180; %Hight above the tracking panel to reserve
    fullw_sz = [fig_sz(1)-2*boarder-2*sp uH]; %size of a full-width component that stays within the boarders.
    halfw_sz = [halfw-sp-boarder uH];
    axes_pos = [halfw boarder fig_sz(1)-halfw fig_sz(2)-2*boarder]; % To be updated later
    tabs_pos = [sp, sp, halfw_sz(1), fig_sz(2)-boarder]; %initial size will be fixed by alignPanels
    fullPan_pos = [boarder boarder fullw_sz];
    %These color the current phase button
    PHASE_BG_COLORS = [.5 0 0; 1 0.5 0.5; 1 0.5 0.2; 1 .8 0.4; 1 1 0.4; 0.8 1 0.2; 0.2 1 0.2]; 


    %% Create figure and controls
    % Make figure
    guiFig = figure('Units','pixels','Position',[10 0 fig_sz],'Resize','off',...
                    'MenuBar','none','ToolBar','none','NumberTitle','off',...
                    'Name',gui_name,'Visible','on',...
                    'CloseRequestFcn',@close_CB);
    obj.guiFig = guiFig;
    whitebg(guiFig,[0.1, 0.1, 0.1]);

    % Make Panels
    handles.panels.file = uipanel('Parent',guiFig,'Units','Pixels','Position',fullPan_pos,'Title','File Associations');
    handles.tabG = uitabgroup('Parent',guiFig,'Units','Pixels','Position',tabs_pos,'TabLocation','top');
    handles.tabs.tracking = uitab('Parent',handles.tabG,'Title','Tracking');
    handles.tabs.stats = uitab('Parent',handles.tabG,'Title','Stats');
    populateFilePanel(handles.panels.file);    
    populateStatsPanel(handles.tabs.stats);    
    populateTrackingPanel(handles.tabs.tracking);
    populateFindMaximaTab(handles.tabs.findMaxima);
    populateFilterMaximaTab(handles.tabs.filterMaxima);
    populateLocalizeEmittersTab(handles.tabs.localizeEmitters);
    populateFilterEmittersTab(handles.tabs.filterEmitters);
    populateTrackEmittersTab(handles.tabs.trackEmitters);

    %% initialization
    alignPanels();
    createMenus();
    initializeAxes(); %now that we know where the panels are we can assign the axes to the remaining space
    updateControls();

    function populateFilePanel(panH)
        % File Panel
        % Bottom row of buttons
        bot_row.pos = [sp, sp, fullw_sz];
        bot_row.names = {'Load','Save','SaveAs', 'SetDefaultProps', 'BatchProcess', 'OpenSPData'};
        bot_row.CBs = {@load_CB,@save_CB,@saveAs_CB, @setDefaultProps_CB, @batchProcess_CB, @openSPData_CB};
        GUIBuilder.buttonRow(panH, bot_row.pos, but_sz+[20 0], bot_row.names, bot_row.CBs);
        %Quit is out to the side so we don't include it.
        handles.quitButton = uicontrol('Parent',panH,'Style','pushbutton','String','Quit',...
            'Position',[fullw_sz(1)-but_sz(1) sp but_sz],...
            'Callback',@close_CB);
        edits.pos = bot_row.pos + [0, but_sz(2)+sp, 0, 2*sp];
        edits.hNames = {'spdata','saveFile'}; %names in the handles structure returned
        edits.labels = {'SPData (.spdata) Path', 'Save File (.rpt) Path'};
        edits.values = {obj.getFilePath('data'), obj.saveFilePath};
        edits.CBs = {0,0};
        handles.file = GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
    end

    function populateStatsPanel(panH)
        col_names = {'Property', 'Value'};
        pos = [sp sp halfw_sz(1)-2*boarder sp];
        col_width = {150, pos(3)-150-sp};
        handles.StatsTable = uitable(panH,'ColumnName',col_names,'Position', pos,...
                                          'ColumnWidth',col_width, 'RowName',[]);       
    end
   
    function populateTrackingPanel(panH)
        % Bottom row of buttons
        TGpos = [sp,sp, halfw_sz(1)-2*sp, fig_sz(2)-reserveH];
        buttons.pos = [sp, TGpos(2)+TGpos(4)+sp, halfw_sz(1), uH];
        buttons.names = {'Auto Track','Run Next'};
        buttons.CBs = {@autoTrack_CB,@runNextPhase_CB};
        butH = GUIBuilder.buttonRow(panH, buttons.pos, but_sz, buttons.names, buttons.CBs);
        but_pos = get(butH{end},'Position');
        butEnd_h = but_pos(1)+but_pos(3)+sp;
        phaselabel_pos = [butEnd_h, buttons.pos(2), 120, uH*0.75];
        phase_pos = [phaselabel_pos(1)+phaselabel_pos(3)+sp, buttons.pos(2), halfw_sz(1)-(phaselabel_pos(1)+phaselabel_pos(3)+4*sp), uH];
        uicontrol('Parent',panH,'Style','text', 'String','Current Phase Status:','Position',phaselabel_pos,'HorizontalAlignment','right');
        handles.phaseDesc = uicontrol('Parent',panH,'Style','edit', 'String',obj.phase,'Position',phase_pos, 'HorizontalAlignment','center');
        handles.trackingTG = uitabgroup('Parent',panH,'Units','Pixels','Position',TGpos);
        handles.trackingTG.TabLocation = 'left';
        handles.tabs.findMaxima = uitab('Parent',handles.trackingTG,'Title', makePhaseTitle(3),'Units','Pixels','UserData',3);    
        handles.tabs.filterMaxima = uitab('Parent',handles.trackingTG,'Title', makePhaseTitle(4),'Units','Pixels','UserData',4);    
        handles.tabs.localizeEmitters = uitab('Parent',handles.trackingTG,'Title', makePhaseTitle(5),'Units','Pixels','UserData',5);    
        handles.tabs.filterEmitters = uitab('Parent',handles.trackingTG,'Title', makePhaseTitle(6),'Units','Pixels','UserData',6);    
        handles.tabs.trackEmitters = uitab('Parent',handles.trackingTG,'Title', makePhaseTitle(7),'Units','Pixels','UserData',7);    
        handles.result_buttons = cell(1,obj.N_PHASES);
        handles.inspect_buttons = cell(1,obj.N_PHASES);
        handles.params_table = cell(1,obj.N_PHASES);        
    end

    function populateFindMaximaTab(tabH)
        %Configure the Analysis tab for the Find Maxima phase
        phaseIdx=3;
        
        plot_buts.names = {'FSumImage','RMaxWeight','RMax/Frame','RMax/Scale'};
        plot_buts.CBs = {@plotFilteredSumImage_CB,@plotRawMaximaImage_CB, @plotRawMaximaPerFrame_CB, @plotRawMaximaPerScale_CB};
        
        view_buts.names = {'FiltFrames','RMaxima','RMaximaFilt'};
        view_buts.CBs = {@viewFilteredFrames_CB, @(~,~) obj.viewRawMaximaMovie(),...
                         @(~,~) obj.viewRawMaximaFilteredMovie()};
        
        inspect_buts.names = {'SumImage','View Frames'};
        inspect_buts.CBs = {@plotSumImage_CB, @viewFrames_CB};
        
        params.method.range = {'LoG','DoG'};        
        params.method.disp = 'Filter Method';
        params.method.desc = 'Filter Shape (LoG=Laplacian of Gaussian) (DoG=Difference of Gaussian)';
        params.filterSigmas.disp = 'Filter Kernel Sigma';
        params.filterSigmas.desc = 'Should match PSF of emitter.  Format:[sigmaX sigmaY]';
        params.filterSigmas.range = [];
        params.maximaNeighborhoodSize.range = {int32(3), int32(5),int32(7),int32(9),int32(11),int32(13),int32(15)};
        params.maximaNeighborhoodSize.disp = 'Maxima Neighborhood Size';
        params.maximaNeighborhoodSize.desc = 'Width of neighborhood around each maxima.';
        params.scaleNeighborhoodSize.range = {int32(3), int32(5),int32(7),int32(9),int32(11),int32(13),int32(15)};
        params.scaleNeighborhoodSize.disp = 'Scale Space Neighborhood Size';
        params.scaleNeighborhoodSize.desc = 'Width of neighborhood around each maxima in scale space.';
        
        populatePhaseTab(tabH,params,plot_buts,view_buts,inspect_buts,phaseIdx);
    end

    function populateFilterMaximaTab(tabH)
        %Configure the Analysis tab for the Filter Maxima phase
        phaseIdx=4;
        
        plot_buts.names = {'MaximaWeight','EmitterSumImage', 'Max/Frame', 'Max/Scale', 'BoxesImage'};
        plot_buts.CBs = {@plotMaximaImage_CB,@plotEmitterSumImage_CB, @plotMaximaPerFrame_CB,...
                         @plotMaximaPerScale_CB, @plotBoxesImage_CB};
        
        view_buts.names = {'Maxima', 'FilterMaxima', 'BoxesSlices', 'BoxesMovie',...
                           'BoxesFiltMovie','EmitterImages', 'EmitterMovie'};
        view_buts.CBs = {@(~,~) obj.viewMaximaMovie(),...
                         @(~,~) obj.viewMaximaFilteredMovie(),...
                         @(~,~) obj.viewBoxesSlices(),...
                         @(~,~) obj.viewBoxesMovie(),...
                         @(~,~) obj.viewBoxesFilteredMovie(),...
                         @(~,~) obj.viewEmitterImages(),...
                         @(~,~) obj.viewEmitterMovie()};
        
        inspect_buts.names = {'Threshold'};
        inspect_buts.CBs = {@plotThreshold_CB};
        
        params.maximaThreshold.disp = 'Maxima Threshold';
        params.maximaThreshold.desc = 'The minimum filtered pixel value for a maxima to be considered an emitter.  Set to -1 for auto threshold.';
        params.minimumBoxSize.disp = 'Minimum Box Size';
        params.minimumBoxSize.desc = 'The Minimum Box Size [X Y] that the Boxxer will generate.';
        params.optimalBoxSize.disp = 'Optimal Box Size';
        params.optimalBoxSize.desc = 'The Preferred Box Size [X Y] for the Boxxer to generate.';       
        
        populatePhaseTab(tabH,params,plot_buts,view_buts,inspect_buts,phaseIdx);
    end

    function populateLocalizeEmittersTab(tabH)
        %Configure the Analysis tab for the Localize Emitters phase
        phaseIdx=5;
        
        plot_buts.names = {'RawPosDist','Raw I Dist', 'RawSigma Dist', 'Raw LLH dist'};
        plot_buts.CBs = {@plotRawLocalizationPosDist_CB, @plotRawLocalizationIDist_CB,...
                         @plotRawLocalizationSigmaDist_CB,@plotRawLocalizationLLHDist_CB};
        
        view_buts.names = {'Raw SRImage', 'EmitterFits', 'Raw Loc Movie'};
        view_buts.CBs = {@(~,~) obj.viewRawEmitterSuperResGauss(),...
                         @(~,~) obj.viewEmitterModelComparison(),... 
                         @(~,~) obj.viewRawLocalizationMovie()};
        
        inspect_buts.names = {};
        inspect_buts.CBs = {};
        
        params.model.disp = 'Model';
        params.model.desc = 'The MAPPEL Class name for the approprtiate emitter parameterization model.';
        params.estimator.disp = 'Estimator';
        params.estimator.desc = 'The Estimation method for finding the maximum likelyhood parameters.';
        
        populatePhaseTab(tabH,params,plot_buts,view_buts,inspect_buts,phaseIdx);
    end

    function populateFilterEmittersTab(tabH)
        %Configure the Analysis tab for the Filter Emitters phase
        phaseIdx=6;
        
        plot_buts.names = {'LocSumImage','SRImage','PosDist','I Dist', 'Sigma Dist', 'LLH dist'};
        plot_buts.CBs = {@plotLocalizationSumImage_CB, @plotEmitterSuperResGauss_CB, @plotLocalizationPosDist_CB, @plotLocalizationIDist_CB,...
                         @plotLocalizationSigmaDist_CB,@plotLocalizationLLHDist_CB};
        
        view_buts.names = {'SRView', 'Selected BoxMovie','Loc Movie'};
        view_buts.CBs = {@(~,~) obj.viewEmitterSuperResGauss(),...
                         @(~,~) obj.viewSelectedBoxesMovie(),...
                         @(~,~) obj.viewLocalizationMovie()};
        
        inspect_buts.names = {};
        inspect_buts.CBs = {};
       
        params.minIntensity.disp = 'Min Intensity';
        params.minIntensity.desc = 'Minimum Locazliation Fit Intensity for acceptance';
        params.minSigma.disp = 'Min Sigma';
        params.minSigma.desc = 'Minumum Apparent Gaussian Sigma';
        params.maxSigma.disp = 'Max Sigma';
        params.maxSigma.desc = 'Maximum Apparent Gaussian Sigma';
        params.maxPositionSE.disp = 'Max Position SE';
        params.maxPositionSE.desc = 'Maximum standard error in positions (x or y). [pixels]';
        params.certVsUniformModel.disp = 'Certainty Vs Uniform Background Model';
        params.certVsUniformModel.desc = 'Probability that emitter is explained by simpler 1-parameter uniform background model.  [Set to -1 to disable.]';
        params.certVsNoiseModel.disp = 'Certainty Vs Noise Model';
        params.certVsNoiseModel.desc = 'Probability that emitter is explained by more complex noise model.  [Set to -1 to disable.]';
        params.overlapDistance.disp = 'Min Localization Overlap Distance';
        params.overlapDistance.desc = 'Minimum distance 2 localizations in the same frame can be in pixels.';
        populatePhaseTab(tabH,params,plot_buts,view_buts,inspect_buts,phaseIdx);
    end

    function populateTrackEmittersTab(tabH)
        %Configure the Analysis tab for the Track Emitters phase
        phaseIdx=7;
        
        plot_buts.names = {'Tracks3D'};% 'Track SumImage', 'Track Len Dist'};
        plot_buts.CBs = {@plot3DTrackSequence_CB}; %@plotTrackSumImage_CB, @plotTrackLengthDist_CB};
        
        view_buts.names = {'Tracks3DSeq', 'Tracks3DSpeed', 'Tracks3DTemporal', 'TracksMovie', 'SegmentAnalysis'};
        view_buts.CBs = {@(~,~) obj.view3DTrackSequence(),...
                         @(~,~) obj.view3DTrackSpeed(),...
                         @(~,~) obj.view3DTrackTemporal(),...
                         @(~,~) obj.viewTrackMovie(),...
                         @makeTrackSegmentAnalysis_CB};
        
        inspect_buts.names = {'F2F Debug'};
        inspect_buts.CBs = {@debugFrame2Frame_CB};
       
        params.D.disp = 'Diffusion Constant';
        params.D.desc = '[Units: px^2/frame] The nominal diffusion constant for all tracks';
        params.Kon.disp = 'Kon';
        params.Kon.desc = '[Units: 1/frame] Rate of new track appearence';
        params.Koff.disp = 'Koff';
        params.Koff.desc = '[Units: 1/frame] Rate of track ending';
        params.MaxSpeed.disp = 'Max Speed';
        params.MaxSpeed.desc = '[Units: px/frame] The Maximum Velocity for Track Connections.';
        params.MaxGapCloseFrames.disp = 'Max Gap Close Frames';
        params.MaxGapCloseFrames.desc = 'Maximum number of frames to connect tracks over in gap-closing.';
        params.MinGapCloseTrackLength.disp = 'Min Gap Close Track Length';
        params.MinGapCloseTrackLength.desc = 'The minimum size of track fragments to be used in gap-closing';
        params.MinFinalTrackLength.disp = 'Min Final Track Length';
        params.MinFinalTrackLength.desc = 'The minimum final length of tracks to be returned.';
        populatePhaseTab(tabH,params,plot_buts,view_buts,inspect_buts,phaseIdx);
    end

    function updatePhaseParamPane(phaseIdx)
        %Update the values in the PropertyGrid object using the phase
        %parameters for this phaseIdx.  This is run after a compute phase
        %where the parameters may have changed.
        for phase = phaseIdx
            H = handles.params_table{phase};
            if isempty(H); return; end
            pgf = H.Properties;
            params = obj.getParams(phase);
            fns = fieldnames(params);
            if obj.phaseIdx>=2;
                for i = 1:length(fns)
                    fn = fns{i};
                    pgf.FindByName(fn).Value = params.(fn);
                end
            end
            H.Properties = pgf;
        end
    end

    function updatePhaseTabs()
        for i = 1:obj.N_PHASES
            setEnable(handles.result_buttons{i},obj.phaseIdx>=i);
            setEnable(handles.inspect_buttons{i},obj.phaseIdx>=i-1 && obj.phaseIdx>=2);
            GUIBuilder.setPropertyGridEnabled(handles.params_table{i},obj.phaseIdx<i && obj.phaseIdx>=2);
            if i==obj.phaseIdx && i>2
                handles.trackingTG.SelectedTab = handles.trackingTG.Children(i-2);
            end
        end
        if obj.dirty
            handles.file.saveFile.BackgroundColor = obj.color_unsavedBG;
        else
            handles.file.saveFile.BackgroundColor = obj.color_editBG;
        end
    end

    function setEnable(Hs,state)
        if ~iscell(Hs)
            Hs = {Hs};
        end
        for h = 1:length(Hs)
            if ischar(state)
                Hs{h}.Enable = state;
            elseif state
                Hs{h}.Enable = 'on';
            else
                Hs{h}.Enable = 'off';
            end
        end
    end
    
    function populatePhaseTab(tabH,param_info,plot_buttons, view_buttons, inspect_buttons, phaseIdx)
        phase_title = obj.PHASE_TITLES{phaseIdx};
        pos = tabH.Position;
        w = pos(3)-pos(1)-130;
        pan_pos = [sp, sp, w, sp];
        plotPanH = uipanel('Parent',tabH,'Units','Pixels','Position',pan_pos,'Title','Result Plots'); 
        viewPanH = uipanel('Parent',tabH,'Units','Pixels','Position',pan_pos,'Title','Result Views');
        butHplot = GUIBuilder.buttonRow(plotPanH, pan_pos, but_sz+[20 0], plot_buttons.names, plot_buttons.CBs);
        butHview = GUIBuilder.buttonRow(viewPanH, pan_pos, but_sz+[20 0], view_buttons.names, view_buttons.CBs);
        resultButH = [butHplot butHview];
        GUIBuilder.autoSizePanel(plotPanH);
        plotPanH.Position(3) = w;
        GUIBuilder.autoSizePanel(viewPanH);
        viewPanH.Position(3) = w;
        panels = {viewPanH, plotPanH};
        if ~isempty(inspect_buttons.names)
            inspectPanH = uipanel('Parent',tabH,'Units','Pixels','Position',pan_pos,'Title','Debug/Inspect Plots&Views');  
            inspectButH = GUIBuilder.buttonRow(inspectPanH, pan_pos, but_sz+[20 0], inspect_buttons.names, inspect_buttons.CBs);
            GUIBuilder.autoSizePanel(inspectPanH);
            inspectPanH.Position(3) = w;
            panels = [panels {inspectPanH}];
        else
            inspectButH = [];
        end
        GUIBuilder.align(panels,'Left','Fixed',sp);
        panHeight = sum(cellfun(@(pH) 2*sp+pH.Position(4), panels));
        param_pos = [sp, panHeight, w, pos(4)-2*sp-panHeight];
        pane_title = sprintf('%s Properties',phase_title);
        tableH = GUIBuilder.makePropertyGrid(tabH,pane_title,obj.getParams(phaseIdx),param_info, param_pos);
        GUIBuilder.setPropertyGridMenu(tableH,{@paramsTableMouse_CB, phaseIdx});
        handles.result_buttons{phaseIdx} = resultButH;
        handles.inspect_buttons{phaseIdx} = inspectButH;
        handles.params_table{phaseIdx} = tableH;        
    end

    function createMenus()
        %File Menu
        %File Menu
        labels = {'Load ...','Save','Save As ...','Reset Object',[],...
                  'Set Default Properties ...', 'Batch Process...','Open SPData...',[],...
                  'Quit'};
        CBs = {@load_CB, @save_CB, @saveAs_CB, @resetObject_CB,[],...
               @setDefaultProps_CB,@batchProcess_CB,@openSPData_CB,[],...
               @close_CB};
        GUIBuilder.makeFigureMenu(guiFig,'File',labels,CBs);
        
        %This menu will be reused for all Track objects
        labels = {'Auto Track', 'Run Next Phase', 'Reset Tracking'};
        CBs = {@autoTrack_CB,  @runNextPhase_CB, @resetTracking_CB};
        GUIBuilder.makeFigureMenu(guiFig,'Track',labels,CBs);
        
        %This menu will be reused for all Track objects
        labels = {'Export Localizations...', 'Export Tracks...', 'Track Segment Analysis'};
        CBs = {@exportLocalizations_CB,  @exportTracks_CB, @makeTrackSegmentAnalysis_CB};
        handles.menus.export = GUIBuilder.makeFigureMenu(guiFig,'Export',labels,CBs);
        if obj.phaseIdx==7
            handles.menus.export.Enable = 'on';
        else
            handles.menus.export.Enable = 'off';
        end
        labels = {'Run'};
        CBs = {@runPhase_CB};
        handles.contextmenus.phase_incomplete = GUIBuilder.makeContextMenu(labels, CBs);
        handles.contextmenus.phase_incomplete.Parent = guiFig; %Sometimes this is set to gcf which might not be the correct figure.
        handles.contextmenus_java.phase_incomplete = GUIBuilder.makeJavaContextMenu(labels, CBs);
        
        
        labels = {'ResetPhase', 'RerunPhase'};
        CBs = {@resetPhase_CB, @rerunPhase_CB};
        handles.contextmenus.phase_complete = GUIBuilder.makeContextMenu(labels, CBs);
        handles.contextmenus.phase_complete.Parent = guiFig;
        handles.contextmenus_java.phase_complete = GUIBuilder.makeJavaContextMenu(labels, CBs);
    end

    function panels_pos = alignPanels()
        pans = fieldnames(handles.panels);
        for n = 1:length(pans)
            GUIBuilder.autoSizePanel(handles.panels.(pans{n}));
        end
        panels_pos = GUIBuilder.align(struct2cell(handles.panels),'Left','Fixed',sp);
        panH = panels_pos(2)+panels_pos(4)+sp;
        tabs_pos = [sp, panH, halfw_sz(1), fig_sz(2)-panH-boarder];
        handles.tabG.Position = tabs_pos;
        handles.StatsTable.Position = [sp, sp, tabs_pos(3)-3*sp, tabs_pos(4)-3*sp-uH];
    end

    function initializeAxes()
        %compute the global variable giving the axes position
        % Also set up axes manipulation buttons
        % Make image axes
        handles.axes = axes('Parent', guiFig, 'Units', 'Pixels','Visible','off');

        tk_pos = round(handles.tabG.Position);
        left_pos = tk_pos(1)+tk_pos(3);
        axes_control_pos = [left_pos, tk_pos(2), fig_sz(1)-left_pos-boarder, uH];
        axes_pos = [left_pos, tk_pos(2)+uH+sp, fig_sz(1)-left_pos-boarder, fig_sz(2)-boarder-tk_pos(2)-uH-sp];
        handles.axes.Position = axes_pos;
        c_pos = axes_control_pos;
        c_pos(3) = 80;
        b_pos = [c_pos(1:2) 28 24];
        colormap('gray');

        %These are the buttons to change the plot mode that emulate the toolbar
        H{1} = uicontrol('Parent',guiFig,'Style','text', 'String','Colormap:',...
                         'Position',c_pos,'HorizontalAlignment','right','BackgroundColor',guiFig.Color);
        H{2} = uicontrol('Parent',guiFig,'Style','popupmenu', 'String',GUIBuilder.colormaps,...
                         'Position',c_pos,'HorizontalAlignment','center','Callback',@changeColormap_CB);
        H{3} = uicontrol('Parent',guiFig,'Style','togglebutton','CData',GUIBuilder.readMatlabIcon('tool_rotate_3d.png'),...
                         'Position',b_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('rotate3d',a,b));
        H{4} = uicontrol('Parent',guiFig,'Style','togglebutton','CData',GUIBuilder.readMatlabIcon('tool_zoom_in.png'),...
                         'Position',b_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('zoom',a,b));
        H{5} = uicontrol('Parent',guiFig,'Style','togglebutton','CData',GUIBuilder.readMatlabIcon('tool_hand.png'),...
                         'Position',b_pos,'HorizontalAlignment','center','Callback',@(a,b) axesMode_CB('pan',a,b));
        GUIBuilder.align(H,'Fixed',sp,'Bottom');
        H{1}.Position = H{1}.Position - [0 0 0 5];
        handles.axes_colormap=H{2};
        handles.axes_buttons.rotate3d = H{3};
        handles.axes_buttons.zoom = H{4};
        handles.axes_buttons.pan = H{5};
    end

    function populateStatsTableData()
        stats = obj.getStats();
        handles.StatsTable.Data = makeParamTableData(stats);
    end

    function data = makeParamTableData(st)
        % [IN] st - a Structure
        data = [fieldnames(st) cellmap(@arr2str, struct2cell(st))];
    end

    function reloadAllControls()
        % This is a hard update.  Everything is redone
        for i = 3:7
            updatePhaseParamPane(i);
        end
        updateFileAssociationControls();
        updateControls();
    end

    function updateControls()
        % This is a soft update
        populateStatsTableData();
        updatePhase()
    end

    function title = makePhaseTitle(idx)
        %Color the phase title appropriately
        if obj.phaseIdx<idx
            color = '#FF8888';
        else
            color = '#88FF88';
        end
        title = sprintf('<html><div style = "background:%s;">%s</div></html>',color,obj.PHASE_TITLES{idx});
    end

    function ctxH = getPhaseUIContextMenu(idx)
        %Color the phase title appropriately
        if obj.phaseIdx<idx
            ctxH = handles.contextmenus.phase_incomplete;
        else
            ctxH = handles.contextmenus.phase_complete;
        end
    end


    function paramsTableMouse_CB(H,eData,idx)
        %This makes the context menus work for the properties tables
        if eData.isMetaDown
            if obj.phaseIdx<idx
                menuH = handles.contextmenus_java.phase_incomplete;
            else
                menuH = handles.contextmenus_java.phase_complete;
            end
            menuH.show(H,eData.getX(),eData.getY());
            menuH.repaint();
        end
    end

    function updatePhase()
        %Set colored titles and assign the correct context menu
        handles.tabs.findMaxima.Title = makePhaseTitle(3);
        handles.tabs.findMaxima.UIContextMenu = getPhaseUIContextMenu(3);
        handles.tabs.filterMaxima.Title = makePhaseTitle(4);
        handles.tabs.filterMaxima.UIContextMenu = getPhaseUIContextMenu(4);
        handles.tabs.localizeEmitters.Title = makePhaseTitle(5);
        handles.tabs.localizeEmitters.UIContextMenu = getPhaseUIContextMenu(5);
        handles.tabs.filterEmitters.Title = makePhaseTitle(6);
        handles.tabs.filterEmitters.UIContextMenu = getPhaseUIContextMenu(6);
        handles.tabs.trackEmitters.Title = makePhaseTitle(7);
        handles.tabs.trackEmitters.UIContextMenu = getPhaseUIContextMenu(7);

       
        populateStatsTableData();
        updatePhaseTabs();
        plotPhase(); %updatre plot to show most relevent phase
        handles.phaseDesc.String = obj.phase;
        handles.phaseDesc.BackgroundColor = PHASE_BG_COLORS(obj.phaseIdx,:);
        if obj.phaseIdx==7
            handles.menus.export.Enable = 'on';
        else
            handles.menus.export.Enable = 'off';
        end

        drawnow();
    end

    function updateFileAssociationControls()
        %Update the file panel when needed
        handles.file.spdata.String = obj.getFilePath('data');
        handles.file.saveFile.String = obj.saveFilePath;
        if obj.dirty
            handles.file.saveFile.BackgroundColor = obj.color_unsavedBG;
        else
            handles.file.saveFile.BackgroundColor = obj.color_editBG;
        end
    end

    function load_CB(~, ~)
        if isempty(lastPath)
            lastPath = obj.workingDir;
        end
        [filename, pathname] = uigetfile(obj.LoadableDataFormats,'Select data source', lastPath);
        if ~filename; return; end
        obj.load(fullfile(pathname,filename));
        obj.inGui=true;
        lastPath = obj.workingDir;
        reloadAllControls();
    end

    function save_CB(~,~)
        obj.save();
        updateFileAssociationControls();
        updateControls();
    end

    function saveAs_CB(~,~)       
        obj.saveas();
        updateControls();
    end

    function resetObject_CB(~,~)
        if obj.initialized
            lastPath = obj.workingDir; %Only change the lastPath if we actually have a valid workingDir
        end
        obj.resetObject();
        cla(handles.axes);
        handles.axes.Visible='off';
        reloadAllControls();
    end
   
    function setDefaultProps_CB(~,~)
        openDir = obj.workingDir;
        if ~exist(openDir, 'dir')
            openDir = obj.data.workingDir;
        end
        defaultFile = Pickle.selectExistingFileName(openDir,obj.saveFileExt,obj.SaveableDataFormats,'Select file to copy parameters from');
        if isempty(defaultFile); return; end
        if strcmp(defaultFile, obj.saveFilePath)
            error('RPT:setDefaultProps_CB','Must choose different file for a default parameters file');
        end
        obj.setDefaultParams(RPT(defaultFile));
        reloadAllControls();
    end

    function autoTrackInteractive(endPhase)
        try
            nPhases = length(obj.PHASES);
            if nargin<1
                endPhase = nPhases;
            end
            if obj.phaseIdx>=nPhases
               error('RPT:autoTrack','Desired phase %i > Maximum phase %i',endPhase,nPhases);
            end
            if obj.phaseIdx==1
                error('RPT:autoTrack','Object is invalid (phaseIdx=1)');
            end
            if obj.phaseIdx<3 && endPhase>=3
                obj.setParams(3, handles.params_table{3}.GetPropertyValues());
                obj.findMaxima();
                updatePhaseParamPane(obj.phaseIdx);
                updatePhase();                
            end
            if obj.phaseIdx<4 && endPhase>=4
                obj.setParams(4, handles.params_table{4}.GetPropertyValues());
                obj.filterMaxima();
                updatePhaseParamPane(obj.phaseIdx);
                updatePhase();
            end
            if obj.phaseIdx<5 && endPhase>=5
                obj.setParams(5, handles.params_table{5}.GetPropertyValues());
                obj.localizeEmitters();
                updatePhaseParamPane(obj.phaseIdx);
                updatePhase();
            end
            if obj.phaseIdx<6 && endPhase>=6
                obj.setParams(6, handles.params_table{6}.GetPropertyValues());
                obj.filterEmitters();
                updatePhaseParamPane(obj.phaseIdx);
                updatePhase();
            end
            if obj.phaseIdx<7 && endPhase==7
                obj.setParams(7, handles.params_table{7}.GetPropertyValues());
                obj.trackEmitters();
                updatePhaseParamPane(obj.phaseIdx);
                updatePhase();   
            end
            if obj.phaseIdx~=endPhase;
                error('RPT:autoTrack','Could not complete tracking to desired endPhase:%s.  Current phase: %s',obj.PHASES{endPhase},obj.phase);
            end
        catch err
            updatePhase();   
            disp(getReport(err));
            errordlg(err.message,err.identifier);
        end    
    end

    function plotPhase()
        % The default plot for a particular phase
        switch obj.phaseIdx
            case 2; plotSumImage_CB(0,0);
            case 3; plotRawMaximaImage_CB(0,0);
            case 4; plotMaximaImage_CB(0,0);
            case 5; plotSumImage_CB(0,0);
            case 6; plotLocalizationSumImage_CB(0,0);
            case 7; plot3DTrackSequence_CB(0,0);
        end
    end

    function autoTrack_CB(~,~)
        autoTrackInteractive();
    end

    function runNextPhase_CB(~,~)
        autoTrackInteractive(obj.phaseIdx+1);
    end   
    
    function tabPhase=getCurrentTabPhase()
        tabPhase = handles.trackingTG.SelectedTab.UserData;
    end

    function runPhase_CB(~,~)
        autoTrackInteractive(getCurrentTabPhase());
    end

    function resetPhase_CB(~,~)
        tabPhase = getCurrentTabPhase();
        try
            obj.setPhase(tabPhase-1);
        catch err
            updatePhase();   
            disp(getReport(err));
            errordlg(err.message,err.identifier);
        end    
        updatePhase();
        handles.trackingTG.SelectedTab = handles.trackingTG.Children(tabPhase-2);
    end

    function rerunPhase_CB(~,~)
        tabPhase = getCurrentTabPhase();
        try
            obj.setPhase(tabPhase-1);
        catch err
            updatePhase();   
            disp(getReport(err));
            errordlg(err.message,err.identifier);
        end    
        updatePhase();
        autoTrackInteractive(tabPhase)
    end

    function resetTracking_CB(~,~)
        obj.setPhase(2); %initialized
        updatePhase();
    end

    function close_CB(~,~)
        % Most of these items must be cleared or there will be a memory leak since we are passing
        % callbacks to Java.
        try
            GUIBuilder.clearPropertyGridMenus(handles.params_table(3:7));
            GUIBuilder.clearJavaContextMenu(handles.contextmenus_java.phase_incomplete);
            GUIBuilder.clearJavaContextMenu(handles.contextmenus_java.phase_complete);
        catch err
            disp(getReport(err));
        end
        obj.closeGUI();
    end

    function batchProcess_CB(~,~)
        if obj.dirty
            error('SPData:gui:batchProcess_CB','Object is dirty. Try saving first.');
        end
        title = sprintf('Select SPData Files To Track ROI Name: "%s" for:',obj.ROIname);
        filepaths = Pickle.selectBatchProccesingFileNames(obj.data.workingDir,'*.spdata',SPData.SaveableDataFormats,title);
        if isempty(filepaths); return; end
        overwriteFlag = 1; % ask for overwrite
        rpt_files = RPT.batchProcess([], filepaths, obj.ROIname, obj, overwriteFlag);
        if isempty(rpt_files); return; end
        count = sum(cellfun(@ischar, rpt_files));
        msgbox(sprintf('Succesfully processed %i/%i .rpt files for ROI: "%s" ',count,numel(rpt_files), obj.ROIname));
    end

    function openSPData_CB(~,~)
        obj.data.gui();
    end

    %% Viewing
    function changeColormap_CB(hObj,~)
        colormap(hObj.String{hObj.Value});
    end

    function axesMode_CB(mode,~,~)
        modes = {'rotate3d','pan','zoom'};
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

    function viewFrames_CB(~,~)
        obj.viewMaximizedDipFig(dip_image(obj.getFrames()),'Name','Calibrated Frames');
    end

    function viewFilteredFrames_CB(~,~)
        name = sprintf('Filtered Frames sigma = %s',arr2str(obj.ParamsFindMaxima.filterSigmas));
        obj.viewMaximizedDipFig(dip_image(obj.getFilteredFrames()),'Name',name);
    end

    function plotMainAxes(plotFun, varargin)
        %plot something into the main axes
        axes(handles.axes);
        cla(handles.axes,'reset');
        plotFun(varargin{:});
        cmapObj = handles.axes_colormap;
        colormap(cmapObj.String{cmapObj.Value});
        axesMode_CB('rotate3d',[],[]);
    end

    function plotSumImage_CB(~,~)
        plotMainAxes(@obj.plotSumImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    function plotFilteredSumImage_CB(~,~)
        plotMainAxes(@obj.plotFilteredSumImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    function plotEmitterSumImage_CB(~,~)
        plotMainAxes(@obj.plotEmitterSumImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    function plotRawMaximaImage_CB(~,~)
        plotMainAxes(@obj.plotRawMaximaImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    function plotMaximaImage_CB(~,~)
        plotMainAxes(@obj.plotMaximaImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    function plotRawMaximaPerFrame_CB(~,~)
        plotMainAxes(@obj.plotRawMaximaPerFrame,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotRawMaximaPerScale_CB(~,~)
        plotMainAxes(@obj.plotRawMaximaPerScale,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotMaximaPerFrame_CB(~,~)
        plotMainAxes(@obj.plotMaximaPerFrame,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotMaximaPerScale_CB(~,~)
        plotMainAxes(@obj.plotMaximaPerScale,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotBoxesImage_CB(~,~)
        plotMainAxes(@obj.plotBoxesImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    %% Phase 5 Plots and Views

    function plotRawLocalizationPosDist_CB(~,~)
        plotMainAxes(@obj.plotRawLocalizationPosDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotRawLocalizationIDist_CB(~,~)
        plotMainAxes(@obj.plotRawLocalizationIDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotRawLocalizationSigmaDist_CB(~,~)
        plotMainAxes(@obj.plotRawLocalizationSigmaDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotRawLocalizationLLHDist_CB(~,~)
        plotMainAxes(@obj.plotRawLocalizationLLHDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end
    
    %% Phase 6 Plots and Views
    function plotLocalizationSumImage_CB(~,~)
        plotMainAxes(@obj.plotLocalizationSumImage,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos,[10, 10, 0, 0]);
    end

    function plotLocalizationPosDist_CB(~,~)
        plotMainAxes(@obj.plotLocalizationPosDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotLocalizationIDist_CB(~,~)
        plotMainAxes(@obj.plotLocalizationIDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotLocalizationSigmaDist_CB(~,~)
        plotMainAxes(@obj.plotLocalizationSigmaDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotLocalizationLLHDist_CB(~,~)
        plotMainAxes(@obj.plotLocalizationLLHDist,handles.axes);
        GUIBuilder.positionAxes(handles.axes, axes_pos,[10, 10, 0, 0]);
    end

    function plotThreshold_CB(~,~)
       plotMainAxes(@obj.plotThreshold,handles.axes); 
       % We position this axes as if it were a square image.  
       % This will keep the axes square so the threshold line looks correct
       GUIBuilder.positionImageAxes(handles.axes, [1,1], axes_pos,[10, 10, 0, 0]);
    end

    function plot3DTrackSequence_CB(~,~)
        axes(handles.axes);
        cla(handles.axes,'reset');
        obj.plot3DTrackSequence();
        axesMode_CB('rotate3d',[],[]);
        GUIBuilder.positionImageAxes(handles.axes,[1 1], axes_pos,[10, 10, 80, 0]);
    end

    function plotEmitterSuperResGauss_CB(~,~)
        plotMainAxes(@obj.plotEmitterSuperResGauss,handles.axes);
        GUIBuilder.positionImageAxes(handles.axes, flip(obj.frameSize), axes_pos);
        colormap(hot)
    end

    function makeTrackSegmentAnalysis_CB(~,~)
        tsaObj=TrackSegmentAnalysis(obj);
        name=saveGlobalVar('tsaObj',tsaObj);
        fprintf('Saved New TSA object as workspace variable "%s"\n',name);
        tsaObj.gui()
    end
end
