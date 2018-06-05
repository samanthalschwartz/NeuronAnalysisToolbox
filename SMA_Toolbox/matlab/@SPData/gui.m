% SPData.gui
function guiFig=gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end
    gui@GUIBuilder(obj); %Call GUIBuilder initialization
    gui_name='[SPData] Single Particle EMCCD Data';
    
    %% Configure shared local variables.
    lastPath=obj.workingDir; %save this for file open dialogs after a reset
    hJTable=[]; %Handle to java object for ROI table
    
    uH=GUIBuilder.default_unitHeight; % unit height for elements
    boarder=GUIBuilder.default_boarder;%Boarder width around the outside of the gui
    sp=GUIBuilder.default_spacing; %spacing between elements.
    but_sz=GUIBuilder.default_buttonSize; %Button size
%     
    fig_sz=[910 535]; %figure size
    fullw_sz=[fig_sz(1)-2*boarder-2*sp uH]; %size of a full-width component that stays within the boarders.
    halfw=404; %vertical divider position
    halfw_sz=[halfw-2*sp uH];
    roi_table_H = 120; % Height of ROI table

    %% Create figure and controls
    % Make figure
    guiFig = figure('Units','pixels','Position',[200 50 fig_sz],'Resize','off',...
                    'MenuBar','none','ToolBar','none','NumberTitle','off',...
                    'Name',gui_name,'Visible','on',...
                    'CloseRequestFcn',@close_CB);
    obj.guiFig=guiFig;

    % Make Panels 
    fullPan_pos=[boarder boarder fullw_sz];
    halfPan_pos=[boarder boarder halfw_sz];
    handles.panels.file=uipanel('Parent',guiFig,'Units','Pixels','Position',fullPan_pos,'Title','Files');
    handles.panels.cal=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Calibration');
    handles.panels.physical=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Physical Parameters');
    handles.panels.ROI=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','ROI');
    populateFilePanel(handles.panels.file);
    populateCalPanel(handles.panels.cal);
    populatePhysicalPanel(handles.panels.physical);
    populateROIPanel(handles.panels.ROI);

    % Make image axes
    handles.axes = axes('Units', 'Pixels','Visible','off');

    %% initialization
    alignPanels();
    drawnow();
    createMenus();
    positionAxes(handles.axes);
    setROITableData();
    obj.makeTableRowSelectable(handles.ROITable, @(~,~) displayROI);
    updateControls();
    positionAxes(handles.axes);
    pause(0.05); %Sometimes we have a draw error maybe we need to wait for Java EDT???
    drawnow();



    %% Panel creation functions
    function populateFilePanel(panH)
        % File Panel
        % Bottom row of buttons
        bot_row.pos = [sp, sp, fullw_sz];
        bot_row.names = {'Load','Save','SaveAs','SetDefaultProps', 'BatchProcess'};
        bot_row.CBs = {@load_CB,@save_CB,@saveAs_CB, @setDefaultProps_CB, @batchProcess_CB};
        GUIBuilder.buttonRow(panH, bot_row.pos, but_sz+[20 0], bot_row.names, bot_row.CBs);
        %Quit is out to the side so we don't include it.
        handles.quitButton = uicontrol('Parent',panH,'Style','pushbutton','String','Quit',...
            'Position',[fullw_sz(1)-but_sz(1) sp but_sz],...
            'Callback',@close_CB);
        edits.pos = bot_row.pos + [0, bot_row.pos(2)+but_sz(2)+sp, 0, bot_row.pos(4)+2*sp];
        edits.hNames = {'rawDataPath','saveFilePath'}; %names in the handles structure returned
        edits.labels = {'Raw Data Path', 'SPData Path'};
        edits.values = {obj.rawDataPath, obj.saveFilePath};
        edits.CBs = {0,0};
        
        handles.file = GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
    end

    function populateCalPanel(panH)
        % Bottom row of buttons
        bot_row.pos = [sp, sp, halfw_sz];
        bot_row.names = {'Calibrate'};
        bot_row.CBs = {@gainCal_CB};
        GUIBuilder.buttonRow(panH, bot_row.pos, but_sz, bot_row.names, bot_row.CBs);
        
        edits.pos = [sp+but_sz(1)+sp, sp, halfw_sz-but_sz(1)-sp];
        edits.hNames = {'CCDGain','CCDBackground'}; %names in the handles structure returned
        edits.labels = {'CCD Gain [e-/ADU]', 'CCD Background'};
        edits.values = {obj.CCDGain, obj.CCDBackground};
        edits.CBs = {@setCCDGain_CB, @setCCDBackground_CB};
        handles.cal = GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
    end

    function populatePhysicalPanel(panH)
        edits.pos = [sp, sp, halfw_sz];
        edits.hNames = {'size','nFrames','globalTBounds', 'pixelSize','psf','frameT'};
        edits.labels = {'Size[X,Y]:', 'Num Frames:', 'Global Time Bounds:', 'Pixel Size (um)',...
            'Point Spread Function (px)','Frame Time (s)'};
        edits.values = {obj.frameSize, obj.nFrames, obj.globalTBounds, obj.pixelSize, obj.psf, obj.frameT};
        edits.CBs = {[],[],@setGlobalTBounds_CB, @setPixelSize_CB, @setPSF_CB,@setFrameT_CB};
        handles.physical = GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels, edits.values, edits.CBs);
    end

    function populateROIPanel(panH)
        % File Panel
        % Bottom row of buttons
        column_width = [125 cellfun(@(~) 40,cell(1,6),'Uniform',0)];
        ROIcols = {'ROI Name','xmin','xmax','ymin','ymax','tmin','tmax'};
        handles.ROITable = uitable(panH,'ColumnName',ROIcols,...
                                        'Position', [sp sp halfw_sz(1) roi_table_H],...
                                        'ColumnWidth', column_width,...
                                        'ColumnEditable',true(1,7),...
                                        'CellEditCallback',@ROIEdit_CB);
        
    end

    function createMenus()
        %File Menu
        file_labels = {'Load ...','Save','Save As ...','Reset Object',[],...
                       'Save Data As ...','Convert -> MAT ...',[],...
                       'Set Default Properties ...', 'Batch Process...',[],...
                       'Quit'};
        file_CBs = {@load_CB, @save_CB, @saveAs_CB, @resetObject_CB,[],...
                    @saveDataAs_CB, @convertToMat_CB,[],...
                    @setDefaultProps_CB,@batchProcess_CB,[],...
                    @close_CB};
        obj.makeFigureMenu(guiFig, 'File',file_labels, file_CBs);
        
        %View Menu
        view_labels = {'View Full Frames','View Full Raw Frames'};
        view_CBs = {@viewFrames_CB, @(~,~) obj.viewRawFrames()};
        GUIBuilder.makeFigureMenu(guiFig,'View',view_labels, view_CBs);
        
        %ROI Menu
        roi_labels = {'View ROI Frames','Track ROI','Batch Track ROI',[],'Add ROI', 'Add ROI DIP', 'Delete ROI', 'Clear ROI'};
        roi_CBs = {@viewROI_CB, @roiRPT_CB, @roiBatchRPT_CB,[],@addROI_CB,@addROIdip_CB,@delROI_CB,@clearROI_CB};
        GUIBuilder.makeFigureMenu(guiFig,'ROI',roi_labels, roi_CBs);
        
        %Analysis Menu
        analysis_labels = {'RPT Track ROI','SPT Track ROI','TrackSegmentAnalysis for ROI','View Figures for ROI'};
        analysis_CBs = {@roiRPT_CB,@roiSPT_CB,@roiTSA_CB,@viewFigs_CB};
        GUIBuilder.makeFigureMenu(guiFig,'Analysis',analysis_labels, analysis_CBs);
        
        % ROITable menu
        labels = {'View ROI','Track ROI','Batch Track ROI',[],'Add ROI','Add ROI DIP','Delete ROI'};
        CBs = {@viewROI_CB, @roiRPT_CB, @roiBatchRPT_CB, [], @addROI_CB,@addROIdip_CB, @delROI_CB};
        handles.ROITable.UIContextMenu = GUIBuilder.makeContextMenu(labels, CBs);
        
        %This menu will be reused for all ROI patch objects
        labels = {'View ROI','Track ROI','Batch Track ROI',[], 'Add ROI', 'Delete ROI'};
        CBs = {@viewROI_CB, @roiRPT_CB, @roiBatchRPT_CB, [], @addROI_CB, @delROI_CB};
        handles.contextmenus.roicm_h = GUIBuilder.makeContextMenu(labels, CBs);
        
        %This menu will be reused for all image context menus
        labels = {'Add ROI','Add ROI DIP'};
        CBs = {@addROI_CB, @addROIdip_CB};
        handles.contextmenus.roiimage_h = GUIBuilder.makeContextMenu(labels, CBs);
    end


    %% Helper functions for Alignment and updating of elements
    function panels_pos=alignPanels()
        structfun(@GUIBuilder.autoSizePanel, handles.panels,'Uniform',0); %Auto size each panel in the handles.panels struct
        panels_pos=GUIBuilder.align(struct2cell(handles.panels),'Left','Fixed',sp);
    end
    
    function positionAxes(axesH)
        %Calculate the axes area for drawing into and then fit an image of given frame size inside
        cal_pos = round(handles.panels.cal.Position);
        axes_pos = [cal_pos(1)+cal_pos(3), cal_pos(2), 0 0];
        axes_pos = axes_pos+[0, 0, fig_sz(1)-axes_pos(1), fig_sz(2)-axes_pos(2)];
        GUIBuilder.positionImageAxes(axesH,obj.frameSize, axes_pos,[10 10 0 0]);
    end
   
    function setROITableData()
        handles.ROITable.Data = [obj.ROIname', num2cell(cell2mat(obj.ROI'))];
        drawnow(); % Need to call this for the JTable java object to be created
        hJTable = obj.getJTableHandle(handles.ROITable);
    end
    
    function updateControls()
        handles.file.saveFilePath.String = obj.saveFilePath;
        handles.file.rawDataPath.String = obj.rawDataPath;
        handles.physical.size.String = arr2str(obj.frameSize);
        handles.physical.nFrames.String = arr2str(obj.nFrames);
        obj.setDefaultableControl(handles.physical.globalTBounds,'globalTBounds');
        obj.setDefaultableControl(handles.physical.frameT,'frameT');
        obj.setDefaultableControl(handles.physical.psf,'psf');
        obj.setDefaultableControl(handles.physical.pixelSize,'pixelSize');
        obj.setDefaultableControl(handles.cal.CCDBackground,'CCDBackground');
        obj.setDefaultableControl(handles.cal.CCDGain,'CCDGain');
        
        handles.SPT.sptParamPath.String = obj.getFilePath('SPTParams');

        if obj.dirty
            handles.file.saveFilePath.BackgroundColor = obj.color_unsavedBG;
        else
            handles.file.saveFilePath.BackgroundColor = obj.color_editBG;
        end
        if ~isempty(obj.ROI) && getSelectedROI()==0
            setSelectedROI(1);
        end
        displayROI();
    end

    function displayROI()
        if ~obj.initialized 
            return
        end

        axes(handles.axes);
        imH=obj.plotSumImage();
        imH.ButtonDownFcn=@roiButtonDown_CB;
        imH.UIContextMenu=handles.contextmenus.roiimage_h;
        if ~isempty(obj.ROI)
            selectedIdx=getSelectedROI();
            for n=1:length(obj.ROI)
                roi=obj.ROI{n};
                roi = roi -0.5; %this 0.5 pixel difference is between world and intrinsic coordinates.
                lc=obj.colors(mod(n-1,size(obj.colors,1))+1,:);
                rw=roi(2)-roi(1)+1;
                rh=roi(4)-roi(3)+1;
                buf_w=obj.sizeX/20;
                th=text(roi(1)+0.5*rw-buf_w, roi(3)+0.5*rh-3, obj.ROIname{n},...
                     'Color',lc,'FontSize',10,'FontWeight','bold');
                h=patch([roi(1) roi(2) roi(2) roi(1)], [roi(3) roi(3) roi(4) roi(4)], lc,...
                        'LineWidth',1,'EdgeColor',lc,'FaceAlpha',0,'HitTest','on');
                if n==selectedIdx
                    set(h,'FaceColor',lc,'LineWidth',2,'FaceAlpha',0.1);
                end
                h.ButtonDownFcn = @roiButtonDown_CB;
                h.UIContextMenu = handles.contextmenus.roicm_h;
                th.ButtonDownFcn = @roiButtonDown_CB;
                th.UIContextMenu = handles.contextmenus.roicm_h;
            end
        end
        set(handles.axes,'Visible','on','HitTest','on');
        drawnow();
    end
    
    function roi_idx=getSelectedROI()
        obj.assertInitialized()
        roi_idx = hJTable.getSelectedRow + 1;
        if roi_idx==0 && ~isempty(obj.ROI)
            setSelectedROI(1);
        end
        roi_idx = hJTable.getSelectedRow + 1;
        if roi_idx==0 || roi_idx>length(obj.ROI)
            error('SPData:gui', 'Unable to select an ROI because no ROI have been defined.\nCreate an ROI first.');
        end
    end

    function setSelectedROI(roi_idx)
        if roi_idx>length(obj.ROI)
            roi_idx=-1;
        else
            roi_idx=roi_idx-1;
        end
        hJTable.changeSelection(roi_idx,-1, false, false);
    end
      
    %% Callbacks -- File Panel
    function load_CB(~, ~)
        if ~isempty(obj.workingDir)
            lastPath = obj.workingDir;
        end
        [filename, pathname]=uigetfile(obj.LoadableDataFormats,'Select data source', lastPath);
        if ~filename; return; end
        obj.load(fullfile(pathname,filename))
        obj.inGui=true;
        setROITableData();
        positionAxes(handles.axes)
        updateControls();
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
        setROITableData();
        updateControls();
    end

    function saveDataAs_CB(~,~)
        obj.saveDataAs();
        updateControls();
    end

    function convertToMat_CB(~, ~)
        obj.convertToMat();
        updateControls();
    end
    
    function setDefaultProps_CB(~,~)
        [filename, pathname]=uigetfile(obj.SaveableDataFormats,'Select default file', obj.saveFilePath);
        if ~filename; return; end
        defaultFile = fullfile(pathname,filename);
        if strcmp(defaultFile, obj.saveFilePath)
            error('SPData:gui:setDefaultProps_CB','Must choose different file for a default parameters file');
        end
        spd = SPData(defaultFile);
        obj.setPreservedProperties(spd);
        setROITableData();
        updateControls();
    end

    function batchProcess_CB(~,~)
        if obj.dirty
            error('SPData:gui:batchProcess_CB','Object is dirty. Try saving first.');
        end
        filepaths = Pickle.selectBatchProccesingFileNames(obj.workingDir,obj.RawDataFormats{1},obj.LoadableRawDataFormats, 'Select data to batch process as .spdata');
        if isempty(filepaths); return; end
        overwriteFlag = true; % will still ask for confirmation if necessary
        new_files = SPData.batchProcess(path,filepaths,obj,overwriteFlag);
        if isempty(new_files); return; end
        count = sum(cellfun(@ischar, new_files));
        msgbox(sprintf('Succesfully processed %i/%i %s files.',count,numel(new_files),obj.saveFileExt));
    end

    function close_CB(~,~)
        try
            obj.clearJTableCallbacks(handles.ROITable);
        catch err
            disp(getReport(err))
        end
        obj.closeGUI();
    end

    %% Callbacks -- Calibration Panel
    function setCCDBackground_CB(H,~)
        bg = str2double(H.String);
        if ~isfinite(bg) || bg<0
            H.String = obj.CCDBackground;
            warning('SPData:gui','Invalid CCD Background: %s', H.String);
        elseif ~isempty(obj.CCDGain)
            obj.recalibrateFrames(bg, obj.CCDGain);
            updateControls();
        else
            obj.CCDBackground = bg; %Just set the background because gain is not set yet.
        end
    end

    function setCCDGain_CB(H,~)
        gain = str2double(H.String);
        if ~isfinite(gain) || gain<=0
            H.String = obj.CCDGain;
            warning('SPData:gui','Invalid CCD Gain: %s', H.String);
        elseif ~isempty(obj.CCDBackground)
            obj.recalibrateFrames(obj.CCDBackground, gain);
            updateControls();
        else
            obj.CCDGain = gain; %Just set the gain because background is not set yet.
        end
    end

    function gainCal_CB(~,~)
        obj.gainCalGui(@updateControls);
    end

    %% Callbacks -- Physical Panel
    function setGlobalTBounds_CB(H,~)
        try
            obj.setGlobalTBounds(str2num(H.String));
        catch err
            H.String = obj.globalTBounds;
            throw(err);
        end
        updateControls();    
        setROITableData();
    end
    
    function setFrameT_CB(H,~)
        frameT = str2double(H.String);
        if isfinite(frameT) && frameT>0
            obj.frameT = frameT;
            obj.dirty = true;
        else
            H.String = obj.frameT;
            warning('SPData:gui','Invalid frame time: %s',H.String);
        end
    end

    function setPSF_CB(H,~)
        psf=str2num(H.String); %#ok<*ST2NM>
        if (isscalar(psf) || length(psf)==2) && all(isfinite(psf)) && all(psf>0)
            obj.psf=psf;
            obj.dirty=true;
        else
            H.String = arr2str(obj.psf);
            warning('SPData:gui','Invalid psf: %s',H.String);
        end
    end

    function setPixelSize_CB(H,~)
        pixelsize = str2num(H.String);
        if (isscalar(pixelsize) || length(pixelsize)==2) && all(isfinite(pixelsize)) && all(pixelsize>0)
            obj.pixelSize=pixelsize;
            obj.dirty=true;
        else
            H.String = arr2str(obj.pixelSize); %Reset the string to the original
            warning('SPData:gui','Invalid pixel size: %s',H.String);
        end
    end
    
    %% Callbacks -- ROI Panel
    function addROIdip_CB(~,~)
        if obj.calibrated
            frames=obj.getFrames();
        else
            frames=obj.getRawFrames();
        end
        if isempty(frames)
            error('SPData:gui','Unable to load raw frames');
        end
        h = obj.viewMaximizedDipFig(frames);
        [~,bbox] = dipcrop(h);
        roi = [bbox(1) bbox(1)+bbox(2) bbox(3) bbox(3)+bbox(4)]+1;
        obj.addROI(roi);
        setROITableData();
        setSelectedROI(length(obj.ROI));
        updateControls();
        close(h);
    end

    function addROI_CB(~,~)
        rect=getrect(handles.axes);
        roi= round([rect(1), rect(1)+rect(3), rect(2), rect(2)+rect(4)]);
        obj.addROI(roi);
        setROITableData();
        setSelectedROI(length(obj.ROI));
        updateControls();
    end

    function delROI_CB(~,~)
        roi_idx=getSelectedROI();
        obj.deleteROI(roi_idx);
        setROITableData();        
        setSelectedROI(min(length(obj.ROI),roi_idx));
        updateControls();
    end

    function clearROI_CB(~,~)
        if strcmp('Yes',questdlg('Really clear all ROI?'))
            obj.clearROI();
            setROITableData();
        end
        updateControls();
    end

    function viewROI_CB(~,~)
        roi_idx = getSelectedROI();
        if obj.calibrated
            obj.viewFrames(roi_idx);
        else
            warning('SPData:gui:viewROI','Object is not calibrated.  Viewing Raw frames.');
            obj.viewRawFrames(roi_idx);
        end
    end

    function ROIEdit_CB(~,ev)
        roi_idx = ev.Indices(1);
        edit_col = ev.Indices(2);
        roi=obj.ROI{roi_idx};
        roi_name=obj.ROIname{roi_idx};
        if edit_col==1
            roi_name=ev.NewData;
        else
            roi(edit_col-1)=ev.NewData;
        end
        obj.modifyROI(roi_idx, roi, roi_name);
        setROITableData();
        setSelectedROI(roi_idx);
        updateControls();
    end

    function roiButtonDown_CB(hObj, ~)
        nROI=length(obj.ROI);
        if nROI<2
            return
        end
        select_type=get(gcbf(),'SelectionType');
        switch select_type
            case {'normal','open'}
                cp = hObj.Parent.CurrentPoint;
                cp = cp(1,1:2);
                selROI = getSelectedROI();
                for N=1:nROI
                    n=mod(N+selROI-1,nROI)+1; % n will iterate through all ROI starting at selected
                    roi=obj.ROI{n};
                    if cp(1)>=roi(1) && cp(1)<=roi(2) && cp(2)>=roi(3) && cp(2)<=roi(4)
                        setSelectedROI(n);
                        updateControls();
                        return
                    end
                end
        end
    end

    function viewFrames_CB(~,~)
        obj.viewFrames();
    end

    %% Callbacks - RPT Analysis   
    function roiRPT_CB(~,~)
        rpt=obj.trackRPT(getSelectedROI());
        rpt.gui();
        name=saveGlobalVar('rptObj',rpt);
        fprintf('Saved RPT object as workspace variable "%s"\n',name);
    end

    function roiBatchRPT_CB(~,~)
        roi_idx = getSelectedROI();
        roi_name = obj.ROIname{roi_idx};
        selectSingle = true; % Ask if there are more than one tracked file for this data
        default_rpt_file = obj.getROIFiles(roi_idx, 'RPT', selectSingle);
        if isempty(default_rpt_file)
            default_rpt_file = Pickle.selectExistingFileName(obj.workingDir,RPT.saveFileExt,RPT.SaveableDataFormats,...
                                        'Select .RPT with desired Default Properties');
            if isempty(default_rpt_file); return; end
        end
        title = sprintf('Select SPData files for batch processing of ROI: "%s"',roi_name);
        spdata_files = obj.selectBatchProccesingFileNames(obj.workingDir, '*.spdata',SPData.SaveableDataFormats, title);
        if isempty(spdata_files); return; end
        overwriteFlag=1;
        rpt_files = RPT.batchProcess(obj.workingDir, spdata_files, roi_name, default_rpt_file, overwriteFlag);
        if isempty(rpt_files); return; end
        count = sum(cellfun(@ischar, rpt_files));
        msgbox(sprintf('Succesfully processed %i/%i .rpt files for ROI: "%s" ',count,numel(rpt_files), roi_name));
    end

    function viewFigs_CB(~,~)
        roi_path=obj.getROIPath(getSelectedROI());
        if ~exist(roi_path,'dir')
            errordlg(sprintf('No analysis folder for ROI at "%s"',roi_path),'Error viewing Figures');
            return
        end
        [filen, pathn]=uigetfile({'*.fig','Matlab Figure (.fig)'}, 'Open figure to view', roi_path);
        if ~filen; return; end
        obj.appendOpenFigs(openfig(fullfile(pathn,filen)));
    end
end


