function gui(obj)
    % HSData GUI
    %
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end
    gui@GUIBuilder(obj); %Call GUIBuilder initialization

    %% Configure shared local variables.
    lastPath=obj.workingDir; %save this for file open dialogs after a reset
    roi_hJTable=[]; %Handle to java object for ROI table
    uH=25; % unit height for elements
    fig_sz=[1000 685]; %figure size
    boarder=GUIBuilder.default_boarder;%Boarder width around the outside of the gui
    but_sz=[100 uH]; %Button size
    sp=3; %spacing between elements.
    fullw_sz=[fig_sz(1)-2*boarder-2*sp uH]; %size of a full-width component that stays within the boarders.
    halfw=485;
    halfw_sz=[halfw-2*sp uH];
    bg_color=[0.2 0.2 0.2];
    roi_table_H = 120; % Height of ROI table

    %% Create figure and controls
    % Make figure
    guiFig = figure('Units','pixels','Position',[200 50 fig_sz],'Resize','off',...
                    'MenuBar','none','ToolBar','none','NumberTitle','off',...
                    'Name','HSData GUI','Visible','on',...
                    'CloseRequestFcn',@close_CB);
    obj.guiFig=guiFig;
    %Change the background colors
    whitebg(obj.guiFig);

    % Make Panels 
    fullPan_pos=[boarder boarder fullw_sz];
    halfPan_pos=[boarder boarder halfw_sz];
    handles.panels.file=uipanel('Parent',guiFig,'Units','Pixels','Position',fullPan_pos,'Title','Files');
    handles.panels.cal=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Calibration');
    handles.panels.physical=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','Physical Parameters');
    handles.panels.ROI=uipanel('Parent',guiFig,'Units','Pixels','Position',halfPan_pos,'Title','ROI');
    populateFilePanel(handles.panels.file);
    populatePhysicalPanel(handles.panels.physical);
    populateCalPanel(handles.panels.cal);
    populateROIPanel(handles.panels.ROI);

    % Display axes
    handles.axes.sum= axes('Units', 'Pixels','Visible','off');
    handles.axes.spectra= axes('Units', 'Pixels','Visible','off');
    handles.axes.spectralseries= axes('Units', 'Pixels','Visible','off');

    %% initialization
    alignPanels();
    createMenus();
    positionAxes();
    if obj.initialized
        obj.plotSpectralSeries(handles.axes.spectralseries);
        view(0,90);
    end
    setROITableData();
    positionAxes();
    obj.makeTableRowSelectable(handles.ROITable,@updateROISelection_CB);
    updateControls();

    %% Panel creation functions
    function populateFilePanel(panH)
        % File Panel
        % Bottom row of buttons
        bot_row.pos=[sp, sp, fullw_sz];
        bot_row.names={'Load','Save','SaveAs','Reset'};
        bot_row.CBs={@load_CB,@save_CB,@saveAs_CB,@reset_CB};
        GUIBuilder.buttonRow(panH, bot_row.pos, but_sz, bot_row.names, bot_row.CBs);
        %Quit is out to the side so we don't include it.
        handles.quitButton = uicontrol('Parent',panH,'Style','pushbutton','String','Quit',...
            'Position',[fullw_sz(1)-but_sz(1) sp but_sz],...
            'Callback',@close_CB);
        edits.pos=bot_row.pos;
        edits.pos(2)=edits.pos(2)+but_sz(2)+sp;
        edits.pos(4)=edits.pos(4)+2*sp;
        edits.hNames={'hsdataPath'}; %names in the handles structure returned
        edits.labels={'HSData Path:'};
        edits.values={obj.saveFilePath};
        edits.CBs={@saveAs_CB};
        
        handles.file=GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels,...
            edits.values, edits.CBs);
    end

    function populateCalPanel(panH)
        % File Panel
        % Bottom row of buttons
        edits.pos=[sp, sp, halfw_sz];
        edits.hNames={'CCDGain','CCDBackground', 'CCDCalDate','WVCalDate' }; %names in the handles structure returned
        edits.labels={'CCD Gain [e-/ADU]', 'Mean CCD Background', 'CCD Calibration Date:', 'Wavelength Calibration Date:'};
        edits.values={obj.CCDGain, mean(obj.CCDBackground(:)),...
                      datestr(obj.AcquisitionParams.gainCalDateTime,31),...
                      datestr(obj.AcquisitionParams.wvCalDateTime,31)};
        edits.CBs={[], [], [], []};
        
        handles.cal=GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels,...
            edits.values, edits.CBs);
    end

    function populatePhysicalPanel(panH)
        edits.pos=[sp, sp, halfw_sz];
        edits.hNames={'acqDate','size', 'nFrames', 'physicalSize', 'nFramePixels','sensorT','frameT',...
                      'globalFrameTime','pixelSize','psf',};
        edits.labels={'Acquisition Date:', 'Size [LYX] (px):','Num Frames:',...
                      'PhysicalSize [XY] (micron):','# Frame Pixels:',...
                      'Sensor Time [LY] (s):', 'Frame Time [LYX] (s):',...
                      'Global Time Bounds (frames):', 'Pixel Size (um):',...
                      'Point Spread Function  [XY] (px):'};
        edits.values={datestr(obj.AcquisitionParams.acqDateTime,31), [obj.sizeX, obj.sizeY, obj.sizeL],...
                      obj.nFrames, obj.physicalSize, obj.NFramePixels,...
                      obj.sensorT, obj.frameT, obj.globalTBounds, obj.pixelSize, obj.psf};
        edits.CBs={[],[],[],[],[],[],[],@setGlobalTimeBounds_CB, @setPixelSize_CB, @setPSF_CB};
        handles.physical=GUIBuilder.labeledHEdits(panH, edits.pos, uH, edits.hNames, edits.labels,...
                                                  edits.values, edits.CBs);
    end

    function populateROIPanel(panH)
        % File Panel
        % Bottom row of buttons
        column_width=[125 cellfun(@(~) 40, cell(1,8),'Uniform',0)];
        ROIcols={'ROI Name','xmin','xmax','ymin','ymax','Lmin','Lmax','tmin','tmax'};
        handles.ROITable=uitable(panH,'ColumnName',ROIcols,...
                                      'Position', [sp sp halfw_sz(1) roi_table_H],...
                                      'ColumnWidth', column_width,...
                                      'ColumnEditable',true(1,9),...
                                      'CellEditCallback',@ROIEdit_CB);    
    end

    function createMenus()
        %File Menu
        labels = {'Load...', 'Save', 'Save As ...', 'Reset Object',[], 'Quit'};
        CBs = {@load_CB, @save_CB, @saveAs_CB, @resetObject_CB,[], @close_CB};
        GUIBuilder.makeFigureMenu(guiFig,'File',labels,CBs);
        
        %View Menu
        labels = {'View Full Frames RGB (3D) ...','View Full Frames Grayscale (4D) ...',  'View Full Raw Frames Grayscale (4D) ...'};
        CBs = {@(~,~) obj.viewFramesRGB(), @(~,~) obj.viewFrames(),@(~,~) obj.viewRawFrames() };
        GUIBuilder.makeFigureMenu(guiFig,'View',labels,CBs); 

        %ROI Menu
        labels = {'Track ROI','Batch Track ROI',[],...
                  'View Frames RGB (3D) [XYT] ...', 'View Frames Grayscale (4D) [LYXT] ...',...
                  'View Frames Volumes Intensity (4D) [XYLT] ...', 'View Frames Volumes Wavelength (4D) [XYLT] ...',[],...
                  'View Sum Image RGB (2D) [XY] ...','View Sum Image Gray (2D) [XY] ...',...
                  'View Sum Image Slices Intensity (3D) [XYL] ...','View Sum Image Slices Wavelength (3D) [XYL]',...
                  'View Sum Image Volumes Intensity (3D) [XYL]', 'View Sum Image Volumes Wavelength (3D) [XYL]', [],...
                  'View ROI Spectra (1D) [L] ...', 'View ROI Spectral Series (2D) [LT] ...',[],...
                  'Add ROI','Delete ROI','Clear ALL ROI'};
        CBs = {@trackHSRPT_ROI_CB, @roiBatchTrackHSRPT_CB,[],...
               @viewROI_RGB_CB,@viewROI_4D_CB,...
               @viewROI_Frames_Volumes_Intensity_CB, @viewROI_Frames_Volumes_Wavelength_CB,[],...
               @viewROI_Sum_RGB_CB,@viewROI_Sum_Gray_CB,...
               @viewROI_Sum_Slices_Intensity_CB, @viewROI_Sum_Slices_Wavelength_CB, @viewROI_Sum_Volumes_Intensity_CB, @viewROI_Sum_Volumes_Wavelength_CB,[],...
               @viewROI_Spectra_CB,@viewROI_SpectralSeries_CB,[],...
               @addROI_CB, @delROI_CB, @clearROI_CB};
        GUIBuilder.makeFigureMenu(guiFig,'ROI',labels,CBs);

        % ROITable menu
        labels = {'View ROI RGB (3D) [XYT] ...',[],...
                  'Track ROI','Batch Track ROI',[],'Select Spectra','Add ROI','Delete ROI', 'Clear ALL ROI'};
        CBs = {@viewROI_RGB_CB,[],@trackHSRPT_ROI_CB, @roiBatchTrackHSRPT_CB, [],@selectSpectra_CB, @addROI_CB, @delROI_CB, @clearROI_CB};
        handles.ROITable.UIContextMenu = GUIBuilder.makeContextMenu(labels, CBs);

        %This menu will be reused for all ROI patch objects
        labels = {'Track ROI','Batch Track ROI',[],...
                  'View Frames RGB (3D) [XYT] ...', 'View Frames Grayscale (4D) [LYXT] ...',...
                  'View Frames Volumes Intensity (4D) [XYLT] ...', 'View Frames Volumes Wavelength (4D) [XYLT] ...',[],...
                  'View Sum Image RGB (2D) [XY] ...','View Sum Image Gray (2D) [XY] ...',...
                  'View Sum Image Slices Intensity (3D) [XYL] ...','View Sum Image Slices Wavelength (3D) [XYL]',...
                  'View Sum Image Volumes Intensity (3D) [XYL]', 'View Sum Image Volumes Wavelength (3D) [XYL]', [],...
                  'View ROI Spectra (1D) [L] ...', 'View ROI Spectral Series (2D) [LT] ...',[],...
                  'Add ROI','Delete ROI'};
        CBs = {@trackHSRPT_ROI_CB, @roiBatchTrackHSRPT_CB,[],...
               @viewROI_RGB_CB,@viewROI_4D_CB,...
               @viewROI_Frames_Volumes_Intensity_CB, @viewROI_Frames_Volumes_Wavelength_CB,[],...
               @viewROI_Sum_RGB_CB,@viewROI_Sum_Gray_CB,...
               @viewROI_Sum_Slices_Intensity_CB, @viewROI_Sum_Slices_Wavelength_CB, @viewROI_Sum_Volumes_Intensity_CB, @viewROI_Sum_Volumes_Wavelength_CB,[],...
               @viewROI_Spectra_CB,@viewROI_SpectralSeries_CB,[],...
               @addROI_CB, @delROI_CB};
        handles.contextmenus.roicm_h = GUIBuilder.makeContextMenu(labels, CBs);


        %Context menu for blank parts of ROI image
        labels = {'Add ROI','Clear ALL ROI'};
        CBs = {@addROI_CB, @clearROI_CB};
        handles.contextmenus.roiimage_h=GUIBuilder.makeContextMenu(labels, CBs);

        %Context menu for spectral patches
        labels = { 'Select Spectra',[],...
                   'View Frames RGB (3D) [XYT] ...', 'View Frames Grayscale (4D) [LYXT] ...',...
                   'View Frames Volumes Intensity (4D) [XYLT] ...', 'View Frames Volumes Wavelength (4D) [XYLT] ...',[],...
                   'View Sum Image RGB (2D) [XY] ...','View Sum Image Gray (2D) [XY] ...',...
                   'View Sum Image Slices Intensity (3D) [XYL] ...','View Sum Image Slices Wavelength (3D) [XYL]',...
                   'View Sum Image Volumes Intensity (3D) [XYL]', 'View Sum Image Volumes Wavelength (3D) [XYL]', [],...
                   'View ROI Spectra (1D) [L] ...', 'View ROI Spectral Series (2D) [LT] ...'};
        CBs = {@selectSpectra_CB,[], ...
               @viewROI_RGB_CB,@viewROI_4D_CB,...
               @viewROI_Frames_Volumes_Intensity_CB, @viewROI_Frames_Volumes_Wavelength_CB,[],...
               @viewROI_Sum_RGB_CB,@viewROI_Sum_Gray_CB,...
               @viewROI_Sum_Slices_Intensity_CB, @viewROI_Sum_Slices_Wavelength_CB, @viewROI_Sum_Volumes_Intensity_CB, @viewROI_Sum_Volumes_Wavelength_CB,[],...
               @viewROI_Spectra_CB,@viewROI_SpectralSeries_CB};
        handles.contextmenus.spectra_h = GUIBuilder.makeContextMenu(labels, CBs);
    end


    %% Helper functions for Alignment and updating of elements
    function panels_pos=alignPanels()
        pans=fieldnames(handles.panels);
        for n=1:length(pans)
            GUIBuilder.autoSizePanel(handles.panels.(pans{n}));
        end
        panels_pos=GUIBuilder.align(struct2cell(handles.panels),'Left','Fixed',sp);
    end

    function positionAxes()
        file_pos=round(get(handles.panels.file,'Position'));
        phys_pos=round(get(handles.panels.physical,'Position'));
        minH=file_pos(2)+file_pos(4)+sp;
        maxH=fig_sz(2)-boarder;
        minW=phys_pos(1)+phys_pos(3)+sp;
        maxW=fig_sz(1)-sp;
        W=maxW-minW;
        H=maxH-minH;
        Hsp=round(H/12);
        axes_pos=[minW minH W 3*Hsp;
                  minW minH+3*Hsp W 3*Hsp;
                  minW minH+6*Hsp W 6*Hsp];
        GUIBuilder.positionImageAxes(handles.axes.sum,[obj.sizeX obj.sizeY], axes_pos(3,:),[20 20 0 0]);
        GUIBuilder.positionAxes(handles.axes.spectra,axes_pos(2,:),[20 20 0 0]);
        GUIBuilder.positionAxes(handles.axes.spectralseries,axes_pos(1,:),[0 0 0 0]);
        rotate3d('off');
        handles.rotateButton = uicontrol('Style','togglebutton','CData',GUIBuilder.readMatlabIcon('tool_rotate_3d.png'),...
            'Position',[axes_pos(1,1:2) 28 24],...
            'Callback',@rotateSpectra_CB);
    end

    
    function setROITableData()
        table_data=[obj.ROIname', num2cell(cell2mat(obj.ROI'))];
        set(handles.ROITable, 'Data',table_data);
        drawnow;
        roi_hJTable = obj.getJTableHandle(handles.ROITable);
    end

    function updateControls()
        obj.setDefaultableControl(handles.file.hsdataPath,'saveFilePath');
        obj.setDefaultableControl(handles.physical.globalFrameTime,'globalTBounds');
        obj.setDefaultableControl(handles.physical.pixelSize,'pixelSize');
        obj.setDefaultableControl(handles.physical.psf,'psf');
        if obj.dirty
            set(handles.file.hsdataPath,'BackgroundColor',obj.color_unsavedBG);
        else
            set(handles.file.hsdataPath,'BackgroundColor',obj.color_editBG);
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
        axs=handles.axes.sum;
        selectedIdx=getSelectedROI();
        axes(axs);
        handles.sumimage=image([0.5, obj.sizeX-0.5], [0.5 obj.sizeY-0.5], obj.sumImageRGB,'AlphaData',0.99,'ButtonDownFcn',@roiButtonDown_CB);
        xlabel('X (px)');
        ylabel('Y (px)');
        if ~isempty(obj.ROI)
            idxs = [setdiff(1:numel(obj.ROI), selectedIdx), selectedIdx]; % make selected idx draw last
        else
            idxs = [];
        end
        for n=idxs
            roi=obj.ROI{n};
            lc=obj.colors(mod(n-1,size(obj.colors,1))+1,:);
            lw=2.0;
            rx=roi(2)-roi(1);
            ry=roi(4)-roi(3);
            
            text(roi(1)+0.5*rx-5,roi(3)+0.5*ry-5,obj.ROIname{n},...
                 'Color',lc,'FontSize',10,'FontWeight','bold','HitTest','off','PickableParts','none');
            h=patch([roi(1)-1 roi(2) roi(2) roi(1)-1],[roi(3)-1 roi(3)-1 roi(4) roi(4)],lc,...
                    'LineWidth',1,'EdgeColor',lc,'FaceAlpha',0,'HitTest','on','Clipping','off');
            if n==selectedIdx
                set(h,'FaceColor',lc,'LineWidth',lw,'FaceAlpha',0.1);
            end
            set(h,'ButtonDownFcn', @roiButtonDown_CB);
            set(h,'uicontextmenu',handles.contextmenus.roicm_h);
        end
        set(axs,'Visible','on','HitTest','on');
%         set(axs,'uicontextmenu',handles.contextmenus.roiimage_h);      
        set(handles.sumimage,'uicontextmenu',handles.contextmenus.roiimage_h);
        set(handles.sumimage,'Visible','on','HitTest','on','PickableParts','all');

        obj.plotSpectra(handles.axes.spectra);
        for n=1:length(obj.ROI)
            if n==selectedIdx
                roi=obj.ROI{n};
                lc=obj.colors(mod(n-1,size(obj.colors,1))+1,:);
                lw=2.0;
                rL=roi(6)-roi(5);
                maxL=obj.lambda(roi(5));
                minL=obj.lambda(roi(6));
                yL=get(gca,'YLim');
                height=diff(yL);
                text(minL+(maxL-minL)*.5-10, yL(1)+height*.5 ,obj.ROIname{n},...
                 'Color',lc,'FontSize',10,'FontWeight','bold','HitTest','off');
                h=patch([minL maxL maxL minL],[yL(1) yL(1) yL(1)+height yL(1)+height],...
                        lc,'LineWidth',1,'EdgeColor',lc,'FaceAlpha',0.05,'HitTest','on','Clipping','off');
                set(h,'uicontextmenu',handles.contextmenus.spectra_h);
            end
        end
    end


    
    function row=getSelectedROI()
        row = roi_hJTable.getSelectedRow + 1;
        if ~isempty(obj.ROI) && row==0
            setSelectedROI(1);
        end
    end

    function setSelectedROI(roi_idx)
        if roi_idx>length(obj.ROI)
            roi_idx=-1;
        else
            roi_idx=roi_idx-1;
        end
        roi_hJTable.changeSelection(roi_idx,-1, false, false);
    end
    
    
    %% Callbacks -- File Panel
    function load_CB(~, ~)
        curr_hsdataPath=get(handles.file.hsdataPath,'String');
        if isempty(curr_hsdataPath)
            curr_hsdataPath=lastPath;
        end
        [filename, pathname]=uigetfile(obj.LoadableDataFormats, 'Load A Data Source', curr_hsdataPath);
        if ~ischar(filename) || ~ischar(pathname)
            return
        end
        try
            obj.load(fullfile(pathname,filename));
        catch err
            disp(getReport(err));
            errordlg(err.message,err.identifier);
        end
        setROITableData();
        positionAxes();
        updateControls();
    end

    function save_CB(~,~)
        if ~obj.initialized
            errordlg('Object is not yet Initialized');
            return
        end
        curr_dataPath=get(handles.file.hsdataPath,'String');
        if ~isempty(obj.saveFilePath) && ~strcmp(curr_dataPath,obj.saveFilePath)
            obj.saveas(curr_dataPath);
        else
            obj.save();
        end
        obj.dirty=false;
        updateControls();
    end

    function reset_CB(~,~)
        lastPath=obj.workingDir;
        obj.reset();
        setROITableData();
        updateControls();
        cla(handles.axes);
        set(handles.axes,'Visible','off');
    end

    
    function close_CB(~,~)
        try
            obj.clearJTableCallbacks(handles.ROITable);
        catch err
            disp(getReport(err))
        end
        obj.closeGUI();
    end

    %% Callbacks -- Physical Panel

    function setPSF_CB(~,~)
        psf=str2num(get(gcbo(),'String')); %#ok<*ST2NM>
        if isscalar(psf)
            if ~isfinite(psf) || psf<=0
               warning('HSData:gui','Invalid psf: %g',psf);
            else
                obj.psf=psf;
                obj.dirty=true;
            end
        elseif length(psf)==2
            if any(~isfinite(psf)) || any(psf<=0)
               warning('HSData:gui','Invalid psf: %s',mat2str(psf));
            else
                obj.psf=psf;
                obj.dirty=true;
            end
        else
            warning('HSData:gui','Invalid psf: "%s"',get(gcbo(),'String'));
        end
        updateControls();
    end

    function setPixelSize_CB(~,~)
        pixelsize=str2num(get(gcbo(),'String'));
        if isscalar(pixelsize)
            if ~isfinite(pixelsize) || pixelsize<=0
               warning('HSData:gui','Invalid pixel size: %g',pixelsize);
            else
                obj.pixelSize=pixelsize;
                obj.dirty=true;
            end
        elseif length(pixelsize)==2
            if any(~isfinite(pixelsize)) || any(pixelsize<=0)
               warning('HSData:gui','Invalid pixel size: %s',mat2str(pixelsize));
            else
                obj.pixelSize=pixelsize;
                obj.dirty=true;
            end
        else
            warning('HSData:gui','Invalid pixel size: "%s"',get(gcbo(),'String'));
        end
        updateControls();
    end

%% Callbacks -- ROI Panel
   
    function addROI_CB(~,~)
        rect=getrect(handles.axes.sum);
        roi=round([rect(1), rect(1)+rect(3), rect(2), rect(2)+rect(4)]);
        roi(1)=max(1,roi(1));
        roi(2)=min(obj.sizeX,roi(2));
        roi(3)=max(1,roi(3));
        roi(4)=min(obj.sizeY,roi(4));
        obj.addROI(roi);
        setROITableData();
        setSelectedROI(length(obj.ROI));
        updateControls();
    end

    function selectSpectra_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx==0
            return
        end
        rect=getrect(handles.axes.spectra);
        roi=obj.ROI{roi_idx};
        roi(6)=find(obj.lambda<=max(obj.lambda(end),rect(1)),1,'first');
        roi(5)=find(obj.lambda>=min(obj.lambda(1), rect(1)+rect(3)),1,'last');
        obj.modifyROI(roi_idx,roi);
        setROITableData();
        setSelectedROI(roi_idx);
        updateControls();
    end

    function delROI_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0
            obj.deleteROI(roi_idx);
            setROITableData();        
            setSelectedROI(min(length(obj.ROI),roi_idx));
            updateControls();
        end
    end

    function clearROI_CB(~,~)
        if strcmp('Yes',questdlg('Really clear all ROI?'))
            obj.clearROI();
            setROITableData();
            updateControls();
        end
    end

    function trackHSRPT_ROI_CB(~,~)
        rpt=obj.trackHSRPT(getSelectedROI());
        rpt.gui();
        name=saveGlobalVar('hsrpt',rpt);
        fprintf('Saved RPT object as workspace variable "%s"\n',name);
    end

    function viewROI_RGB_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewFramesRGB(roi_idx);
        end
    end

    function viewROI_Sum_RGB_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewSumImageRGB(roi_idx);
        end
    end

    function viewROI_Sum_Gray_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewSumImageGray(roi_idx);
        end
    end

    function viewROI_4D_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewFrames(roi_idx);
        end
    end

    function viewROI_Sum_Slices_Intensity_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewSurfaceSliceIntensity(0,roi_idx);
        end
    end
    function viewROI_Sum_Slices_Wavelength_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewSurfaceSliceWavelength(0,roi_idx);
        end
    end
    function viewROI_Sum_Volumes_Intensity_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewVolumeIntensity(0,roi_idx);
        end
    end
    function viewROI_Frames_Volumes_Intensity_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            roi=obj.ROI{roi_idx};
            obj.viewVolumeIntensity(roi(7),roi_idx);
        end
    end
    function viewROI_Frames_Volumes_Wavelength_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            roi=obj.ROI{roi_idx};
            obj.viewVolumeWavelength(roi(7),roi_idx);
        end
    end
    function viewROI_Sum_Volumes_Wavelength_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewVolumeWavelength(0,roi_idx);
        end
    end

    function viewROI_Spectra_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewSpectra(roi_idx);
        end
    end

    function viewROI_SpectralSeries_CB(~,~)
        roi_idx=getSelectedROI();
        if roi_idx>0 && roi_idx<=length(obj.ROI)
            obj.viewSpectralSeries(roi_idx);
        end
    end


    function updateROISelection_CB(~,~)
        displayROI();
    end

    function ROIEdit_CB(~,ev)
        roi_idx=ev.Indices(1);
        edit_col=ev.Indices(2);
        if roi_idx>=0 && roi_idx<=length(obj.ROI)
            if edit_col==1
                roi_name=ev.NewData;
                roi=obj.ROI{roi_idx};
            else
                roi_name=obj.ROIname{roi_idx};
                roi=obj.ROI{roi_idx};
                roi(edit_col-1)=ev.NewData;
            end
            try
                obj.modifyROI(roi_idx, roi, roi_name);
            catch
            end
        end
        setSelectedROI(roi_idx);
        updateControls();
    end

    function roiContextMenu_CB(hObj,~)
        nROI=length(obj.ROI);
        cp=get(handles.axes.sum,'CurrentPoint');
        cp=cp(1,1:2);
        selROI=getSelectedROI();
        for N=1:nROI
            n=mod(N-1+selROI-1,nROI)+1; % n will iterate through all ROI starting at selected
            roi=obj.ROI{n};
            if cp(1)>=roi(1) && cp(1)<=roi(2) && cp(2)>=roi(3) && cp(2)<=roi(4)
                setSelectedROI(n);
                updateControls();
                return
            end
        end
    end

    function roiButtonDown_CB(hObj, ~)
        nROI=length(obj.ROI);
        if nROI<2
            return
        end
        select_type=get(gcbf(),'SelectionType');
        switch select_type
            case {'normal','open'}
                cp=get(get(hObj,'Parent'),'CurrentPoint');
                cp=cp(1,1:2);
                selROI=getSelectedROI();
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

    function rotateSpectra_CB(~,~)
        rotate3d();
%         set(h,'ButtonDownFilter', @(~,~) false);
    end
end



