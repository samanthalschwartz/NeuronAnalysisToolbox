% to still do/figure out:
% set cell mask cleanup button
% set save to excel button
% -- need to properly set the load function for ashley files to fix memory
% issues of saving things as dipimages
function guiFig = gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end
%--- create structures for all guis associated with this gui----
    app = struct();
    imloader = struct();
    Masking = struct();
    % --- global vars for surf sign accum
    filenames = [];
    updatefilenames();
    maskchangebool = 0;
    savenamechangebool = 0;
    % --- global vars related to imloader
    defaultdropdowntxt = 'Load Open Image';
    %------------------------------------
%--- now actually make the base gui    
    createComponents();
    obj.guiFig = app.SurfaceSignalAccumulationBase_UIFigure;
    %----
    % -------------  gui makers -------------------
    function createComponents()  
        % Create SurfaceSignalAccumulationBase_UIFigure and hide until all components are created
            app.SurfaceSignalAccumulationBase_UIFigure = uifigure('Visible', 'off');
            app.SurfaceSignalAccumulationBase_UIFigure.Position = [100 100 715 705];
            app.SurfaceSignalAccumulationBase_UIFigure.Name = 'Surface Signal Accumulation';

            % Create MethodsPanel
            app.MethodsPanel = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
            app.MethodsPanel.Title = 'Methods';
            app.MethodsPanel.Position = [553 567 148 126];

            % Create ImageLoaderButton
            app.ImageLoaderButton = uibutton(app.MethodsPanel, 'push');
            app.ImageLoaderButton.Position = [21 72 100 27];
            app.ImageLoaderButton.Text = 'Image Loader';

            % Create MaskingButton
            app.MaskingButton = uibutton(app.MethodsPanel, 'push');
            app.MaskingButton.Position = [21 42 100 26];
            app.MaskingButton.Text = 'Masking';

            % Create SelectRegionsButton
            app.SelectRegionsButton = uibutton(app.MethodsPanel, 'push');
            app.SelectRegionsButton.Position = [21 9 100 28];
            app.SelectRegionsButton.Text = 'Select Regions';

            % Create ViewerPanel
            app.ViewerPanel = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
            app.ViewerPanel.Title = 'Viewer';
            app.ViewerPanel.Position = [15 377 526 316];

            % Create Tree
            app.Tree = uitree(app.ViewerPanel);
            app.Tree.Multiselect = 'on';
            app.Tree.Position = [16 72 157 214];

            % Create CellMarkerNode
            app.CellMarkerNode = uitreenode(app.Tree);
            app.CellMarkerNode.Text = 'Cell Marker';

            % Create CellMarker_ImageNode
            app.CellMarker_ImageNode = uitreenode(app.CellMarkerNode);
            app.CellMarker_ImageNode.Text = 'Cell Image';

            % Create CellMarker_MaskNode
            app.CellMarker_MaskNode = uitreenode(app.CellMarkerNode);
            app.CellMarker_MaskNode.Text = 'Cell Mask';

            % Create SurfaceSignalNode
            app.SurfaceSignalNode = uitreenode(app.Tree);
            app.SurfaceSignalNode.Text = 'Surface Signal';

            % Create SurfSign_ImageNode
            app.SurfSign_ImageNode = uitreenode(app.SurfaceSignalNode);
            app.SurfSign_ImageNode.Text = 'Surface Signal Image';

            % Create SurfSign_MaskNode
            app.SurfSign_MaskNode = uitreenode(app.SurfaceSignalNode);
            app.SurfSign_MaskNode.Text = 'Surface Signal Mask';

            % Create CellMarker_DistanceMaskNode
            app.CellMarker_DistanceMaskNode = uitreenode(app.Tree);
            app.CellMarker_DistanceMaskNode.Text = 'Distance Mask';

            % Create AISRegionNode
            app.AISRegionNode = uitreenode(app.Tree);
            app.AISRegionNode.Text = 'AIS Region';

            % Create SomaRegionNode
            app.SomaRegionNode = uitreenode(app.Tree);
            app.SomaRegionNode.Text = 'Soma Region';

            % Create ViewSelectedButton
            app.ViewSelectedButton = uibutton(app.ViewerPanel, 'push');
            app.ViewSelectedButton.Position = [45 36 100 22];
            app.ViewSelectedButton.Text = 'View Selected';

            % Create UIAxes_ImageViewer
            app.UIAxes_ImageViewer = uiaxes(app.ViewerPanel);
            title(app.UIAxes_ImageViewer, '')
            xlabel(app.UIAxes_ImageViewer, '')
            ylabel(app.UIAxes_ImageViewer, '')
            app.UIAxes_ImageViewer.Position = [176 10 333 284];

            % Create SaveNameEditField
            app.SaveNameEditField = uieditfield(app.SurfaceSignalAccumulationBase_UIFigure, 'text');
            app.SaveNameEditField.Position = [98 331 460 22];

            % Create SaveNameButton
            app.SaveNameButton = uibutton(app.SurfaceSignalAccumulationBase_UIFigure, 'push');
            app.SaveNameButton.Position = [15 331 78 22];
            app.SaveNameButton.Text = 'Save Name';

            % Create Panel_FrameRates
            app.Panel_FrameRates = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
            app.Panel_FrameRates.Title = 'Frame Rate (min/frame)';
            app.Panel_FrameRates.Position = [553 482 148 75];

            % Create PreReleaseEditFieldLabel
            app.PreReleaseEditFieldLabel = uilabel(app.Panel_FrameRates);
            app.PreReleaseEditFieldLabel.HorizontalAlignment = 'right';
            app.PreReleaseEditFieldLabel.Position = [23 30 75 22];
            app.PreReleaseEditFieldLabel.Text = 'Pre-Release:';

            % Create PreReleaseEditField
            app.PreReleaseEditField = uieditfield(app.Panel_FrameRates, 'numeric');
            app.PreReleaseEditField.Position = [103 30 27 22];

            % Create PostReleaseEditFieldLabel
            app.PostReleaseEditFieldLabel = uilabel(app.Panel_FrameRates);
            app.PostReleaseEditFieldLabel.HorizontalAlignment = 'right';
            app.PostReleaseEditFieldLabel.Position = [15 4 81 22];
            app.PostReleaseEditFieldLabel.Text = 'Post-Release:';

            % Create PostReleaseEditField
            app.PostReleaseEditField = uieditfield(app.Panel_FrameRates, 'numeric');
            app.PostReleaseEditField.Position = [103 4 27 22];

            % Create Panel_ReleaseFrame
            app.Panel_ReleaseFrame = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
            app.Panel_ReleaseFrame.Position = [553 445 148 33];

            % Create ReleaseFrameEditFieldLabel
            app.ReleaseFrameEditFieldLabel = uilabel(app.Panel_ReleaseFrame);
            app.ReleaseFrameEditFieldLabel.HorizontalAlignment = 'right';
            app.ReleaseFrameEditFieldLabel.Position = [5 5 91 22];
            app.ReleaseFrameEditFieldLabel.Text = 'Release Frame:';

            % Create ReleaseFrameEditField
            app.ReleaseFrameEditField = uieditfield(app.Panel_ReleaseFrame, 'numeric');
            app.ReleaseFrameEditField.Position = [103 5 28 22];

            % Create SAVEButton
            app.SAVEButton = uibutton(app.SurfaceSignalAccumulationBase_UIFigure, 'push');
            app.SAVEButton.Position = [574 325 118 35];
            app.SAVEButton.Text = 'SAVE';

            % Create PlottingPanel
            app.PlottingPanel = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
            app.PlottingPanel.Title = 'Plotting';
            app.PlottingPanel.Position = [15 13 686 300];

            % Create NormalizationLabel
            app.NormalizationLabel = uilabel(app.PlottingPanel);
            app.NormalizationLabel.HorizontalAlignment = 'right';
            app.NormalizationLabel.Position = [5 214 78 22];
            app.NormalizationLabel.Text = 'Normalization:';

            % Create NormalizationDropDown
            app.NormalizationDropDown = uidropdown(app.PlottingPanel);
            app.NormalizationDropDown.Position = [99 214 123 22];

            % Create ImageRegionLabel
            app.ImageRegionLabel = uilabel(app.PlottingPanel);
            app.ImageRegionLabel.HorizontalAlignment = 'right';
            app.ImageRegionLabel.Position = [1 183 83 22];
            app.ImageRegionLabel.Text = 'Image Region:';

            % Create ImageRegionDropDown
            app.ImageRegionDropDown = uidropdown(app.PlottingPanel);
            app.ImageRegionDropDown.Position = [99 183 123 22];

            % Create EndFramefornormalizationEditFieldLabel
            app.EndFramefornormalizationEditFieldLabel = uilabel(app.PlottingPanel);
            app.EndFramefornormalizationEditFieldLabel.HorizontalAlignment = 'right';
            app.EndFramefornormalizationEditFieldLabel.Position = [9 150 168 22];
            app.EndFramefornormalizationEditFieldLabel.Text = 'End Frame (for normalization):';

            % Create EndFramefornormalizationEditField
            app.EndFramefornormalizationEditField = uieditfield(app.PlottingPanel, 'numeric');
            app.EndFramefornormalizationEditField.Position = [187 150 27 22];

            % Create UIAxes_ResultsPlot
            app.UIAxes_ResultsPlot = uiaxes(app.PlottingPanel);
            title(app.UIAxes_ResultsPlot, 'Post-Release Surface Accumulation')
            xlabel(app.UIAxes_ResultsPlot, 'Time (mins)')
            ylabel(app.UIAxes_ResultsPlot, 'Relative Intensity Density (AU)')
            app.UIAxes_ResultsPlot.Position = [229 10 448 240];
            
             % Create TemporalHeatMapButton
            app.TemporalHeatMapButton = uibutton(app.PlottingPanel, 'push');
            app.TemporalHeatMapButton.Position = [65 104 118 30];
            app.TemporalHeatMapButton.Text = 'Temporal HeatMap';

            % Create ViewResultsArrayButton
            app.ViewResultsArrayButton = uibutton(app.PlottingPanel, 'push');
            app.ViewResultsArrayButton.Position = [65 67 117 29];
            app.ViewResultsArrayButton.Text = 'View Results Array';
            
            % Create RawResultsVal_CheckBox
            app.RawResultsVal_CheckBox = uicheckbox(app.PlottingPanel);
            app.RawResultsVal_CheckBox.Text = '';
            app.RawResultsVal_CheckBox.Position = [26 66 25 22];

            % Create RawResultsValRawValsLabel
            app.RawResultsValRawValsLabel = uilabel(app.PlottingPanel);
            app.RawResultsValRawValsLabel.Position = [8 83 55 22];
            app.RawResultsValRawValsLabel.Text = 'Raw Vals';

            % Create SavetoExcelButton
            app.SavetoExcelButton = uibutton(app.PlottingPanel, 'push');
            app.SavetoExcelButton.Position = [66 29 118 29];
            app.SavetoExcelButton.Text = 'Save to Excel';

            % Create DistancespixelsLabel
            app.DistancespixelsLabel = uilabel(app.PlottingPanel);
            app.DistancespixelsLabel.HorizontalAlignment = 'right';
            app.DistancespixelsLabel.Position = [1 249 103 22];
            app.DistancespixelsLabel.Text = 'Distances (um):';

            % Create DistancespixelsEditField
            app.DistancespixelsEditField = uieditfield(app.PlottingPanel, 'text');
            app.DistancespixelsEditField.Position = [119 249 124 22];
            
            % Create ViewDistancesButton
            app.ViewDistancesButton = uibutton(app.PlottingPanel, 'push');
            app.ViewDistancesButton.Position = [252 249 100 22];
            app.ViewDistancesButton.Text = 'View Distances';

            % Create Panel_ReleaseFrame_2
            app.Panel_ReleaseFrame_2 = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
            app.Panel_ReleaseFrame_2.Position = [553 380 148 61];

            % Create PixelSizeumEditFieldLabel
            app.PixelSizeumEditFieldLabel = uilabel(app.Panel_ReleaseFrame_2);
            app.PixelSizeumEditFieldLabel.HorizontalAlignment = 'right';
            app.PixelSizeumEditFieldLabel.Position = [0 33 90 22];
            app.PixelSizeumEditFieldLabel.Text = 'Pixel Size (um):';

            % Create PixelSizeumEditField
            app.PixelSizeumEditField = uieditfield(app.Panel_ReleaseFrame_2, 'numeric');
            app.PixelSizeumEditField.Position = [95 33 44 22];

            % Create PixelSizeDropDown
            app.PixelSizeDropDown = uidropdown(app.Panel_ReleaseFrame_2);
            app.PixelSizeDropDown.Position = [8 7 131 22];
            
            % Create UpdatePlotButton
            app.UpdatePlotButton = uibutton(app.PlottingPanel, 'push');
            app.UpdatePlotButton.Enable = 'off';
            app.UpdatePlotButton.Position = [577 249 100 22];
            app.UpdatePlotButton.Text = 'Update Plot';

            % Show the figure after all components are created
            app.SurfaceSignalAccumulationBase_UIFigure.Visible = 'on';
        
        %------ setting callbacks etc -------
        app.SurfaceSignalAccumulationBase_UIFigure.WindowButtonMotionFcn = @(object, eventdata) callback_CheckPlotButtonState();
        expand(app.Tree);
        app.ImageLoaderButton.ButtonPushedFcn = @(btn,event) callback_ImageLoaderButton();
        app.ViewSelectedButton.ButtonPushedFcn = @(btn,event) callback_ViewSelectedButton();
        app.SelectRegionsButton.ButtonPushedFcn= @(btn,event) callback_SelectRegionsButton();
        makeTempHeatMapNode(); % if there is a heatmap then add the node - otherwise it will be added when made with the temporal heatmap button
        app.Tree.SelectionChangedFcn =  @(btn,event) callback_TreeNodeChange();
        app.Tree.Interruptible = 'on';
        app.MaskingButton.ButtonPushedFcn = @(btn,event) callback_MaskingButton();
        app.PixelSizeDropDown.Items = {'Bin 2x2 = 0.228 um','Bin 1x1 = 0.114 um'};
        initializePixelSize();
        app.PixelSizeumEditField.ValueChangedFcn = @(btn,event) callback_PixelFieldChange();
        app.PixelSizeDropDown.ValueChangedFcn = @(btn,event) callback_PixelSizeDropdownChange();
        app.DistancespixelsEditField.Value = '4,10,200';
        app.DistancespixelsEditField.ValueChangedFcn =  @(btn,event) callback_DistancesChanged();
        app.ViewDistancesButton.ButtonPushedFcn = @(btn,event) callback_ViewDistancesButton();
        app.NormalizationDropDown.Items = {'Per Distance Region','Compare to First Region'};
        app.NormalizationDropDown.ValueChangedFcn = @(btn,event) makePlot();
        app.ImageRegionDropDown.ValueChangedFcn = @(btn,event) makePlot();
        app.EndFramefornormalizationEditField.ValueChangedFcn = @(btn,event) makePlot(); 
        app.EndFramefornormalizationEditField.RoundFractionalValues = 'on';
        initializePlottingComponents();
        initializeImagingParams();
        app.UpdatePlotButton.ButtonPushedFcn = @(btn,event) callback_UpdatePlot();
        app.PreReleaseEditField.ValueChangedFcn = @(btn,event) setImagingParams(); 
        app.PostReleaseEditField.ValueChangedFcn = @(btn,event) setImagingParams(); 
        app.ReleaseFrameEditField.ValueChangedFcn = @(btn,event) setImagingParams(); 
        app.ViewResultsArrayButton.ButtonPushedFcn = @(btn,event) makePlot(1,0);
        app.SavetoExcelButton.ButtonPushedFcn = @(btn,event) makePlot(0,1);
        app.TemporalHeatMapButton.ButtonPushedFcn = @(btn,event) callback_TempHeatMapButton();
        app.SAVEButton.ButtonPushedFcn = @(btn,event) callback_saveAshleyFile();
        app.SaveNameButton.ButtonPushedFcn = @(btn,event) callback_selectSaveName();
        app.SaveNameEditField.ValueChangedFcn = @(btn,event) setUpdatedSaveName();
        app.SAVEButton.BackgroundColor = [1 .7 .7];
        makePlot();
        if ~isempty(obj.savename)
            app.SaveNameEditField.Value = obj.savename;
        end
            % start with cell fill image - but if there's nothing than have
            % blank axes (if started gui with empty object)
        app.Tree.SelectedNodes = app.CellMarker_ImageNode;
        if ~isempty(obj.cellFill.image)
            callback_TreeNodeChange();
        end
    end
    function createImageLoader()
        updatefilenames();
        imloader.ImageLoader = uifigure('WindowButtonDownFcn', @(src,event)callback_updateFileDropdown);
        imloader.ImageLoader.Position = [100 100 390 275];
        imloader.ImageLoader.Name = 'Image Loader';
        
        % Create TabGroup
        imloader.TabGroup = uitabgroup(imloader.ImageLoader,'SelectionChangedFcn', @(src,event)callback_updateFileDropdown);
        imloader.TabGroup.Position = [18 141 357 124];
        
        % Create CellMarkerTab
        imloader.CellMarkerTab = uitab(imloader.TabGroup,'ButtonDownFcn', @(src,event)callback_updateFileDropdown);
        imloader.CellMarkerTab.Title = 'Cell Marker';
        
        % Create CellMarker_ImageNameEditFieldLabel
        imloader.CellMarker.EditFieldLabel = uilabel(imloader.CellMarkerTab);
        imloader.CellMarker.EditFieldLabel.HorizontalAlignment = 'right';
        imloader.CellMarker.EditFieldLabel.Position = [5 39 56 22];
        imloader.CellMarker.EditFieldLabel.Text = 'Image:';
        
        % Create CellMarker_ImageNameEditField
        imloader.CellMarker.ImageNameEditField = uieditfield(imloader.CellMarkerTab, 'text');
        imloader.CellMarker.ImageNameEditField.Position = [76 41 272 18];
        imloader.CellMarker.ImageNameEditField.Editable = 'off';
        
        % Create CellMarker_LoadImageDropDown
        imloader.CellMarker.LoadImageDropDown = uidropdown(imloader.CellMarkerTab,'ValueChangedFcn',@(btn,event) callback_CellMarker_LoadImageDropDown());
        imloader.CellMarker.LoadImageDropDown.Items = [{defaultdropdowntxt},filenames];
        imloader.CellMarker.LoadImageDropDown.Position = [197 64 151 21];
        imloader.CellMarker.LoadImageDropDown.Value = defaultdropdowntxt;
        
        % Create CellMarker_ORtext
        imloader.CellMarker_ORtext = uilabel(imloader.CellMarkerTab);
        imloader.CellMarker_ORtext.Position = [175 63 23 22];
        imloader.CellMarker_ORtext.Text = 'or';
        
        % Create CellMarker_LoadImageFromFileButton
        imloader.CellMarker.LoadImageFromFileButton = uibutton(imloader.CellMarkerTab, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_CellMarker_LoadImage());
        imloader.CellMarker.LoadImageFromFileButton.Position = [13 63 149 22];
        imloader.CellMarker.LoadImageFromFileButton.Text = 'Load Image From File';
        
        % Create CellMarker_ViewImageButton
        imloader.CellMarker.ViewImageButton = uibutton(imloader.CellMarkerTab, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_CellMarker_ViewImage());
        imloader.CellMarker.ViewImageButton.Position = [121 11 99 22];
        imloader.CellMarker.ViewImageButton.Text = 'View Image';
        imloader.CellMarker.ViewImageButton.Enable = 'off';
        
        % Create ChannelSpinnerLabel
        imloader.CellMarker.ChannelSpinnerLabel = uilabel(imloader.CellMarkerTab);
        imloader.CellMarker.ChannelSpinnerLabel.HorizontalAlignment = 'right';
        imloader.CellMarker.ChannelSpinnerLabel.Enable = 'off';
        imloader.CellMarker.ChannelSpinnerLabel.Position = [13 11 45 22];
        imloader.CellMarker.ChannelSpinnerLabel.Text = 'Channel';
        
        % Create CellMarker_ChannelSpinner
        imloader.CellMarker.ChannelSpinner = uispinner(imloader.CellMarkerTab);
        imloader.CellMarker.ChannelSpinner.Enable = 'off';
        imloader.CellMarker.ChannelSpinner.Position = [66 11 39 22];
        
        % Create SurfaceSignalTab
        imloader.SurfaceSignalTab = uitab(imloader.TabGroup,'ButtonDownFcn', @(src,event)callback_updateFileDropdown);
        imloader.SurfaceSignalTab.Title = 'Surface Signal';
         
         % Create SurfSign_ImageNameEditFieldLabel
        imloader.SurfSign.EditFieldLabel = uilabel(imloader.SurfaceSignalTab);
        imloader.SurfSign.EditFieldLabel.HorizontalAlignment = 'right';
        imloader.SurfSign.EditFieldLabel.Position = [5 39 56 22];
        imloader.SurfSign.EditFieldLabel.Text = 'Image:';
        
        % Create SurfSign_ImageNameEditField
        imloader.SurfSign.ImageNameEditField = uieditfield(imloader.SurfaceSignalTab, 'text');
        imloader.SurfSign.ImageNameEditField.Position = [76 41 272 18];
        imloader.SurfSign.ImageNameEditField.Editable = 'off';
        
        % Create SurfSign_LoadImageDropDown
        imloader.SurfSign.LoadImageDropDown = uidropdown(imloader.SurfaceSignalTab,'ValueChangedFcn',@(btn,event) callback_SurfaceSignal_LoadImageDropDown());
        imloader.SurfSign.LoadImageDropDown.Items = [{defaultdropdowntxt},filenames];
        imloader.SurfSign.LoadImageDropDown.Position = [197 64 151 21];
        imloader.SurfSign.LoadImageDropDown.Value = defaultdropdowntxt;
        
        % Create SurfSign_ORtext
        imloader.SurfSign_ORtext = uilabel(imloader.SurfaceSignalTab);
        imloader.SurfSign_ORtext.Position = [175 63 23 22];
        imloader.SurfSign_ORtext.Text = 'or';
        
        % Create SurfSign_LoadImageFromFileButton
        imloader.SurfSign.LoadImageFromFileButton = uibutton(imloader.SurfaceSignalTab, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_SurfSignal_LoadImage());
        imloader.SurfSign.LoadImageFromFileButton.Position = [13 63 149 22];
        imloader.SurfSign.LoadImageFromFileButton.Text = 'Load Image From File';
        
        % Create SurfSign_ViewImageButton
        imloader.SurfSign.ViewImageButton = uibutton(imloader.SurfaceSignalTab, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_SurfSignal_ViewImage());
        imloader.SurfSign.ViewImageButton.Position = [121 11 99 22];
        imloader.SurfSign.ViewImageButton.Text = 'View Image';
        imloader.SurfSign.ViewImageButton.Enable = 'off';
        
        % Create SurfSign ChannelSpinner_Label
        imloader.SurfSign.ChannelSpinnerLabel = uilabel(imloader.SurfaceSignalTab);
        imloader.SurfSign.ChannelSpinnerLabel.HorizontalAlignment = 'right';
        imloader.SurfSign.ChannelSpinnerLabel.Enable = 'off';
        imloader.SurfSign.ChannelSpinnerLabel.Position = [13 11 45 22];
        imloader.SurfSign.ChannelSpinnerLabel.Text = 'Channel';
        
        % Create SurfSign_ChannelSpinner
        imloader.SurfSign.ChannelSpinner = uispinner(imloader.SurfaceSignalTab);
        imloader.SurfSign.ChannelSpinner.Enable = 'off';
        imloader.SurfSign.ChannelSpinner.Position = [66 11 39 22];
        
        % Create SetPanel
        imloader.SetPanel = uipanel(imloader.ImageLoader);
        imloader.SetPanel.Position = [18 71 357 71];
        
        % Create DriftCorrectImageCheckBox
        imloader.DriftCorrectImageCheckBox = uicheckbox(imloader.SetPanel,'ValueChangedFcn',@(btn,event) callback_DriftCorrectBoxChanged());
        imloader.DriftCorrectImageCheckBox.Text = 'Drift Correct Images';
        imloader.DriftCorrectImageCheckBox.Position = [150 30 124 23];
        imloader.DriftCorrectImageCheckBox.Enable = 'off';
        
        % Create DrifCorrText
        imloader.DrifCorrText = uilabel(imloader.SetPanel);
        imloader.DrifCorrText.Position = [150 14 135 22];
        imloader.DrifCorrText.Text = '(using Cell Marker drift)';
        imloader.DrifCorrText.Enable = 'off';
        % Create SetImagesButton
        imloader.SetImagesButton = uibutton(imloader.SetPanel, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_SetImageButton());
        imloader.SetImagesButton.Position = [12 14 132 40];
        imloader.SetImagesButton.Text = 'Set Images';
        imloader.SetImagesButton.Enable = 'off';
        
        % Create CropPanel
        imloader.CropPanel = uipanel(imloader.ImageLoader);
        imloader.CropPanel.Position = [18 16 357 56];
        
        % Create CropActiveImageRegionButton
        imloader.CropActiveImageRegionButton = uibutton(imloader.CropPanel, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_CropImage());
        imloader.CropActiveImageRegionButton.Position = [12 11 149 35];
        imloader.CropActiveImageRegionButton.Text = 'Crop Active Image Region';
        
        % Create ResetButton
        imloader.ResetButton = uibutton(imloader.CropPanel, 'state');
        imloader.ResetButton.Text = 'Reset';
        imloader.ResetButton.Position = [305 11 43 28];
        if obj.cellFill.ROI_trim
            imloader.ResetButton.Value = true; %if a trim region has been set than have this option
        else
            imloader.ResetButton.Value = false;
        end
        %----------------- setup callbacks
        imloader.ResetButton.ValueChangedFcn =  @(src,event) callback_ResetButton();
        
    end
    function createMasking()
        % Create UIFigure and hide until all components are created
        Masking.UIFigure = uifigure('Visible', 'off');
        Masking.UIFigure.Position = [100 100 305 331];
        Masking.UIFigure.Name = 'Masking';
        
        % Create TabGroup
        Masking.TabGroup = uitabgroup(Masking.UIFigure);
        Masking.TabGroup.Position = [22 21 258 294];
        
        % Create Cell_Tab
        Masking.Cell_Tab = uitab(Masking.TabGroup);
        Masking.Cell_Tab.Title = 'Cell Mask';
        
        % Create Cell_MaskPanel
        Masking.Cell_MaskPanel = uipanel(Masking.Cell_Tab);
        Masking.Cell_MaskPanel.Position = [20 158 218 99];
        
        % Create Cell_MaskButton
        Masking.Cell_MaskButton = uibutton(Masking.Cell_MaskPanel, 'push');
        Masking.Cell_MaskButton.Position = [42 17 130 42];
        Masking.Cell_MaskButton.Text = 'Mask Image';
        
        % Create MaskingMethodDropDown
        Masking.MaskingMethodDropDownLabel = uilabel(Masking.Cell_MaskPanel);
        Masking.MaskingMethodDropDownLabel.HorizontalAlignment = 'right';
        Masking.MaskingMethodDropDownLabel.Position = [6 67 94 22];
        Masking.MaskingMethodDropDownLabel.Text = 'Masking Method';
        
        % Create Cell_MaskDropDown
        Masking.Cell_MaskDropDown = uidropdown(Masking.Cell_MaskPanel);
        Masking.Cell_MaskDropDown.Position = [107 67 100 22];
        
        % Create Cell_CleanPanel
        Masking.Cell_CleanPanel = uipanel(Masking.Cell_Tab);
        Masking.Cell_CleanPanel.Position = [20 40 218 99];
        
        % Create CleanUpMethodDropDown_2Label
        Masking.CleanUpMethodDropDownLabel = uilabel(Masking.Cell_CleanPanel);
        Masking.CleanUpMethodDropDownLabel.HorizontalAlignment = 'right';
        Masking.CleanUpMethodDropDownLabel.Position = [6 69 99 22];
        Masking.CleanUpMethodDropDownLabel.Text = 'Clean Up-Method';
        
        % Create Cell_CleanUpMethodDropDown
        Masking.Cell_CleanUpMethodDropDown = uidropdown(Masking.Cell_CleanPanel);
        Masking.Cell_CleanUpMethodDropDown.Position = [112 69 100 22];
        
        % Create Cell_CleanUpMaskButton
        Masking.Cell_CleanUpMaskButton = uibutton(Masking.Cell_CleanPanel, 'push');
        Masking.Cell_CleanUpMaskButton.Position = [44 14 130 44];
        Masking.Cell_CleanUpMaskButton.Text = 'Clean Up Mask';
        
        % Create Cell_ResetCheckBox
        Masking.Cell_ResetCheckBox = uicheckbox(Masking.Cell_CleanPanel);
        Masking.Cell_ResetCheckBox.Text = '';
        Masking.Cell_ResetCheckBox.Position = [16 20 25 22];
        
        % Create Cell_ResetLabel
        Masking.Cell_ResetLabel = uilabel(Masking.Cell_CleanPanel);
        Masking.Cell_ResetLabel.Position = [6 36 37 22];
        Masking.Cell_ResetLabel.Text = 'Reset';
        
        % Create Surface_Tab
        Masking.Surface_Tab = uitab(Masking.TabGroup);
        Masking.Surface_Tab.Title = 'Surface Signal Mask';
        
        
        % Create Surface_MaskPanel
        Masking.Surface_MaskPanel = uipanel(Masking.Surface_Tab);
        Masking.Surface_MaskPanel.Position = [20 158 218 99];
        
        % Create Surface_MaskButton
        Masking.Surface_MaskButton = uibutton(Masking.Surface_MaskPanel, 'push');
        Masking.Surface_MaskButton.Position = [34 14 130 42];
        Masking.Surface_MaskButton.Text = 'Mask Image';
        
        % Create MaskingMethodDropDownLabel
        Masking.MaskingMethodDropDownLabel = uilabel(Masking.Surface_MaskPanel);
        Masking.MaskingMethodDropDownLabel.HorizontalAlignment = 'right';
        Masking.MaskingMethodDropDownLabel.Position = [6 67 94 22];
        Masking.MaskingMethodDropDownLabel.Text = 'Masking Method';
        
        % Create Surface_MaskDropDown
        Masking.Surface_MaskDropDown = uidropdown(Masking.Surface_MaskPanel);
        Masking.Surface_MaskDropDown.Position = [107 67 100 22];
        
        % Create Surface_CleanPanel
        Masking.Surface_CleanPanel = uipanel(Masking.Surface_Tab);
        Masking.Surface_CleanPanel.Position = [20 10 218 139];
        
        % Create CleanUpMethodDropDownLabel
        Masking.CleanUpMethodDropDownLabel = uilabel(Masking.Surface_CleanPanel);
        Masking.CleanUpMethodDropDownLabel.HorizontalAlignment = 'right';
        Masking.CleanUpMethodDropDownLabel.Position = [6 109 99 22];
        Masking.CleanUpMethodDropDownLabel.Text = 'Clean Up-Method';
        
        % Create Surface_CleanUpMethodDropDown
        Masking.Surface_CleanUpMethodDropDown = uidropdown(Masking.Surface_CleanPanel);
        Masking.Surface_CleanUpMethodDropDown.Position = [112 109 100 22];
        
        % Create Surface_CleanUpMaskButton
        Masking.Surface_CleanUpMaskButton = uibutton(Masking.Surface_CleanPanel, 'push');
        Masking.Surface_CleanUpMaskButton.Position = [34 50 130 44];
        Masking.Surface_CleanUpMaskButton.Text = 'Clean Up Mask';
        
        % Create Surface_ResetCheckBox
        Masking.Surface_ResetCheckBox = uicheckbox(Masking.Surface_CleanPanel);
        Masking.Surface_ResetCheckBox.Text = '';
        Masking.Surface_ResetCheckBox.Position = [48 9 25 21];
        
        % Create Surface_ResetLabel
        Masking.Surface_ResetLabel = uilabel(Masking.Surface_CleanPanel);
        Masking.Surface_ResetLabel.Position = [38 26 37 21];
        Masking.Surface_ResetLabel.Text = 'Reset';
        
        % Create Surface_UseCellMaskLabel
        Masking.Surface_UseCellMaskLabel = uilabel(Masking.Surface_CleanPanel);
        Masking.Surface_UseCellMaskLabel.Position = [99 24 83 24];
        Masking.Surface_UseCellMaskLabel.Text = 'Use Cell Mask';
        
        % Create UseCellmaskCheckBox
        Masking.UseCellmaskCheckBox = uicheckbox(Masking.Surface_CleanPanel);
        Masking.UseCellmaskCheckBox.Text = '';
        Masking.UseCellmaskCheckBox.Position = [127 8 25 22];
        
        % Show the figure after all components are created
        Masking.UIFigure.Visible = 'on';
        %------------------- set the callbacks for these items
        Masking.Cell_MaskButton.ButtonPushedFcn = @(btn,event) callback_CellMaskButton();
        Masking.Cell_CleanUpMaskButton.ButtonPushedFcn = @(btn,event) callback_CellMaskCleanUpButton();
        %           Masking.Cell_ResetCheckBox;
        Masking.Surface_MaskButton.ButtonPushedFcn = @(btn,event) callback_SurfaceMaskButton();
        Masking.Surface_CleanUpMaskButton.ButtonPushedFcn = @(btn,event) callback_SurfaceMaskCleanUpButton();
        %           Masking.Surface_ResetCheckBox.ButtonPushedFcn;
        Masking.Cell_MaskDropDown.Items ={'Medium'};
        Masking.Cell_MaskDropDown.Value = 'Medium';
        Masking.Cell_CleanUpMethodDropDown.Items = {'Manual'};
        Masking.Cell_CleanUpMethodDropDown.Value = 'Manual';
        Masking.Surface_MaskDropDown.Items = {'High','Medium','Low'};
        Masking.Surface_MaskDropDown.Value = 'Medium';
        Masking.Surface_CleanUpMethodDropDown.Items = {'Manual','By Frame'};
        Masking.Surface_CleanUpMethodDropDown.Value = 'Manual';
        Masking.UseCellmaskCheckBox.Value = 1;
    end
% ----------- call backs for Masking
    function callback_CellMaskButton()
        switch Masking.Cell_MaskDropDown.Value
            case 'Medium'
                obj.cellFill.mask_img();
                obj.distmask = [];
        end
        app.Tree.SelectedNodes = [app.CellMarker_MaskNode,app.CellMarker_ImageNode];
        callback_TreeNodeChange();
        %callback_ViewSelectedButton();
        maskchangebool = 1;
         obj.distmask = [];
    end
    function callback_CellMaskCleanUpButton()
        app.Tree.SelectedNodes = [app.CellMarker_MaskNode,app.CellMarker_ImageNode];
        %callback_ViewSelectedButton();
        maskchangebool = 1;
    end
    function callback_SurfaceMaskButton()
        switch Masking.Surface_MaskDropDown.Value
            case 'High'
                obj.surfaceCargo.mask_img_highsens();
            case 'Medium'
                obj.surfaceCargo.mask_img();
            case 'Low'
                obj.surfaceCargo.mask_img_lowsens();
        end
        obj.cleanedcargomask = obj.surfaceCargo.mask;
        app.Tree.SelectedNodes = [app.SurfSign_MaskNode,app.SurfSign_ImageNode];
        callback_TreeNodeChange();
        %callback_ViewSelectedButton();
        initializePlottingComponents();
        maskchangebool = 1;
    end
    function callback_SurfaceMaskCleanUpButton()
        if Masking.Surface_ResetCheckBox.Value
            reset = 1;
        else
            reset = 0;
        end
        if Masking.UseCellmaskCheckBox.Value
            cellmaskbool = 1;
        else
            cellmaskbool = 0;
        end
        switch Masking.Surface_CleanUpMethodDropDown.Value
            case 'Manual'
                obj.cleanSurfaceCargoMask_Manual(reset,cellmaskbool);
            case 'By Frame'
                obj.cleanSurfaceCargoMaskbyFrame_Manual(reset,cellmaskbool);
        end
        app.Tree.SelectedNodes = [app.SurfSign_MaskNode,app.SurfSign_ImageNode];
        callback_TreeNodeChange();
        %callback_ViewSelectedButton();
        maskchangebool = 1;
    end
% ----------- call backs for imloader
    function callback_CellMarker_LoadImage()
        ImageLoader_loadfile('CellMarker');
    end
    function callback_SurfSignal_LoadImage()
        ImageLoader_loadfile('SurfSign');
    end
    function checkIfBothImages()
        if strcmp(imloader.SurfSign.ViewImageButton.Enable,'on') && strcmp(imloader.CellMarker.ViewImageButton.Enable,'on')
            imloader.SetImagesButton.Enable = 'on';
            imloader.DriftCorrectImageCheckBox.Enable = 'on';
            imloader.DrifCorrText.Enable = 'on';
        else
            imloader.SetImagesButton.Enable = 'off';
            imloader.DriftCorrectImageCheckBox.Enable = 'off';
            imloader.DrifCorrText.Enable = 'off';
        end
    end
    function callback_CellMarker_ViewImage()
        ImageLoader_viewImage('CellMarker');
    end
    function callback_SurfSignal_ViewImage()
        ImageLoader_viewImage('SurfSign');
    end
    function imp = ImageLoader_viewImage(channelstr)
        if strcmp(imloader.(channelstr).ChannelSpinner.Enable,'on')
            imp = MatIJ.showImage(imloader.(channelstr).image(:,:,:,imloader.(channelstr).ChannelSpinner.Value));
        else
            imp = MatIJ.showImage(imloader.(channelstr).image);
        end
    end
    function callback_SurfaceSignal_LoadImageDropDown()
        ijname = imloader.SurfSign.LoadImageDropDown.Value;
        if strcmp(ijname,defaultdropdowntxt)
            return
        end
        ImageLoader_loadfile('SurfSign',ijname);
         callback_updateFileDropdown();
    end
    function callback_CellMarker_LoadImageDropDown()
        ijname = imloader.CellMarker.LoadImageDropDown.Value;
        if strcmp(ijname,defaultdropdowntxt)
            return
        end
        ImageLoader_loadfile('CellMarker',ijname);
        callback_updateFileDropdown();
    end
    function ImageLoader_loadfile(channelstr,ijname)
        % this function brings an image into matlab to the gui two possible ways:
        % channelstr: str, either 'CellMarker' or 'SurfaceSign' specifying
        %   whether to set the image for which channel (correlated with the
        %   imloader structure
        % ijname: str with the name of the image window name open in the
        % instance of Image J running
        % if nargin<2 then prompt the user to load a tif file. if an ijname
        % is specified then attempt to open the image from image j
        if nargin<2 || isempty(ijname)
            [image,name] = loadFile();
        else
            try
                image = squeeze(MatIJ.pullimage(ijname));
                name = ijname;
            catch
            end
        end
        if ~image
            return
        end
        imloader.(channelstr).image= image;
        imloader.(channelstr).ImageNameEditField.Value = name;
        if ndims(image)>3 % need to guess which channel is the color channel
            if size(image,3)<size(image,4) % 3rd channel is color channel
                colch = 3;
            else
                colch = 4;
            end
            imloader.(channelstr).ChannelSpinner.Enable = 'on';
            imloader.(channelstr).ChannelSpinnerLabel.Enable = 'on';
            imloader.(channelstr).ChannelSpinner.Limits = [1 size(image,colch)];
        else
            imloader.(channelstr).ChannelSpinner.Enable = 'off';
            imloader.(channelstr).ChannelSpinnerLabel.Enable = 'off';
        end
        imloader.(channelstr).ViewImageButton.Enable = 'on';
        checkIfBothImages(); %enable set button if both images have been loaded
    end
    function [image,name] = loadFile()
        [filename, pathname] = uigetfile({'*.tif;*.tiff;*.TIF;*.TIFF','Tiff Files'},'Pick a file');
        if filename %make sure user actually selected a file
            name = fullfile(pathname,filename);
            image = loadtiff(name); % uiopen should load tif file and save variable as image to base workspace
        else
            image = 0;
            name = 0;
        end
    end
    function callback_SetImageButton()
        obj.clearimages(); %first re-set everything
        maskchangebool = 1;
        % ----first set the raw images to the class
        if strcmp(imloader.CellMarker.ChannelSpinner.Enable,'on')
            obj.cellFill.setimage(imloader.CellMarker.image(:,:,:,imloader.CellMarker.ChannelSpinner.Value));
        else
            obj.cellFill.setimage(imloader.CellMarker.image);
        end
        if strcmp(imloader.SurfSign.ChannelSpinner.Enable,'on')
            obj.surfaceCargo.setimage(imloader.SurfSign.image(:,:,:,imloader.SurfSign.ChannelSpinner.Value));
        else
            obj.surfaceCargo.setimage(imloader.SurfSign.image);
        end
        % ---- then drift correct
        if imloader.DriftCorrectImageCheckBox.Value
            [cfim,sv_arr] = GeneralAnalysis.timedriftCorrect(obj.cellFill.rawimage);
            sfim = GeneralAnalysis.applydriftCorrect(obj.surfaceCargo.rawimage,sv_arr);
            obj.cellFill.setimage(cfim);
            obj.surfaceCargo.setimage(sfim);
        end
        % --- also set the image values for the class
        obj.cellFill.image = obj.cellFill.rawimage;
        obj.surfaceCargo.image = obj.surfaceCargo.rawimage;
        % ---- now turn inactive because already set
        imloader.SetImagesButton.Enable = 'off';
        initializePlottingComponents();
        app.Tree.SelectedNodes = app.CellMarker_ImageNode;
        callback_TreeNodeChange();
    end
    function callback_DriftCorrectBoxChanged()
        imloader.SetImagesButton.Enable = 'on';
    end
    function callback_CropImage()
        obj.cropimage();
        imloader.ResetButton.Value = false;
    end
    function callback_ResetButton()
        if imloader.ResetButton.Value
        obj.resetimages();
        else
            imloader.ResetButton.Value = true;
        end
    end
    function callback_updateFileDropdown()
        updatefilenames();
        if ~isempty(filenames)
            imloader.CellMarker.LoadImageDropDown.Items = [{'Load Open Image'},filenames];
            imloader.SurfSign.LoadImageDropDown.Items = [{'Load Open Image'},filenames];
        end
    end
    function updatefilenames()
        try
            filenames = cell(ij.WindowManager.getImageTitles())';
        catch
            filenames = [];
        end
    end
% ---------- call backs for app.SurfaceSignalAccumulationBase_UIFigure ----
    function callback_ImageLoaderButton()
        if ~isfield(imloader, 'ImageLoader') || ~isvalid(imloader.ImageLoader)
            createImageLoader();
        else
            imloader.ImageLoader.Visible = 'off';
            imloader.ImageLoader.Visible = 'on';
        end
    end
    function callback_ViewSelectedButton()
        if isempty(app.Tree.SelectedNodes)
            return;
        end
        if numel(app.Tree.SelectedNodes)>1  
            image = [];
            imagename = 'Overlay';
        for n = 1: numel(app.Tree.SelectedNodes)
            if n==1
                [loadedimage,~] = checkViewInput(app.Tree.SelectedNodes(n).Text);
                image = loadedimage;
            else
                [loadedimage,~] = checkViewInput(app.Tree.SelectedNodes(n).Text);
                % check that images are the same size - if not, then copy
                % the smaller sized image (ie repmat for the number of
                % frames)
                if size(loadedimage,3)>size(image,3)
                    repvec = ones(1,ndims(loadedimage));
                    repvec(3) = size(loadedimage,3);
                    image = repmat(image,repvec);
                elseif size(image,3)>size(loadedimage,3)
                    repvec = ones(1,ndims(image));
                    repvec(3) = size(image,3);
                    loadedimage = repmat(loadedimage,repvec);
                end
                image = cat(4,image,loadedimage);
            end
        end 
        else
            [image,imagename] = checkViewInput(app.Tree.SelectedNodes.Text);
        end
        if ~isempty(image)
        imp = MatIJ.showImage(image,'YXTCZ');
        imp.setTitle(imagename);
        end
    end
    function callback_MaskingButton()
        % check if the window is already open, if so then just bring it to
        % the front
        if ~isfield(Masking, 'UIFigure') || ~isvalid(Masking.UIFigure)
        createMasking();
        else
            Masking.UIFigure.Visible = 'off';
            Masking.UIFigure.Visible = 'on';
        end
    end
    function callback_TreeNodeChange()
        if isempty(app.Tree.SelectedNodes)
            return;
        end
        selectedtext = app.Tree.SelectedNodes.Text;
        [currimage,~] = checkViewInput(selectedtext);
        switch ndims(currimage)
            case 1
                errordlg('Need to have an image (at least 2-D array)');
                return;
            case 2
                image_out = single(currimage);
            case 3
                image_out = single(max(currimage,[],3));
            case 4 % figure out which dimension is prob the time dimension and max project over that - just choose the first channel
                if size(currimage,3)<=size(currimage,4)
                    timedim = 4;
                else
                    timedim = 3;
                end
                currimage = squeeze(single(max(currimage,[],timedim)));
                image_out = currimage(:,:,1);
            otherwise
                currimage = squeeze(currimage); % reduce excess dims
                currimage = currimage(:,:,:,ones(1,ndims(currimage)-3)); % slice just into the first three dims (take first of all others)
                image_out = single(max(currimage,[],3));
        end
        if ~isempty(image_out)
            imshow(image_out, [0 prctile(image_out(:),99)],'Parent', app.UIAxes_ImageViewer);
        end
    end
    function callback_SelectRegionsButton()
        if isempty(obj.cellFill.mask)
            f = errordlg(['First Mask the Cell Image Channel']);
            return;
        end
        if isempty(obj.cellFill.soma_mask) || strcmp(app.Tree.SelectedNodes.Text,'Soma Region')
            obj.cellFill.selectSoma();
            try
                obj.makeDistanceMask();
                initializePlottingComponents();
                maskchangebool = 1;
            catch
                f = errordlg(['Error calculating distance mask from current Soma Region and Cell Marker Mask (Use Select Regions or Masking Button)']);
            end
        end
        if isempty(obj.cellFill.AIS_mask) || strcmp(app.Tree.SelectedNodes.Text,'AIS Region')
            obj.cellFill.selectAIS();
            initializePlottingComponents();
            maskchangebool = 1;
        end
        callback_TreeNodeChange()
    end
    function callback_PixelSizeDropdownChange()
        switch app.PixelSizeDropDown.Value
            case 'Bin 1x1 = 0.114 um'
                app.PixelSizeumEditField.Value = 0.114;
            case 'Bin 2x2 = 0.228 um'
                app.PixelSizeumEditField.Value = 0.228;
        end
        obj.pxsize = app.PixelSizeumEditField.Value;
        makePlot();
    end
    function callback_PixelFieldChange()
        obj.pxsize = app.PixelSizeumEditField.Value;
        makePlot();
    end
    function initializePixelSize()
        if isempty(obj.pxsize)
            app.PixelSizeumEditField.Value = 0.114;
            obj.pxsize = app.PixelSizeumEditField.Value;
        else
            app.PixelSizeumEditField.Value = obj.pxsize;
        end
    end
    function initializePlottingComponents()
        if isempty(obj.M)
            if isempty(obj.cleanedcargomask)
                app.TemporalHeatMapButton.Enable = 'off';
            else
                app.TemporalHeatMapButton.Enable = 'on';
            end
            app.EndFramefornormalizationEditField.Value = 0;
            app.TemporalHeatMapButton.Enable = 'off';
            app.SavetoExcelButton.Enable = 'off';
            app.ViewResultsArrayButton.Enable = 'off';
            app.RawResultsVal_CheckBox.Enable = 'off';
            app.RawResultsValRawValsLabel.Enable = 'off';
        else
            app.TemporalHeatMapButton.Enable = 'on';
            app.SavetoExcelButton.Enable = 'on';
            app.ViewResultsArrayButton.Enable = 'on';
            app.RawResultsVal_CheckBox.Enable = 'on';
            app.RawResultsValRawValsLabel.Enable = 'on';
            if app.EndFramefornormalizationEditField.Value == 0
                app.EndFramefornormalizationEditField.Value = size(obj.surfaceCargo.mask,3);
            end
        end
        if isempty(obj.M_AIS) || isempty(obj.M_noAIS)
            app.ImageRegionDropDown.Items = {'Total'}; %no AIS calculations so don't have it as an option
        else
            app.ImageRegionDropDown.Items = {'Total','no AIS','AIS'};
        end
    end
    function initializeImagingParams()
        if ~isempty(obj.imagingparams.baseline.framerate)
            app.PreReleaseEditField.Value = obj.imagingparams.baseline.framerate;
        end
        if ~isempty(obj.imagingparams.postrelease.framerate)
            app.PostReleaseEditField.Value = obj.imagingparams.postrelease.framerate;
        end
        if ~isempty(obj.imagingparams.releaseframe)
            app.ReleaseFrameEditField.Value = obj.imagingparams.releaseframe;
        end
    end
    function setImagingParams()
        obj.imagingparams.baseline.framerate = app.PreReleaseEditField.Value;
        obj.imagingparams.postrelease.framerate = app.PostReleaseEditField.Value;
        obj.imagingparams.releaseframe = app.ReleaseFrameEditField.Value;
        makePlot();
    end
    function timevals = calculatePlotTimeAxis()
        % returns time array (in minutes) corresponding to the values in
        % the object   
       lasttime = size(obj.cleanedcargomask,3);
       releaseframe  = app.ReleaseFrameEditField.Value;
       postreleaseframerate = app.PostReleaseEditField.Value;
       prereleaseframerate = app.PreReleaseEditField.Value;
       if isempty(releaseframe) || isempty(postreleaseframerate) || isempty(prereleaseframerate) || ~releaseframe || ~postreleaseframerate || ~prereleaseframerate
           timevals = 1:lasttime;
       else
           postrelease = (1:(lasttime - releaseframe + 1)).*postreleaseframerate;
           prerelease = (-(releaseframe-2):0).*prereleaseframerate;
           timevals  = [prerelease,postrelease];
       end
    end
    function genPlotValues()
        % this function calculates the surface intensity density for all
        % appropriate regions and saves them as part of the object - this
        % is slow...
        dists_um = getDistances_fromstring();
        dists = ceil(dists_um./obj.pxsize);
        if ~isvector(dists)
            errordlg('Enter distances as a list of comma separated values, ie: 18,70,200');
            return;
        end
        obj.calcDensityperTime(dists);
    end
    function [intensvals, areas] = getImageRegionDropDown()
       switch app.ImageRegionDropDown.Value
            case 'Total'
                if isempty(obj.M)
                    intensvals = [];
                else
                intensvals = obj.M.areanormintensity';  
                areas = obj.M.mask_area;
                end
           case 'AIS'
               if isempty(obj.M_AIS)
                   intensvals = [];
               else
                   intensvals = obj.M_AIS.areanormintensity';
                   areas = obj.M_AIS.mask_area;
               end
           case 'no AIS'
               if isempty(obj.M_noAIS)
                   intensvals = [];
               else
                   intensvals = obj.M_noAIS.areanormintensity';
                   areas = obj.M_noAIS.mask_area;
               end
        end 
    end
    function plotvals = getNormalizationANDImageRegion_DropDown()
        intensvals = getImageRegionDropDown();
        endframe = app.EndFramefornormalizationEditField.Value;
        if isempty(intensvals)
            plotvals = [];
        else
            switch app.NormalizationDropDown.Value
                case 'Per Distance Region'
                    plotvals = intensvals./intensvals(endframe,:);
                case 'Compare to First Region'
                    plotvals = intensvals./intensvals(endframe,1);
            end
        end
    end
    function makePlot(tableval,excelval)
        % this function takes data from obj.M/obj.M_AIS/obj.M_noAIS depending on 
        % dropdown selections and plots it in the UIFigure Axes
        % checks if there isn't a mask then just clears the axes and sets end frame to 0
        endframe = app.EndFramefornormalizationEditField.Value;
        if isempty(endframe) || endframe == 0
            if isempty(obj.surfaceCargo.mask)
                endframe = 0;
            else
                endframe = size(obj.surfaceCargo.mask,3);
            end
            app.EndFramefornormalizationEditField.Value = endframe;
        end
       plotvals  = getNormalizationANDImageRegion_DropDown();
        timevals = calculatePlotTimeAxis();
        if isempty(plotvals)
            cla(app.UIAxes_ResultsPlot,'reset');
            return
        end
        plot(app.UIAxes_ResultsPlot,timevals(1:endframe)',plotvals(1:endframe,:));
        app.UIAxes_ResultsPlot.YLim = [0 1];
        distum = getDistances_fromstring();
        dist_strs = arrayfun(@(x) [num2str(x) ' um from Soma'], distum,'UniformOutput', false);
        legend(app.UIAxes_ResultsPlot,dist_strs,'Location','northwest');
        % now make some things or save some things
        if nargin>0
            if tableval
                if app.RawResultsVal_CheckBox.Value %if raw values are checked
                    [intensvals, areas] = getImageRegionDropDown();
                    f_table = figure('Position', [20 20 600 500]);
                    
                    areas_table = uitable(f_table,'Position', [20 20 560 460]);
                    areas_table.Data = areas;
                    areas_table.ColumnName = dist_strs;
                    areas_table.RowName = 'MaskArea_(pixels)';
                    areas_table.ColumnEditable = false;
                    areas_table.Parent.Name = 'Total Intensity in Region (AUs)';
                    
                    results_table = uitable(f_table,'Position', [20 20 500 400]);
                    results_table.Data = [timevals',intensvals];
                    results_table.ColumnName = ['Time (mins)',dist_strs];
                    results_table.ColumnEditable = false;
                else
                f_table = figure('Position', [20 20 600 500]);
                results_table = uitable(f_table,'Position', [20 20 500 400]);
                results_table.Data = [timevals',plotvals(1:endframe,:)];
                results_table.ColumnName = ['Time (mins)',dist_strs];
                results_table.ColumnEditable = false;
                results_table.Parent.Name = 'Surface Signal Density';
                end
            end
            if excelval
            end
        end
    end
    function callback_UpdatePlot()
        disp('Updating Plot.....');
        genPlotValues();
        makePlot();
        maskchangebool = 0;
        initializePlottingComponents();
        disp('Finished Updating Plot.....');
    end
    function callback_CheckPlotButtonState()
        if maskchangebool %if the mask has been changed then need to update plot and remove heatmap node until it is remade with button
            app.UpdatePlotButton.Enable = 'on';
            cla(app.UIAxes_ResultsPlot,'reset');
            %remove HeatmapNode
            obj.cargo_heatmap = [];
            makeTempHeatMapNode()
            app.SAVEButton.BackgroundColor = [1 .7 .7]; %-- also file will need to be resaved
        else
            app.UpdatePlotButton.Enable = 'off';
        end
        if savenamechangebool
            app.SAVEButton.BackgroundColor = [1 .7 .7];
        end       
    end
    function callback_ViewDistancesButton()
%         try
        dist_um = getDistances_fromstring(); %these are in microns
        dists = dist_um/obj.pxsize;
        distanceim = makedistancefigure(dists);
        imp = MatIJ.showImage(distanceim,'YXTCZ');
        imp.setTitle(['Distances: ' app.DistancespixelsEditField.Value]);
%         catch
%             errordlg('Enter distances as a list of comma separated values, ie: 18,70,200');
%         end
    end
    function [h,im] = makeTempHeatMap()
        if isempty(obj.imagingparams.postrelease.framerate)
            f = errordlg(['Need to set the frame rate in the imaging parameters first']);
            return;
        end
        % get max time info here
        default = 160;
        imgparams.maxtime = default;
        if isempty(obj.cargo_heatmap)
            [h,im] = obj.plotCargoHeatMap(1,imgparams);
        else
            [h,im] = obj.plotCargoHeatMap(0,imgparams);
        end
        
    end
    function callback_TempHeatMapButton()
        try
        [~,~] = makeTempHeatMap();
         makeTempHeatMapNode()
        catch
        end
    end
    function imlayerimage = makedistancefigure(dists)
        distim = obj.distmask;
        switch app.ImageRegionDropDown.Value
            case 'Total'
            case 'no AIS'
                if isempty(obj.cellFill.AIS_mask)
                    obj.cellFill.selectAIS();
                end
                distim(obj.cellFill.AIS_mask(:,:,1)) = nan;
            case 'AIS'
                if isempty(obj.cellFill.AIS_mask)
                    obj.cellFill.selectAIS();
                end
                distim(~obj.cellFill.AIS_mask(:,:,1))= nan;
        end
        imlayerimage = zeros([size(distim,2),size(distim,1),size(dists,2)]);
        for dd = 1:size(dists,2)
            if dd==1
                imlayer = distim<=dists(dd);
            else
                imlayer = distim>dists(dd-1) & distim<=dists(dd);
            end
            imlayer = uint16(imlayer);
            imlayerimage(:,:,dd) = imlayer;
        end
    end
    function dists = getDistances_fromstring()
        dists_strs = strsplit(app.DistancespixelsEditField.Value,',');
        dists = cellfun(@(x) str2num(x), dists_strs);
    end
    function callback_DistancesChanged()
        callback_UpdatePlot();
    end
    function makeTempHeatMapNode()
        if ~isempty(obj.cargo_heatmap) % only create temporal heatmap node if there is a obj.cargo_heatmap
            if isfield(app,'TemporalHeatmapNode') && isvalid(app.TemporalHeatmapNode) %if there was a node but it was deleted
                return;
            end
            % Create TemporalHeatmapNode
            app.TemporalHeatmapNode = uitreenode(app.Tree);
            app.TemporalHeatmapNode.Text = 'Temporal Heatmap';
        else
            if isfield(app,'TemporalHeatmapNode') && isvalid(app.TemporalHeatmapNode) %if there was a node delete it
                app.TemporalHeatmapNode.delete;
            end
        end
    end
    function callback_saveAshleyFile()
        if ~savenamechangebool %no namechange just save
            obj.save();
        else % there's been a change
            callback_selectSaveName();
            obj.save();
            savenamechangebool = 0;
        end
    end
    function callback_selectSaveName()
        if isempty(app.SaveNameEditField.Value) %no info from gui
            if isempty(obj.savename) %check if obj has info
                try
                    [file,path] = uiputfile('AshleyFile.m'); %no info from object just prompt
                catch
                end
            else % there is a savename so start with this
                try
                    [file,path] = uiputfile(obj.savename);
                catch
                end
            end
        else %new text in the gui textbox - use this for saving
            try
                [file,path] = uiputfile(app.SaveNameEditField.Value);
            catch
            end
        end
        obj.savename = fullfile(path,file);
        app.SaveNameEditField.Value = obj.savename;
        savenamechangebool = 1;
    end
    function setUpdatedSaveName()
        savenamechangebool = 1;
    end
    function [cellfill_im,cellfill_mask,surfcarg_im,surfcarg_mask] = getImagesANDMasks()
        if isstruct(obj.cellFill.image)
            cellfill_im = obj.cellFill.image.data;
        else
            cellfill_im = uint16(obj.cellFill.image);
        end
        if isstruct(obj.cellFill.mask)
            cellfill_mask =   uint16(obj.cellFill.mask.data);
        else
            cellfill_mask = uint16(obj.cellFill.mask);
        end
        if isstruct(obj.surfaceCargo.image)
            surfcarg_im = obj.surfaceCargo.image.data;
        else
            surfcarg_im = uint16(obj.surfaceCargo.image);
        end
        if isempty(obj.cleanedcargomask)
            if isstruct(obj.surfaceCargo.mask)
                surfcarg_mask =  uint16(obj.surfaceCargo.mask.data);
            else
                surfcarg_mask = uint16(obj.surfaceCargo.mask);
            end
        else
            surfcarg_mask = uint16(obj.cleanedcargomask);
        end
        cellfill_maskval = prctile(cellfill_im(:),60);
        surfcargo_maskval = prctile(surfcarg_im(:),60);
        cellfill_mask = cellfill_mask.*cellfill_maskval;
        surfcarg_mask = surfcarg_mask.*surfcargo_maskval;
    end
    function distmask = getDistanceMask()
        if isstruct(obj.distmask)
            distmask = obj.distmask.data;
        else
            distmask=obj.distmask;
        end
            distmask(isinf(distmask)) = 0;
            distmask = distmask*obj.pxsize;%convert values to microns
            distmask = uint16(distmask);
    end
    function [image,imagename] = checkViewInput(TreeNodeText)
        image = [];
        imagename = '';
        [cellfill_im,cellfill_mask,surfcarg_im,surfcarg_mask] = getImagesANDMasks();
        switch TreeNodeText
            case 'Cell Marker'
                if isempty(obj.cellFill.rawimage)
                    f = errordlg(['Need to Load in a Cell Fill image first (Use Image Loader Button)']);
                    return;
                end
                if ~isempty(cellfill_mask)
                    image = cat(4,cellfill_im,uint16(cellfill_mask));
                else
                    image = cellfill_im;
                end
                imagename = 'CellFill';
                
            case 'Cell Image'
                if isempty(obj.cellFill.rawimage)
                    f = errordlg(['Need to Load in a Cell Fill image first (Use Image Loader Button)']);
                    return;
                end
                image = cellfill_im;
                imagename = 'CellFillImage';
                
            case 'Cell Mask'
                if isempty(cellfill_mask)
                    f = errordlg(['Need to create a Cell Image Mask First (Use Masking Button)']);
                    return;
                end
                image = cellfill_mask;
                imagename = 'CellFillMask';
                
            case 'Surface Signal'
                if isempty(obj.surfaceCargo.rawimage)
                    f = errordlg(['Need to Load in a Surface Cargo image first (Use Image Loader Button)']);
                    return;
                end
                if ~isempty(surfcarg_mask)
                    image = cat(4,surfcarg_im,uint16(surfcarg_mask));
                else
                    image = surfcarg_im;
                end
                imagename = 'SurfaceSignal';
                
            case 'Surface Signal Image'
                if isempty(obj.surfaceCargo.rawimage)
                    f = errordlg(['Need to Load in a Surface Signal Image first (Use Image Loader Button)']);
                    return;
                end
                image = surfcarg_im;
                imagename = 'SurfaceImage';
                
            case 'Surface Signal Mask'
                if isempty(surfcarg_mask)
                    f = errordlg(['Need to create a Surface Signal Mask First (Use Masking Button)']);
                    return;
                end
                image = surfcarg_mask;
                imagename = 'SurfaceMask';
                
            case 'AIS Region'
                if isempty(obj.cellFill.mask)
                    f = errordlg(['First Mask the Cell Image Channel']);
                    return;
                end
                if isempty(obj.cellFill.AIS_mask)
                    obj.cellFill.selectAIS();
                    initializePlottingComponents();
                    maskchangebool = 1;
                end
                image = uint16(obj.cellFill.AIS_mask);
                imagename = 'AIS_mask';
            case 'Soma Region'
                if isempty(obj.cellFill.mask)
                    f = errordlg(['First Mask the Cell Image Channel']);
                    return;
                end
                if isempty(obj.cellFill.soma_mask)
                    obj.cellFill.selectSoma();
                    try
                        obj.makeDistanceMask();
                        initializePlottingComponents();
                        maskchangebool = 1;
                    catch
                        f = errordlg(['Error calculating distance mask from current Soma Region and Cell Marker Mask (Use Select Regions or Masking Button)']);
                    end
                end
                image = uint16(obj.cellFill.soma_mask);
                imagename = 'soma_mask';
            case 'Temporal Heatmap'
                [h,im] = makeTempHeatMap();
                close(h);
                image = uint16(im);
                imagename = 'SurfaceSignal_TemporalHeatMap';
                
            case 'Distance Mask'
                if isempty(obj.cellFill.soma_mask)
                    f = errordlg(['Need to create Select the Soma Region First (Use Select Regions Button)']);
                    return;
                end
                if isempty(obj.distmask)
                    try
                        obj.makeDistanceMask()
                    catch
                        f = errordlg(['Need to create Select the Soma Region or Remask Cell Marker (Use Select Regions or Masking Button)']);
                        return;
                    end
                end
                image = getDistanceMask();
                imagename = 'DistanceMask';
        end
        
    end
end