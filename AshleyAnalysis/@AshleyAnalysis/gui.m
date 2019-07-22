function guiFig = gui(obj)
    if ishandle(obj.guiFig)
        figure(obj.guiFig);
        return
    end
%--- create structures for all guis associated with this gui----
    app = struct();
    imloader = struct();
    % --- global vars related to imloader
    
    %------------------------------------
%--- now actually make the base gui    
    createComponents();
    obj.guiFig = app.SurfaceSignalAccumulationBase_UIFigure;
%----
% -------------  gui makers -------------------
    function createComponents()
        
        % Create SurfaceSignalAccumulationBase_UIFigure
        app.SurfaceSignalAccumulationBase_UIFigure = uifigure;
        app.SurfaceSignalAccumulationBase_UIFigure.Position = [100 100 391 425];
        app.SurfaceSignalAccumulationBase_UIFigure.Name = 'Surface Signal Accumulation';
        
        % Create Tree
        app.Tree = uitree(app.SurfaceSignalAccumulationBase_UIFigure,'Multiselect','on');
        app.Tree.Position = [28 26 183 355];
        
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
             
        % Create AISRegionNode
        app.AISRegionNode = uitreenode(app.Tree);
        app.AISRegionNode.Text = 'AIS Region';
        
        % Create SomaRegionNode
        app.SomaRegionNode = uitreenode(app.Tree);
        app.SomaRegionNode.Text = 'Soma Region';
        
        % Create CellMarker_DistanceMaskNode
        app.DistanceMaskNode = uitreenode(app.Tree);
        app.DistanceMaskNode.Text = 'Distance Mask';
        
        % Create FileInfoNode
        app.FileInfoNode = uitreenode(app.Tree);
        app.FileInfoNode.Text = 'File Info';
        
        % Create ImagingParametersNode
        app.ImagingParametersNode = uitreenode(app.FileInfoNode);
        app.ImagingParametersNode.Text = 'Imaging Parameters';
        
        % Create SaveNameNode
        app.SaveNameNode = uitreenode(app.FileInfoNode);
        app.SaveNameNode.Text = 'Save Name';
        
        % Create OriginalTiffFilesNode
        app.OriginalTiffFilesNode = uitreenode(app.FileInfoNode);
        app.OriginalTiffFilesNode.Text = 'Original Tiff Files';
        
        % Create QuantificationNode
        app.QuantificationNode = uitreenode(app.Tree);
        app.QuantificationNode.Text = 'Quantification';
        
        % Create TemporalHeatmapNode
        app.TemporalHeatmapNode = uitreenode(app.QuantificationNode);
        app.TemporalHeatmapNode.Text = 'Temporal Heatmap';
        
        % Create DensityQuantificationNode
        app.DensityQuantificationNode = uitreenode(app.QuantificationNode);
        app.DensityQuantificationNode.Text = 'Density Quantification';  
               
        % Create SurfaceSignalAccumulationLabel
        app.SurfaceSignalAccumulationLabel = uilabel(app.SurfaceSignalAccumulationBase_UIFigure);
        app.SurfaceSignalAccumulationLabel.Position = [89 393 159 22];
        app.SurfaceSignalAccumulationLabel.Text = 'Surface Signal Accumulation';
        
        % Create ViewSelectedButton
        app.ViewSelectedButton = uibutton(app.SurfaceSignalAccumulationBase_UIFigure, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_ViewSelectedButton());
        app.ViewSelectedButton.Position = [231 107 100 22];
        app.ViewSelectedButton.Text = 'View Selected';
                
        % Create OverlayButton
        app.OverlayButton = uibutton(app.SurfaceSignalAccumulationBase_UIFigure, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_Overlay());
        app.OverlayButton.Position = [231 86 100 22];
        app.OverlayButton.Text = 'Overlay';
        
        % Create MethodsPanel
        app.MethodsPanel = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
        app.MethodsPanel.Title = 'Methods';
        app.MethodsPanel.Position = [231 190 133 191];
        
        % Create ImageLoaderButton
        app.ImageLoaderButton = uibutton(app.MethodsPanel, 'push',...
            'ButtonPushedFcn', @(btn,event) callback_ImageLoaderButton());
        app.ImageLoaderButton.Position = [13 142 100 22];
        app.ImageLoaderButton.Text = 'Image Loader';
        
        % Create MaskingButton
        app.MaskingButton = uibutton(app.MethodsPanel, 'push');
        app.MaskingButton.Position = [13 113 100 22];
        app.MaskingButton.Text = 'Masking';
        
        % Create SelectRegionsButton
        app.SelectRegionsButton = uibutton(app.MethodsPanel, 'push');
        app.SelectRegionsButton.Position = [13 80 100 22];
        app.SelectRegionsButton.Text = 'Select Regions';
        
        % Create DefineImageInfoButton
        app.DefineImageInfoButton = uibutton(app.MethodsPanel, 'push');
        app.DefineImageInfoButton.Position = [8 51 110 22];
        app.DefineImageInfoButton.Text = 'Define Image Info';
        
        % Create QuantifyButton
        app.QuantifyButton = uibutton(app.MethodsPanel, 'push');
        app.QuantifyButton.Position = [13 18 100 22];
        app.QuantifyButton.Text = 'Quantify';
        expand(app.Tree);
    end
    function createImageLoader()
        imloader.ImageLoader = uifigure;
        imloader.ImageLoader.Position = [100 100 390 275];
        imloader.ImageLoader.Name = 'Image Loader';
        
        % Create TabGroup
        imloader.TabGroup = uitabgroup(imloader.ImageLoader);
        imloader.TabGroup.Position = [18 141 357 124];
        
        % Create CellMarkerTab
        imloader.CellMarkerTab = uitab(imloader.TabGroup);
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
        imloader.CellMarker.LoadImageDropDown = uidropdown(imloader.CellMarkerTab);
        imloader.CellMarker.LoadImageDropDown.Items = {'Load Open Image', 'Option 1', 'Option 2', 'Option 3'};
        imloader.CellMarker.LoadImageDropDown.Position = [197 64 151 21];
        imloader.CellMarker.LoadImageDropDown.Value = 'Load Open Image';
        
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
        imloader.SurfaceSignalTab = uitab(imloader.TabGroup);
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
        imloader.SurfSign.LoadImageDropDown = uidropdown(imloader.SurfaceSignalTab);
        imloader.SurfSign.LoadImageDropDown.Items = {'Load Open Image', 'Option 1', 'Option 2', 'Option 3'};
        imloader.SurfSign.LoadImageDropDown.Position = [197 64 151 21];
        imloader.SurfSign.LoadImageDropDown.Value = 'Load Open Image';
        
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
        imloader.DriftCorrectImageCheckBox.Text = 'Drift Correct Image';
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
        imloader.CropActiveImageRegionButton = uibutton(imloader.CropPanel, 'push');
        imloader.CropActiveImageRegionButton.Position = [12 11 149 35];
        imloader.CropActiveImageRegionButton.Text = 'Crop Active Image Region';
        
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
    function ImageLoader_loadfile(channelstr)
        [image,name] = getFile();
        if ~image
            return
        end
        imloader.(channelstr).image= image;
        imloader.(channelstr).ImageNameEditField.Value = name;
        if ndims(image)>3 % assume 4th dimension is channel info
            imloader.(channelstr).ChannelSpinner.Enable = 'on';
            imloader.(channelstr).ChannelSpinnerLabel.Enable = 'on';
            imloader.(channelstr).ChannelSpinner.Limits = [1 size(image,4)];
        else
            imloader.(channelstr).ChannelSpinner.Enable = 'off';
            imloader.(channelstr).ChannelSpinnerLabel.Enable = 'off';
        end
        imloader.(channelstr).ViewImageButton.Enable = 'on';
        checkIfBothImages(); %enable set button if both images have been loaded
    end
    function [image,name] = getFile()
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
            sfim = GeneralAnalysis.applydriftCorrect(cfim,sv_arr);
            obj.cellFill.setimage(cfim);
            obj.surfaceCargo.setimage(sfim);
        end
        % --- also set the image values for the class
        obj.cellFill.image = obj.cellFill.rawimage;
        obj.surfaceCargo.image = obj.surfaceCargo.rawimage;
        % ---- now turn inactive because already set
        imloader.SetImagesButton.Enable = 'off';
    end
    function callback_DriftCorrectBoxChanged()
        imloader.SetImagesButton.Enable = 'on';
    end
    function callback_CropImage()
    end
% ---------- call backs for app.SurfaceSignalAccumulationBase_UIFigure ----
    function callback_ImageLoaderButton()
        createImageLoader();
    end
    function callback_ViewSelectedButton()
        a = [];
        if numel(app.Tree.SelectedNodes)>1
            [image,imagename] = checkViewInput(app.Tree.SelectedNodes(1).Text);
        else
            [image,imagename] = checkViewInput(app.Tree.SelectedNodes.Text);
        end
        imp = MatIJ.showImage(image);
        imp.setTitle(imagename);
    end
    function callback_Overlay()
        [cellfill_im,cellfill_mask,surfcarg_im,surfcarg_mask] = getImagesANDMasks();
        OverlayImage = [];
        for n = 1: numel(app.Tree.SelectedNodes)
            if n==1
                [image,~] = checkViewInput(app.Tree.SelectedNodes(n).Text);
                OverlayImage = image;
            else
                [image,~] = checkViewInput(app.Tree.SelectedNodes(n).Text);
                if size(image,3)>size(OverlayImage,3)
                    repvec = ones(1,ndims(image));
                    repvec(3) = size(image,3);
                    OverlayImage = repmat(OverlayImage,repvec);
                elseif size(OverlayImage,3)>size(image,3)
                    repvec = ones(1,ndims(OverlayImage));
                    repvec(3) = size(OverlayImage,3);
                    image = repmat(image,repvec);
                end
                OverlayImage = cat(4,OverlayImage,image);
            end
        end
         imp = MatIJ.showImage(OverlayImage);
    end 
    function callback_MaskingButton()
    end
    function callback_SelectRegionsButton()
    end
    function callback_DefineImageInfoButton()
    end
    function callback_QuantifyButton()
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
        cellfill_maskval = prctile(cellfill_im(:),0.6);
        surfcargo_maskval = prctile(surfcarg_im(:),0.6);
        cellfill_mask = cellfill_mask.*cellfill_maskval;
        surfcarg_mask = surfcarg_mask.*surfcargo_maskval;
    end
    function distmask = getDistanceMask()
        if isstruct(obj.distmask)
            distmask = obj.distmask.data;
            distmask(isinf(distmask)) = 0;
        else
            distmask=obj.distmask;
            distmask(isinf(distmask)) = 0;
            distmask = uint16(distmask);
        end
    end
    function [image,imagename] = checkViewInput(TreeNodeText)
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
            case 'Soma Region'
            case 'Temporal Heatmap'
                
                
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