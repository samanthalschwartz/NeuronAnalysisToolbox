function guiFig = gui(obj)
    app = struct();
    createComponents();
    obj.guiFig = app.UIFigure;
    function createComponents(app)
        
        % Create SurfaceSignalAccumulationBase_UIFigure
        app.SurfaceSignalAccumulationBase_UIFigure = uifigure;
        app.SurfaceSignalAccumulationBase_UIFigure.Position = [100 100 391 425];
        app.SurfaceSignalAccumulationBase_UIFigure.Name = 'Surface Signal Accumulation';
        
        % Create Tree
        app.Tree = uitree(app.SurfaceSignalAccumulationBase_UIFigure);
        app.Tree.Position = [28 26 183 355];
        
        % Create CellMarkerNode
        app.CellMarkerNode = uitreenode(app.Tree);
        app.CellMarkerNode.Text = 'Cell Marker';
        
        % Create CellMarker_ImageNode
        app.CellMarker_ImageNode = uitreenode(app.CellMarkerNode);
        app.CellMarker_ImageNode.Text = 'Image';
        
        % Create CellMarker_MaskNode
        app.CellMarker_MaskNode = uitreenode(app.CellMarkerNode);
        app.CellMarker_MaskNode.Text = 'Mask';
        
        % Create CellMarker_DistanceMaskNode
        app.CellMarker_DistanceMaskNode = uitreenode(app.CellMarkerNode);
        app.CellMarker_DistanceMaskNode.Text = 'Distance Mask';
        
        % Create SurfaceSignalNode
        app.SurfaceSignalNode = uitreenode(app.Tree);
        app.SurfaceSignalNode.Text = 'Surface Signal';
        
        % Create SurfSign_ImageNode
        app.SurfSign_ImageNode = uitreenode(app.SurfaceSignalNode);
        app.SurfSign_ImageNode.Text = 'Image';
        
        % Create SurfSign_MaskNode
        app.SurfSign_MaskNode = uitreenode(app.SurfaceSignalNode);
        app.SurfSign_MaskNode.Text = 'Mask';
        
        % Create AISRegionNode
        app.AISRegionNode = uitreenode(app.Tree);
        app.AISRegionNode.Text = 'AIS Region';
        
        % Create SomaRegionNode
        app.SomaRegionNode = uitreenode(app.Tree);
        app.SomaRegionNode.Text = 'Soma Region';
        
        % Create QuantificationNode
        app.QuantificationNode = uitreenode(app.Tree);
        app.QuantificationNode.Text = 'Quantification';
        
        % Create TemporalHeatmapNode
        app.TemporalHeatmapNode = uitreenode(app.QuantificationNode);
        app.TemporalHeatmapNode.Text = 'Temporal Heatmap';
        
        % Create DensityQuantificationNode
        app.DensityQuantificationNode = uitreenode(app.QuantificationNode);
        app.DensityQuantificationNode.Text = 'Density Quantification';
        
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
        
        % Create SurfaceSignalAccumulationLabel
        app.SurfaceSignalAccumulationLabel = uilabel(app.SurfaceSignalAccumulationBase_UIFigure);
        app.SurfaceSignalAccumulationLabel.Position = [89 393 159 22];
        app.SurfaceSignalAccumulationLabel.Text = 'Surface Signal Accumulation';
        
        % Create ViewSelectedButton
        app.ViewSelectedButton = uibutton(app.SurfaceSignalAccumulationBase_UIFigure, 'push');
        app.ViewSelectedButton.Position = [231 107 100 22];
        app.ViewSelectedButton.Text = 'View Selected';
        
        % Create OverlayButton
        app.OverlayButton = uibutton(app.SurfaceSignalAccumulationBase_UIFigure, 'push');
        app.OverlayButton.Position = [231 86 100 22];
        app.OverlayButton.Text = 'Overlay';
        
        % Create MethodsPanel
        app.MethodsPanel = uipanel(app.SurfaceSignalAccumulationBase_UIFigure);
        app.MethodsPanel.Title = 'Methods';
        app.MethodsPanel.Position = [231 190 133 191];
        
        % Create ImageLoaderButton
        app.ImageLoaderButton = uibutton(app.MethodsPanel, 'push');
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
    end

    function callback_ViewSelectedButton()
        display(app.Tree.SelectedNodes);
    end
    function callback_Overlay()
    end
    function callback_ImageLoaderButton()
    end
    function callback_MaskingButton()
    end
    function callback_SelectRegionsButton()
    end
    function callback_DefineImageInfoButton()
    end
    function callback_QuantifyButton()
    end









end