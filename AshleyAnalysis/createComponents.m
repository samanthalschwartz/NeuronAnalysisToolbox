function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 391 389];
            app.UIFigure.Name = 'Surface Signal Accumulation';

            % Create Tree
            app.Tree = uitree(app.UIFigure);
            app.Tree.Position = [28 22 183 323];
            app.Tree.Multiselect = 'on';

            % Create CellMarkerNode
            app.CellMarkerNode = uitreenode(app.Tree);
            app.CellMarkerNode.Text = 'Cell Marker';

            % Create ImageNode
            app.ImageNode = uitreenode(app.CellMarkerNode);
            app.ImageNode.Text = 'Image';

            % Create MaskNode
            app.MaskNode = uitreenode(app.CellMarkerNode);
            app.MaskNode.Text = 'Mask';

            % Create DistanceMaskNode
            app.DistanceMaskNode = uitreenode(app.CellMarkerNode);
            app.DistanceMaskNode.Text = 'Distance Mask';

            % Create SurfaceSignalNode
            app.SurfaceSignalNode = uitreenode(app.Tree);
            app.SurfaceSignalNode.Text = 'Surface Signal';

            % Create ImageNode_2
            app.ImageNode_2 = uitreenode(app.SurfaceSignalNode);
            app.ImageNode_2.Text = 'Image';

            % Create MaskNode_2
            app.MaskNode_2 = uitreenode(app.SurfaceSignalNode);
            app.MaskNode_2.Text = 'Mask';

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

            % Create SaveNameNode
            app.SaveNameNode = uitreenode(app.FileInfoNode);
            app.SaveNameNode.Text = 'Save Name';

            % Create OriginalTiffFilesNode
            app.OriginalTiffFilesNode = uitreenode(app.FileInfoNode);
            app.OriginalTiffFilesNode.Text = 'Original Tiff Files';

            % Create SurfaceSignalAccumulationLabel
            app.SurfaceSignalAccumulationLabel = uilabel(app.UIFigure);
            app.SurfaceSignalAccumulationLabel.Position = [89 357 159 22];
            app.SurfaceSignalAccumulationLabel.Text = 'Surface Signal Accumulation';

            % Create ViewSelectedButton
            app.ViewSelectedButton = uibutton(app.UIFigure, 'push');
            app.ViewSelectedButton.Position = [231 71 100 22];
            app.ViewSelectedButton.Text = 'View Selected';

            % Create OverlayButton
            app.OverlayButton = uibutton(app.UIFigure, 'push');
            app.OverlayButton.Position = [231 50 100 22];
            app.OverlayButton.Text = 'Overlay';

            % Create MethodsPanel
            app.MethodsPanel = uipanel(app.UIFigure);
            app.MethodsPanel.Title = 'Methods';
            app.MethodsPanel.Position = [231 154 133 191];

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

            % Create DefineSaveInfoButton
            app.DefineSaveInfoButton = uibutton(app.MethodsPanel, 'push');
            app.DefineSaveInfoButton.Position = [11 51 104 22];
            app.DefineSaveInfoButton.Text = 'Define Save Info';

            % Create QuantifyButton
            app.QuantifyButton = uibutton(app.MethodsPanel, 'push');
            app.QuantifyButton.Position = [13 18 100 22];
            app.QuantifyButton.Text = 'Quantify';
        end