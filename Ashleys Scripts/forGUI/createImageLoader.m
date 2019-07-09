function createImageLoader()
        % Create ImageLoader
      % Create ImageLoader
            imloader.ImageLoader = uifigure;
            imloader.ImageLoader.Position = [100 100 390 275];
            imloader.ImageLoader.Name = 'Image Loader';

            % Create TabGroup
            imloader.TabGroup = uitabgroup(imloader.ImageLoader);
            imloader.TabGroup.Position = [18 141 357 124];

            % Create CellMarkerTab
            imloader.CellMarkerTab = uitab(imloader.TabGroup);
            imloader.CellMarkerTab.Title = 'Cell Marker';

            % Create CellMarker_ImageNameLabel
            imloader.CellMarker_ImageNameLabel = uilabel(imloader.CellMarkerTab);
            imloader.CellMarker_ImageNameLabel.Position = [13 40 365 22];
            imloader.CellMarker_ImageNameLabel.Text = 'Image: ';

            % Create CellMarker_LoadImageDropDown
            imloader.CellMarker_LoadImageDropDown = uidropdown(imloader.CellMarkerTab);
            imloader.CellMarker_LoadImageDropDown.Items = {'Load Open Image', 'Option 1', 'Option 2', 'Option 3'};
            imloader.CellMarker_LoadImageDropDown.Position = [197 64 151 21];
            imloader.CellMarker_LoadImageDropDown.Value = 'Load Open Image';

            % Create CellMarker_ORtext
            imloader.CellMarker_ORtext = uilabel(imloader.CellMarkerTab);
            imloader.CellMarker_ORtext.Position = [175 63 23 22];
            imloader.CellMarker_ORtext.Text = 'or';

            % Create CellMarker_LoadImageFromFileButton
            imloader.CellMarker_LoadImageFromFileButton = uibutton(imloader.CellMarkerTab, 'push');
            imloader.CellMarker_LoadImageFromFileButton.Position = [13 63 149 22];
            imloader.CellMarker_LoadImageFromFileButton.Text = 'Load Image From File';

            % Create CellMarker_ViewImageButton
            imloader.CellMarker_ViewImageButton = uibutton(imloader.CellMarkerTab, 'push');
            imloader.CellMarker_ViewImageButton.Position = [121 11 99 22];
            imloader.CellMarker_ViewImageButton.Text = 'View Image';

            % Create ChannelSpinnerLabel
            imloader.ChannelSpinnerLabel = uilabel(imloader.CellMarkerTab);
            imloader.ChannelSpinnerLabel.HorizontalAlignment = 'right';
            imloader.ChannelSpinnerLabel.Enable = 'off';
            imloader.ChannelSpinnerLabel.Position = [13 11 45 22];
            imloader.ChannelSpinnerLabel.Text = 'Channel';

            % Create CellMarker_ChannelSpinner
            imloader.CellMarker_ChannelSpinner = uispinner(imloader.CellMarkerTab);
            imloader.CellMarker_ChannelSpinner.Enable = 'off';
            imloader.CellMarker_ChannelSpinner.Position = [66 11 39 22];

            % Create SurfaceSignalTab
            imloader.SurfaceSignalTab = uitab(imloader.TabGroup);
            imloader.SurfaceSignalTab.Title = 'Surface Signal';

            % Create SurfSign_ImageNameLabel
            imloader.SurfSign_ImageNameLabel = uilabel(imloader.SurfaceSignalTab);
            imloader.SurfSign_ImageNameLabel.Position = [13 40 365 22];
            imloader.SurfSign_ImageNameLabel.Text = 'Image: ';

            % Create SurfSign_LoadImageDropDown
            imloader.SurfSign_LoadImageDropDown = uidropdown(imloader.SurfaceSignalTab);
            imloader.SurfSign_LoadImageDropDown.Items = {'Load Open Image', 'Option 1', 'Option 2', 'Option 3'};
            imloader.SurfSign_LoadImageDropDown.Position = [197 64 151 21];
            imloader.SurfSign_LoadImageDropDown.Value = 'Load Open Image';

            % Create SurfSign_ORtext
            imloader.SurfSign_ORtext = uilabel(imloader.SurfaceSignalTab);
            imloader.SurfSign_ORtext.Position = [175 63 23 22];
            imloader.SurfSign_ORtext.Text = 'or';

            % Create SurfSign_LoadImageFromFileButton
            imloader.SurfSign_LoadImageFromFileButton = uibutton(imloader.SurfaceSignalTab, 'push');
            imloader.SurfSign_LoadImageFromFileButton.Position = [13 63 149 22];
            imloader.SurfSign_LoadImageFromFileButton.Text = 'Load Image From File';

            % Create SurfSign_CellMarker_ViewImageButton
            imloader.SurfSign_CellMarker_ViewImageButton = uibutton(imloader.SurfaceSignalTab, 'push');
            imloader.SurfSign_CellMarker_ViewImageButton.Position = [121 11 99 22];
            imloader.SurfSign_CellMarker_ViewImageButton.Text = 'View Image';

            % Create ChannelSpinner_2Label
            imloader.ChannelSpinner_2Label = uilabel(imloader.SurfaceSignalTab);
            imloader.ChannelSpinner_2Label.HorizontalAlignment = 'right';
            imloader.ChannelSpinner_2Label.Enable = 'off';
            imloader.ChannelSpinner_2Label.Position = [13 11 45 22];
            imloader.ChannelSpinner_2Label.Text = 'Channel';

            % Create SurfSign_ChannelSpinner
            imloader.SurfSign_ChannelSpinner = uispinner(imloader.SurfaceSignalTab);
            imloader.SurfSign_ChannelSpinner.Enable = 'off';
            imloader.SurfSign_ChannelSpinner.Position = [66 11 39 22];

            % Create SetPanel
            imloader.SetPanel = uipanel(imloader.ImageLoader);
            imloader.SetPanel.Position = [18 71 357 71];

            % Create DriftCorrectImageCheckBox
            imloader.DriftCorrectImageCheckBox = uicheckbox(imloader.SetPanel);
            imloader.DriftCorrectImageCheckBox.Text = 'Drift Correct Image';
            imloader.DriftCorrectImageCheckBox.Position = [150 30 124 23];

            % Create SetImagesButton
            imloader.SetImagesButton = uibutton(imloader.SetPanel, 'push');
            imloader.SetImagesButton.Position = [12 14 132 40];
            imloader.SetImagesButton.Text = 'Set Images';

            % Create DrifCorrText
            imloader.DrifCorrText = uilabel(imloader.SetPanel);
            imloader.DrifCorrText.Position = [150 14 135 22];
            imloader.DrifCorrText.Text = '(using Cell Marker drift)';

            % Create CropPanel
            imloader.CropPanel = uipanel(imloader.ImageLoader);
            imloader.CropPanel.Position = [18 16 357 56];

            % Create CropActiveImageRegionButton
            imloader.CropActiveImageRegionButton = uibutton(imloader.CropPanel, 'push');
            imloader.CropActiveImageRegionButton.Position = [12 11 149 35];
            imloader.CropActiveImageRegionButton.Text = 'Crop Active Image Region';
 
end