function maskingComponents(app)
% Create MaskingUIFigure
            app.MaskingUIFigure = uifigure;
            app.MaskingUIFigure.Position = [100 100 519 212];
            app.MaskingUIFigure.Name = 'Masking';

            % Create SurfaceSignalMaskPanel
            app.SurfaceSignalMaskPanel = uipanel(app.MaskingUIFigure);
            app.SurfaceSignalMaskPanel.Title = 'Surface Signal Mask';
            app.SurfaceSignalMaskPanel.Position = [265 19 243 183];

            % Create CleanUpMethodDropDownLabel
            app.CleanUpMethodDropDownLabel = uilabel(app.SurfaceSignalMaskPanel);
            app.CleanUpMethodDropDownLabel.HorizontalAlignment = 'right';
            app.CleanUpMethodDropDownLabel.Position = [16 104 99 22];
            app.CleanUpMethodDropDownLabel.Text = 'Clean Up Method';

            % Create SurfSign_CleanUpMethodDropDown
            app.SurfSign_CleanUpMethodDropDown = uidropdown(app.SurfaceSignalMaskPanel);
            app.SurfSign_CleanUpMethodDropDown.Position = [130 104 100 22];

            % Create MaskingMethodDropDownLabel
            app.MaskingMethodDropDownLabel = uilabel(app.SurfaceSignalMaskPanel);
            app.MaskingMethodDropDownLabel.HorizontalAlignment = 'right';
            app.MaskingMethodDropDownLabel.Position = [21 136 94 22];
            app.MaskingMethodDropDownLabel.Text = 'Masking Method';

            % Create SurfSign_MaskingMethodDropDown
            app.SurfSign_MaskingMethodDropDown = uidropdown(app.SurfaceSignalMaskPanel);
            app.SurfSign_MaskingMethodDropDown.Position = [130 136 100 22];

            % Create SurfSign_ReMasktheImageButton
            app.SurfSign_ReMasktheImageButton = uibutton(app.SurfaceSignalMaskPanel, 'push');
            app.SurfSign_ReMasktheImageButton.Position = [15 45 130 41];
            app.SurfSign_ReMasktheImageButton.Text = 'Re-Mask the Image';

            % Create SurfSign_MaskResetButton
            app.SurfSign_MaskResetButton = uibutton(app.SurfaceSignalMaskPanel, 'state');
            app.SurfSign_MaskResetButton.Text = 'Reset';
            app.SurfSign_MaskResetButton.Position = [144 54 40 22];

            % Create SurfSign_ViewMaskButton
            app.SurfSign_ViewMaskButton = uibutton(app.SurfaceSignalMaskPanel, 'push');
            app.SurfSign_ViewMaskButton.Position = [73 1 100 22];
            app.SurfSign_ViewMaskButton.Text = 'View Mask';

            % Create CellMarkerMaskPanel
            app.CellMarkerMaskPanel = uipanel(app.MaskingUIFigure);
            app.CellMarkerMaskPanel.Title = 'Cell Marker Mask';
            app.CellMarkerMaskPanel.Position = [14 19 243 183];

            % Create CleanUpMethodDropDown_2Label
            app.CleanUpMethodDropDown_2Label = uilabel(app.CellMarkerMaskPanel);
            app.CleanUpMethodDropDown_2Label.HorizontalAlignment = 'right';
            app.CleanUpMethodDropDown_2Label.Position = [16 104 99 22];
            app.CleanUpMethodDropDown_2Label.Text = 'Clean Up Method';

            % Create CellMarker_CleanUpMethodDropDown
            app.CellMarker_CleanUpMethodDropDown = uidropdown(app.CellMarkerMaskPanel);
            app.CellMarker_CleanUpMethodDropDown.Position = [130 104 100 22];

            % Create MaskingMethodDropDown_2Label
            app.MaskingMethodDropDown_2Label = uilabel(app.CellMarkerMaskPanel);
            app.MaskingMethodDropDown_2Label.HorizontalAlignment = 'right';
            app.MaskingMethodDropDown_2Label.Position = [21 136 94 22];
            app.MaskingMethodDropDown_2Label.Text = 'Masking Method';

            % Create CellMarker_MaskingMethodDropDown
            app.CellMarker_MaskingMethodDropDown = uidropdown(app.CellMarkerMaskPanel);
            app.CellMarker_MaskingMethodDropDown.Position = [130 136 100 22];

            % Create CellMarker_ReMasktheImageButton
            app.CellMarker_ReMasktheImageButton = uibutton(app.CellMarkerMaskPanel, 'push');
            app.CellMarker_ReMasktheImageButton.Position = [15 45 130 41];
            app.CellMarker_ReMasktheImageButton.Text = 'Re-Mask the Image';

            % Create CellMarker_MaskResetButton
            app.CellMarker_MaskResetButton = uibutton(app.CellMarkerMaskPanel, 'state');
            app.CellMarker_MaskResetButton.Text = 'Reset';
            app.CellMarker_MaskResetButton.Position = [144 54 41 22];

            % Create CellMarker_ViewMaskButton
            app.CellMarker_ViewMaskButton = uibutton(app.CellMarkerMaskPanel, 'push');
            app.CellMarker_ViewMaskButton.Position = [72 1 100 22];
            app.CellMarker_ViewMaskButton.Text = 'View Mask';
        end