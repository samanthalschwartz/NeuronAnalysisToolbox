classdef AshleyAnalysisGUI < matlab.apps.AppBase
    
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        SelectCellFillImagePanel       matlab.ui.container.Panel
        orLoadFileButton               matlab.ui.control.Button
        ChSpinnerLabel                 matlab.ui.control.Label
        ChSpinner                      matlab.ui.control.Spinner
        MatlabVariableDropDownLabel    matlab.ui.control.Label
        MatlabVariableDropDown         matlab.ui.control.DropDown
        ViewCellFillImageButton        matlab.ui.control.Button
        OpenWindowsDropDown_3Label     matlab.ui.control.Label
        OpenWindowsDropDown_3          matlab.ui.control.DropDown
        SelectSurfaceCargoImagePanel   matlab.ui.container.Panel
        OpenWindowsDropDown_2Label     matlab.ui.control.Label
        OpenWindowsDropDown_2          matlab.ui.control.DropDown
        orLoadFileButton_2             matlab.ui.control.Button
        ChSpinner_2Label               matlab.ui.control.Label
        ChSpinner_2                    matlab.ui.control.Spinner
        MatlabVariableDropDown_2Label  matlab.ui.control.Label
        MatlabVariableDropDown_2       matlab.ui.control.DropDown
        ViewSurfaceCargoImageButton    matlab.ui.control.Button
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Callback function
        function SelectCellFillDropDownValueChanged(app, event)
            value = app.OpenWindowsDropDown.Value;
            
        end
    end

    % Component initialization
    methods (Access = private)
        %-----added--------------
        function imagelist = updateImageList(app)
            try
                IJfilenames = cell(MIJ.getListImages);
            catch
                IJfilenames = [];
            end
            figHandles = findobj('Type', 'image');
            Mfignames = arrayfun(@(x) ['MatFig Window: ' num2str(x.Parent.Parent.Number)],figHandles,'UniformOutput',false);
            imagelist.numIJfiles = numel(IJfilenames);
            imagelist.numMATfiles = numel(Mfignames);
            imagelist.imagehandles = [IJfilenames;arrayfun(@(x) {x}, figHandles)];
            imagelist.list = [IJfilenames;Mfignames];
            if ~isempty(imagelist.list)
                imagelist.list = {'  ';imagelist.list};
            else
                imagelist.list = {'  ';'Nothing Open'};
            end
        end
        function allvars = getMatvars(app)
            allvars = arrayfun(@(x) x.name,whos,'UniformOutput',false);
            for ii = 1:numel(allvars)
                if ndims(allvars{ii})<2 %only include possible images
                    allvars{ii} = [];
                end
            end
             if ~isempty(allvars)
                allvars = ['  ';allvars];
            else
                allvars = {'  ';'Nothing Open'};
            end
        end
        %----------------------------
        % Create UIFigure and components
        function createComponents(app)
            %---- top functions -----
            %-- if nothing open then just return empty strings
            imglist = app.updateImageList;
            matvars = app.getMatvars;
            %--------------------------
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 534 279];
            app.UIFigure.Name = 'UI Figure';

            % Create SelectCellFillImagePanel
            app.SelectCellFillImagePanel = uipanel(app.UIFigure);
            app.SelectCellFillImagePanel.Title = 'Select Cell Fill Image';
            app.SelectCellFillImagePanel.Position = [25 151 485 111];

            % Create orLoadFileButton
            app.orLoadFileButton = uibutton(app.SelectCellFillImagePanel, 'push');
            app.orLoadFileButton.Position = [223 58 79 22];
            app.orLoadFileButton.Text = 'or Load File';

            % Create ChSpinnerLabel
            app.ChSpinnerLabel = uilabel(app.SelectCellFillImagePanel);
            app.ChSpinnerLabel.HorizontalAlignment = 'right';
            app.ChSpinnerLabel.Position = [223 17 25 22];
            app.ChSpinnerLabel.Text = 'Ch';

            % Create ChSpinner
            app.ChSpinner = uispinner(app.SelectCellFillImagePanel);
            app.ChSpinner.Position = [256 17 42 22];

            % Create MatlabVariableDropDownLabel
            app.MatlabVariableDropDownLabel = uilabel(app.SelectCellFillImagePanel);
            app.MatlabVariableDropDownLabel.HorizontalAlignment = 'right';
            app.MatlabVariableDropDownLabel.Position = [16 17 88 22];
            app.MatlabVariableDropDownLabel.Text = 'Matlab Variable';

            % Create MatlabVariableDropDown - select from MatVar
            app.MatlabVariableDropDown = uidropdown(app.SelectCellFillImagePanel,...
                'Items',matvars,...
                'Value',matvars{1},...
                'ValueChangedFcn',@(dd,event) cellfillselection_matvar(dd,app,imglist));
            app.MatlabVariableDropDown.Position = [114 17 87 22];

            % Create ViewCellFillImageButton
            app.ViewCellFillImageButton = uibutton(app.SelectCellFillImagePanel, 'push');
            app.ViewCellFillImageButton.Position = [336 34 123 42];
            app.ViewCellFillImageButton.Text = {'View Cell Fill'; 'Image'};

            % Create OpenWindowsDropDown_3Label
            app.OpenWindowsDropDown_2Label = uilabel(app.SelectCellFillImagePanel);
            app.OpenWindowsDropDown_2Label.HorizontalAlignment = 'right';
            app.OpenWindowsDropDown_2Label.Position = [15 54 90 22];
            app.OpenWindowsDropDown_2Label.Text = 'Open Windows';

            % Create OpenWindowsDropDown_-- select from Open Window
            app.OpenWindowsDropDown_2 = uidropdown(app.SelectCellFillImagePanel,...
                'Items',imglist.list,...
                'Value',imglist.list{1},...
                'ValueChangedFcn',@(dd,event) cellfillselection_window(dd,app,imglist));
            app.OpenWindowsDropDown_2.Position = [115 54 87 22];
            
            function cellfillselection_window(dd,app,imagelist)
                val = dd.Value;
                
                if val< imagelist.numIJfiles+1
                    assignin('base','cellfill_image', MIJ.getImage(val));
                else
                    try 
                        im = getImage(val);
                        assignin('base','cellfill_image', im);
                    catch
                        disp('woops');
                    end
                end
                app.MatlabVariableDropDown.Value = imagelist.list{1}; %set the other option for retrieving image back to blank
                app.updateImageList;                
            end
            function cellfillselection_matvar(dd,app,imagelist)
                val = dd.Value;
                app.OpenWindowsDropDown_2.Value = imagelist.list{val}; %set the other option for retrieving image back to blank
                app.updateImageList;                
            end

            % Create SelectSurfaceCargoImagePanel
            app.SelectSurfaceCargoImagePanel = uipanel(app.UIFigure);
            app.SelectSurfaceCargoImagePanel.Title = 'Select Surface Cargo Image';
            app.SelectSurfaceCargoImagePanel.Position = [25 28 485 111];

            % Create OpenWindowsDropDown_2Label
            app.OpenWindowsDropDown_3Label = uilabel(app.SelectSurfaceCargoImagePanel);
            app.OpenWindowsDropDown_3Label.HorizontalAlignment = 'right';
            app.OpenWindowsDropDown_3Label.Position = [14 58 90 22];
            app.OpenWindowsDropDown_3Label.Text = 'Open Windows';

            % Create OpenWindowsDropDown_3 _-- select from Open Window
            app.OpenWindowsDropDown_3 = uidropdown(app.SelectSurfaceCargoImagePanel,...
                'Items',imglist.list,...
                'Value',imglist.list{1});
            app.OpenWindowsDropDown_3.Position = [114 58 87 22];

            % Create orLoadFileButton_2
            app.orLoadFileButton_2 = uibutton(app.SelectSurfaceCargoImagePanel, 'push');
            app.orLoadFileButton_2.Position = [223 58 79 22];
            app.orLoadFileButton_2.Text = 'or Load File';

            % Create ChSpinner_2Label
            app.ChSpinner_2Label = uilabel(app.SelectSurfaceCargoImagePanel);
            app.ChSpinner_2Label.HorizontalAlignment = 'right';
            app.ChSpinner_2Label.Position = [223 17 25 22];
            app.ChSpinner_2Label.Text = 'Ch';

            % Create ChSpinner_2
            app.ChSpinner_2 = uispinner(app.SelectSurfaceCargoImagePanel);
            app.ChSpinner_2.Position = [256 17 42 22];

            % Create MatlabVariableDropDown_2Label
            app.MatlabVariableDropDown_2Label = uilabel(app.SelectSurfaceCargoImagePanel);
            app.MatlabVariableDropDown_2Label.HorizontalAlignment = 'right';
            app.MatlabVariableDropDown_2Label.Position = [16 17 88 22];
            app.MatlabVariableDropDown_2Label.Text = 'Matlab Variable';

            % Create MatlabVariableDropDown_2 - select from MatVars
            app.MatlabVariableDropDown_2 = uidropdown(app.SelectSurfaceCargoImagePanel,...
                'Items',matvars,...
                'Value',matvars{1});
            app.MatlabVariableDropDown_2.Position = [114 17 87 22];

            % Create ViewSurfaceCargoImageButton
            app.ViewSurfaceCargoImageButton = uibutton(app.SelectSurfaceCargoImagePanel, 'push');
            app.ViewSurfaceCargoImageButton.Position = [336 34 123 42];
            app.ViewSurfaceCargoImageButton.Text = {'View Surface Cargo'; 'Image'};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AshleyAnalysisGUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end