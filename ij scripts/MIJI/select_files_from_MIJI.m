classdef viewer3D < handle
    properties
        h_figure;
        h_axes;
         s_figure;
        s_axes;
        filepath=[];
        numchannels = 0;
        image = [];
        patch = [];
        endpatch = [];
        scaleval = [];
        sliders = [];
        % vars to hide 
        rawtiff
        default_scaleval = 3;
        colors = [1 1 1; 0 1 0; 1 0 1; lines(10)];
        dimension = [1 1 2.5];
    end
    
    methods
        function loadimagefile(obj,filepath)
        end
        function make_selector(obj)
            fig = uifigure;
            currfiles = cell(MIJ.getListImages);
            dd = uidropdown(fig, 'Items',currfiles,...
                'Value', currfiles{1},'Callback',@obj.ch_select);
        end
        function ch_select(obj)
            % set variable name to what was entered
            
            % determine the size of the image
                % if it's multi-colored then ask which channel
            
            
            
        end
    end
