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
        scaleval = [];
        sliders = [];
        % vars to hide 
        rawtiff
        default_scaleval = 3;
        colors = lines(10);
        dimension = [1 1 2.5];
    end
    
    methods
        function loadimagefile(obj,filepath)
            rawtiff = loadtiff(filepath);
            obj.rawtiff = rawtiff;
            obj.numchannels = size(obj.rawtiff,4);
            for nn = 1:obj.numchannels
               obj.image{nn} = gaussf(obj.rawtiff(:,:,:,nn),[1 1 0]); 
               obj.scaleval{nn} = obj.default_scaleval; 
            end
            obj.filepath  = filepath;
            close all;
        end 
        
        function makepanel(obj)
            obj.s_figure = figure;
            obj.s_axes = obj.h_figure.CurrentAxes;
            for nn = 1:obj.numchannels
                obj.sliders{nn} = uicontrol('Parent',obj.h_figure,'Style','slider','Position',[121,50*nn,410,23],...
                    'value',0.5, 'min',0, 'max',1);
            end
        end
        
        function initialize3Dimage(obj)
            obj.h_figure = figure;
            obj.h_axes = gca;
            obj.makepatches(obj.h_axes);
        end
        function makepatches(obj,ax)
            for nn = 1:obj.numchannels
                obj.patch{nn} = obj.makepatch(ax,obj.image{nn},obj.scaleval{nn},obj.colors(nn,:));
                hold on;
            end
            pbaspect(obj.dimension);
             axis vis3d
        end
%         function sliderupdate(hObject) 
%         end
%         c = uicontrol('Parent',newfig,'Style','slider','Position',[121,150,410,23],...
%     'value',zeta, 'min',0, 'max',1);

    end
    
    methods(Static)
        function p = makepatch(ax,data,scaleval,col)
            szd = size(data);
            colormatrix = ones(szd(2),szd(1),szd(3));
            value = (max(data)+min(data))/scaleval;
            p = patch(ax, isosurface(0:szd(1)-1,0:szd(2)-1,1:szd(3),double(data),value,colormatrix.*double(data)),'parent',ax...
                ,'facelighting','phong','facecolor',col,'edgecolor','none','FaceAlpha',.5);
        end
    end
end
