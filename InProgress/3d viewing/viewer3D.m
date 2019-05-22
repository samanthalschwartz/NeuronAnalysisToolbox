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
            rawtiff = loadtiff(filepath);
            obj.rawtiff = rawtiff;
            obj.numchannels = size(obj.rawtiff,4);
            for nn = 1:obj.numchannels
               obj.image{nn} = gaussf(obj.rawtiff(200:600,200:600,:,nn),[1 1 0]); 
               obj.scaleval{nn} = obj.default_scaleval; 
            end
            obj.filepath  = filepath;
            close all;
        end 
        function loadNDFile(obj,filepath)
            [fullimage,~]=SIM.ndFileloader(filepath);
            obj.rawtiff = fullimage;
            obj.numchannels = size(obj.rawtiff,4);
            for nn = 1:obj.numchannels
               obj.image{nn} = gaussf(obj.rawtiff(200:600,200:600,:,nn),[1 1 0]); 
               obj.scaleval{nn} = obj.default_scaleval; 
            end
            obj.filepath  = filepath;   
        end

        
        function loadSIM(obj,filepath)
            ss = load(filepath);
            obj.numchannels = 3;
            obj.image{1} = dip_image(ss.obj.abeta.image);
            obj.image{2} = dip_image(ss.obj.ch1.image);
            obj.image{3} = dip_image(ss.obj.ch2.image);
            for nn = 1:obj.numchannels
                obj.scaleval{nn} = obj.default_scaleval;
            end
        end
        
        function makepanel(obj)
            obj.s_figure = figure;
            obj.s_axes = obj.s_figure.CurrentAxes;
            for nn = 1:obj.numchannels
                obj.sliders{nn} = uicontrol('Parent',obj.s_figure,'Style','slider','Position',[121,50*nn,410,23],...
                    'value',3, 'min',1, 'max',20,'Callback',@obj.update);
            end
        end
   
        function update(obj,src,eventData)
            for nn = 1:obj.numchannels
                currval = obj.scaleval{nn};
                if ~isequal(currval,get(obj.sliders{nn},'Value'))
                    ch2change = nn;
                    obj.scaleval{nn} = get(obj.sliders{nn},'Value');
                else
                    continue
                end
            end
            obj.updatepatch(ch2change);
        end
        
        function slider_val = slidercallback(hObject)
            slider_val=floor(get(hObject,'Value'));            
        end
        
        function initialize3Dimage(obj)
            obj.h_figure = figure;
            obj.h_axes = gca;
            obj.makepatches();
            axis 'vis3d'
            set(gca,'color','k');
        end
        function makepatches(obj)
            cla(obj.h_axes);
            for nn = 1:obj.numchannels
                [obj.patch{nn},obj.endpatch{nn}] = obj.makepatch(obj.h_axes,obj.image{nn},obj.scaleval{nn},obj.colors(nn,:));
                hold on;
            end
            pbaspect(obj.dimension);
        end
        function updatepatch(obj,ch)
            cla(obj.h_axes);
            nn = ch;
            [obj.patch{nn},obj.endpatch{nn}] = obj.makepatch(obj.h_axes,obj.image{nn},obj.scaleval{nn},obj.colors(nn,:));
            hold on;
            pbaspect(obj.dimension);
        end
        
%         function sliderupdate(hObject) 
%         end
%         c = uicontrol('Parent',newfig,'Style','slider','Position',[121,150,410,23],...
%     'value',zeta, 'min',0, 'max',1);

    end
    
    methods(Static)
        function [p, pe]= makepatch(ax,data,scaleval,col)
            szd = size(data);
            colormatrix = ones(szd(2),szd(1),szd(3));
            value = (max(data)+min(data))/scaleval;
            p = patch(ax, isosurface(0:szd(1)-1,0:szd(2)-1,1:szd(3),double(data),value,colormatrix.*double(data)),'parent',ax...
                ,'facelighting','phong','facecolor',col,'edgecolor','none','FaceAlpha',.5);
            pe = patch(isocaps(0:szd(1)-1,0:szd(2)-1,1:szd(3),double(data),value),'parent',ax,...
                'facecolor',col,'edgecolor','none','FaceAlpha',.5,'facelighting','phong') ;
        end
    end
end
