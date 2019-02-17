classdef CU_FRAP < handle
    properties
        %         fp.ROIs = {[328, 186],[209, 241],[211, 302],[332, 335]};
        % fp.ROIs = {[305, 205],[49, 38],[330, 317],[407, 294]};
        filestr = {'*t*_w0000.tif','*t*_w0001.tif','*.tif'};
        colors = lines(20);
        image = [];
        framerate = []; % frame rate in seconds/frame
%         background = []; % this is background/pixel
        backgroundROI = [];
%         control = [];
        controlROI = [];
        FRAP_ROI = [];
        datadir = [];
        ROIs = []; % list of ROIs to be analyzed
        ROIrect = [];
        %cell array of FRAP_ROI objects
        bleachframes = []; % m x 1 array of frames where photobleaching occured (m<n). First frame = 1. Example: Bleaching from frames 10 - 15 would be the array [10:15];
        inc = 20; %size of bleach square
    end
    
    methods
        function fitRecoveryCurve_1component(obj)

        end
        function load(obj,datadir,uniquestr)
            if nargin<3
                uniquestr = obj.filestr{3};
                if nargin<2
                    datadir = [];
                end
            end
            [obj.image,datadir] = obj.loadtiffSeries(datadir,uniquestr,'sum');
            obj.datadir = datadir;
        end
        function h = showROIs(obj)
           h = dipshow(obj.image,'log');
            diptruesize(h,100);
            for pp = 1:numel(obj.ROIs)
                roi2rectangle(obj,obj.ROIs{pp},obj.colors(pp,:),gca);
            end            
        end
        function adjustROIs(obj)
            for pp = 1:numel(obj.ROIs)
                 h = dipshow(obj.image,'log');
                 diptruesize(h,100);
                 set(h,'Units','pixels')
                 hpos = get(h,'Position');
                 
                 test = rectangle('Position',[obj.ROIs{pp}(1)-obj.inc,obj.ROIs{pp}(2)-obj.inc,obj.inc,obj.inc],...
                     'EdgeColor',obj.colors(pp,:),...
                     'LineWidth',1);
                 % could use this build in matlab function instead of
                 % draggable.mat
                 %                  waitforbuttonpress
                 %                  inputrect =test.Position;
                 %                  inputrect(2) = size(obj.image,2)-inputrect(2)-inputrect(4);
                 %                  obj.ROIs{pp} = dragrect(inputrect);
                 tt = warndlg('Test'); 
                 set(tt,'Units','pixels')
                 ttpos = get(tt,'Position');
                 ttpos(1) = hpos(1)+hpos(3);
                 set(tt,'Position',ttpos);
                 waitforbuttonpress
                 %                  currkey = 0;
                draggable(test);
                waitfor(tt);
%                 waitforbuttonpress
                obj.ROIs{pp} = test.Position(1:2);
%                 currkey= 0;
%                 w = waitforbuttonpress
%                 while currkey~=1
% %                      pause; % wait for a keypress
%                      currkey=get(gcf,'CurrentKey');
%                      if strcmp(currkey, 'enter') % You also want to use strcmp here.
%                          obj.ROIs{pp} = test.Position(1:2);
%                          currkey=1; % Error was here; the "==" should be "="
%                      else
%                          currkey=0; % Error was here; the "==" should be "="
%                      end
%                 end
                 close(h);
            end
        end
        
        function getControlROI(obj)
            h = dipshow(sum(obj.image,[],3),'log');
            diptruesize(h,100);
            for pp = 1:numel(obj.ROIs)
                roi2rectangle(obj,obj.ROIs{pp},'g',gca);
            end
            [~,C] = dipcrop(h);
            obj.controlROI.image = obj.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
            obj.controlROI.coords = C;
            close(h);
        end
        function getBackgroundROI(obj)
            h = dipshow(sum(obj.image,[],3),'log');
            diptruesize(h,100);
            [B,C] = dipcrop(h);
            obj.backgroundROI.image = obj.image(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:);
            obj.backgroundROI.coords = C;
%             obj.background = sum(B)/numel(single(obj.image));
            close(h);  
        end
        function rect = roi2rectangle(obj,roi,cols,ax)
            rect = rectangle('Position',[roi(1)-floor(obj.inc/2),roi(2)-floor(obj.inc/2),obj.inc,obj.inc],...
                    'EdgeColor',cols,...
                    'LineWidth',1);
        end
        function im = roi2image(obj,roi)
            im = obj.image(roi(1)-floor(obj.inc/2):roi(1)+floor(obj.inc/2),roi(2)-floor(obj.inc/2):roi(2)+floor(obj.inc/2),:);
        end
        function h = setROIs(obj)
            for rr = 1:numel(obj.ROIs)
                img_rr = obj.roi2image(obj.ROIs{rr}); 
                obj.FRAP_ROI{rr} = CU_FRAP_ROI();
                obj.FRAP_ROI{rr}.Intensity = img_rr;
                obj.FRAP_ROI{rr}.bleachframes = [obj.bleachframes];
                obj.FRAP_ROI{rr}.Background = obj.backgroundROI.image; %per pixel background
                obj.FRAP_ROI{rr}.Control = obj.controlROI.image;
                time = (1:size(img_rr,3))*obj.framerate;
                obj.FRAP_ROI{rr}.time = time;
            end
            h = obj.showROIs;
        end
        function [h2, h1] = plotRecoveryCurves(obj)
            if isempty(obj.FRAP_ROI)
                h = obj.setROIs;
                close(h);
            end
            h1 = obj.showROIs();
            h2 = figure(); ah2 = gca(); hold on;
           for rr = 1:numel(obj.ROIs)
               obj.FRAP_ROI{rr}.makeRecoveryCurve();
               x = obj.FRAP_ROI{rr}.RecoveryCurve.x;
               y = obj.FRAP_ROI{rr}.RecoveryCurve.y;
               plot(ah2,x,y,'Color',obj.colors(rr,:));  hold on;
               plot(ah2,obj.FRAP_ROI{rr}.RecoveryCurve.fitx,obj.FRAP_ROI{rr}.RecoveryCurve.fity,'--','Color',obj.colors(rr,:)); 
           end
           xlabel('Time (seconds)');
           ylabel('Relative Fluorescence');   
        end
    end
    methods (Static)
        function [image,datadir] = loadtiffSeries(datadir,uniquefilestr,zproject)
            if nargin<3
                zproject = 'sum';
            end
            if nargin < 1 || isempty(datadir)
                datadir = uigetdir('.','Select the folder containing the .tif series for Channel1 and Channel2');
            end
            files = dir(fullfile(datadir,uniquefilestr));
            frame1 = loadtiff(fullfile(datadir,files(1).name));
            image = dip_image(zeros(size(frame1,1),size(frame1,2),numel(files)));
            switch zproject
                case 'sum'
                    for ff = 1:numel(files)
                        frame = loadtiff(fullfile(datadir,files(ff).name));
                        image(:,:,ff-1) = sum(dip_image(frame),[],3);
                    end
                case 'max'
                    for ff = 1:numel(files)
                        frame = loadtiff(fullfile(datadir,files(ff).name));
                        image(:,:,ff-1) = max(dip_image(frame),[],3);
                    end
            end
        end
        function ROImean = calcROImean(in_ROI)
            im = dip_image(in_ROI);
            ROImean = single(squeeze(mean(im,[],[1 2])));            
        end
        function ROIsum = calcROIsum(in_ROI)
            im = dip_image(in_ROI);
            ROIsum = single(squeeze(sum(im,[],[1 2])));  
        end 
    end
end
